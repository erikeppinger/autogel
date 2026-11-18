"""Real-world test set from open-access papers: original, uncropped gel and blot images.

Many journals require the original, uncropped images of every gel and blot (PLOS since 2019,
eLife, MDPI, Springer Nature, Frontiers, ...). This tool finds CC-BY papers that publish them,
lists their raw-image supplements (and data files with the authors' quantification) in a manifest,
fetches them into a local cache and extracts the pictures.

Sources:
  plos       PLOS search API (10 requests/min), one download per supplementary file
  europepmc  Europe PMC REST API: open-access CC-BY articles of many publishers; supplements come
             as one ZIP per article (capped by --max-zip-mb), only the matching files are kept

    python tools/oa_testset.py search  --source europepmc [--kinds western rna dna] [--per-kind 20]
    python tools/oa_testset.py fetch   [--max-mb 40] [--max-zip-mb 200]
    python tools/oa_testset.py extract

The manifest (DOIs, licences, file names, captions) is small and can be shared; the files stay in
the cache (default ``benchmark_data/oa``, not in git) and can be rebuilt any time. Only CC-BY
articles are listed; retracted articles are skipped.
"""
from __future__ import annotations

import argparse
import io
import json
import os
import re
import shutil
import sys
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
import zipfile

HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT_MANIFEST = os.path.join(HERE, "oa_testset.json")
DEFAULT_CACHE = os.path.join(os.path.dirname(HERE), "benchmark_data", "oa")
UA = {"User-Agent": "autogel-testset (open-source gel analysis benchmark)"}
XLINK = "{http://www.w3.org/1999/xlink}href"
RAW_CAPTION = re.compile(r"raw (blot |gel |western )?image|uncropped|original (unedited |uncropped )?"
                         r"(blot|gel|image|western)|full[- ]length (blot|gel|western)|unprocessed (blot|gel|image)|"
                         r"full (blot|gel) image|whole (blot|membrane)", re.I)
DATA_CAPTION = re.compile(r"\b(data|dataset|values|quantification|quantified|densitometr|underlying)\b", re.I)
DATA_EXT = (".xlsx", ".xls", ".csv")
IMG_EXT = (".tif", ".tiff", ".png", ".jpg", ".jpeg", ".jp2", ".bmp")
RETRACTED = re.compile(r"\s*(retracted|retraction|expression of concern|withdrawn)", re.I)

# ---------------------------------------------------------------------------- PLOS
PLOS_API = "https://api.plos.org/search"
PLOS_JOURNAL = {"PLOS ONE": "plosone", "PLOS Biology": "plosbiology", "PLOS Genetics": "plosgenetics",
                "PLOS Pathogens": "plospathogens", "PLOS Computational Biology": "ploscompbiol",
                "PLOS Medicine": "plosmedicine", "PLOS Neglected Tropical Diseases": "plosntds"}
PLOS_QUERIES = {
    "western": 'everything:"western blot" AND everything:"densitometry"',
    "rna": 'everything:"28S" AND everything:"18S" AND everything:"RNA integrity"',
    "dna": 'everything:"agarose gel" AND everything:"ladder" AND everything:"PCR"',
}
# ---------------------------------------------------------------------------- Europe PMC
EPMC = "https://www.ebi.ac.uk/europepmc/webservices/rest"
EPMC_QUERIES = {
    "western": '"uncropped" AND "western blot"',
    "rna": '("28S" AND "18S") AND ("RNA integrity" OR "RNA quality") AND ("uncropped" OR "raw image" OR '
           '"original gel" OR "full gel")',
    "dna": '"agarose gel" AND ("uncropped gel" OR "raw gel image" OR "original gel image" OR "uncropped image")',
}
EPMC_FILTER = ' AND OPEN_ACCESS:y AND LICENSE:"cc by" AND HAS_SUPPL:y AND NOT SRC:PPR AND PUB_YEAR:[2019 TO 2026]'

_last = {"plos": 0.0, "epmc": 0.0}


def _get(url: str, timeout: float = 90, pace: str | None = None, gap: float = 1.0) -> bytes:
    if pace:
        wait = gap - (time.time() - _last[pace])
        if wait > 0:
            time.sleep(wait)
        _last[pace] = time.time()
    with urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=timeout) as r:
        return r.read()


def _caption(s) -> str:
    c = s.find("caption")
    return " ".join("".join(c.itertext()).split())[:300] if c is not None else ""


def _cc_by(lic: str | None) -> bool:
    lic = (lic or "").lower().replace("https://", "http://").strip()
    return lic.startswith("http://creativecommons.org/licenses/by/") or lic in ("cc by", "cc-by")


# ---------------------------------------------------------------------------- search
def plos_search(kind: str, n: int, seen: set) -> list[dict]:
    out, start = [], 0
    while len(out) < n:
        q = urllib.parse.urlencode({"q": f'({PLOS_QUERIES[kind]}) AND everything:"raw images"',
                                    "fl": "id,title_display,journal,publication_date,article_type",
                                    "rows": 20, "start": start, "wt": "json",
                                    "fq": "doc_type:full AND publication_date:[2019-07-01T00:00:00Z TO NOW]"})
        res = json.loads(_get(f"{PLOS_API}?{q}", pace="plos", gap=6.5))["response"]
        if not res["docs"]:
            break
        for d in res["docs"]:
            if len(out) >= n or d["id"] in seen or d.get("article_type") != "Research Article" \
                    or RETRACTED.match(d.get("title_display", "")):
                continue
            path = PLOS_JOURNAL.get(d.get("journal", ""), "plosone")
            try:
                root = ET.fromstring(_get(f"https://journals.plos.org/{path}/article/file?id={d['id']}&type=manuscript",
                                          pace="epmc"))
            except Exception as e:
                print(f"  skip {d['id']}: {e}", file=sys.stderr)
                continue
            lic = root.find(".//license")
            lic = lic.get(XLINK) if lic is not None else None
            files = []
            for s in root.iter("supplementary-material"):
                href = (s.get(XLINK) or "").replace("info:doi/", "")
                cap, label, mime = _caption(s), (s.findtext("label") or "").strip(), s.get("mimetype") or ""
                table = any(t in mime for t in ("sheet", "excel", "csv"))
                kind_f = ("raw" if RAW_CAPTION.search(cap) or RAW_CAPTION.search(label)
                          else "data" if DATA_CAPTION.search(cap) and table else None)
                if kind_f:
                    files.append({"kind": kind_f, "name": href.split("/")[-1], "caption": cap, "mimetype": mime,
                                  "url": f"https://journals.plos.org/{path}/article/file?id={href}&type=supplementary"})
            if not _cc_by(lic) or not any(f["kind"] == "raw" for f in files):
                continue
            seen.add(d["id"])
            out.append({"source": "plos", "kind": kind, "id": d["id"].split("/")[-1], "doi": d["id"],
                        "title": d.get("title_display", ""), "journal": d.get("journal", ""),
                        "date": d.get("publication_date", "")[:10], "license": lic, "files": files})
            _say(out[-1])
        start += len(res["docs"])
        if start >= res["numFound"]:
            break
    return out


def epmc_search(kind: str, n: int, seen: set) -> list[dict]:
    out, cursor = [], "*"
    while len(out) < n:
        q = urllib.parse.urlencode({"query": EPMC_QUERIES[kind] + EPMC_FILTER, "format": "json", "pageSize": 25,
                                    "resultType": "core", "cursorMark": cursor, "sort": "P_PDATE_D desc"})
        res = json.loads(_get(f"{EPMC}/search?{q}", pace="epmc"))
        docs = res.get("resultList", {}).get("result", [])
        if not docs:
            break
        for d in docs:
            doi, pmcid = (d.get("doi") or "").lower(), d.get("pmcid")
            pubtypes = " ".join(d.get("pubTypeList", {}).get("pubType", [])).lower()
            if len(out) >= n or not pmcid or doi in seen or doi.startswith("10.1371/") \
                    or RETRACTED.match(d.get("title", "")) or "retract" in pubtypes or not _cc_by(d.get("license")):
                continue
            try:
                root = ET.fromstring(_get(f"{EPMC}/{pmcid}/fullTextXML", pace="epmc"))
            except Exception as e:
                print(f"  skip {pmcid}: {e}", file=sys.stderr)
                continue
            files, names = [], set()
            for s in root.iter("supplementary-material"):
                media = [m.get(XLINK) for m in s.iter("media")] or [s.get(XLINK)]
                cap, label = _caption(s), (s.findtext("label") or "").strip()
                for m in media:
                    if not m or m in names:
                        continue
                    ext = os.path.splitext(m)[1].lower()
                    kind_f = "raw" if (RAW_CAPTION.search(cap) or RAW_CAPTION.search(label)) and ext not in DATA_EXT \
                        else "data" if DATA_CAPTION.search(cap) and ext in DATA_EXT else None
                    if kind_f:
                        names.add(m)
                        files.append({"kind": kind_f, "name": m, "caption": cap})
            if not any(f["kind"] == "raw" for f in files):
                continue
            seen.add(doi)
            out.append({"source": "europepmc", "kind": kind, "id": pmcid, "doi": doi, "title": d.get("title", ""),
                        "journal": d.get("journalInfo", {}).get("journal", {}).get("title", ""),
                        "date": d.get("firstPublicationDate", ""), "license": d.get("license"), "files": files})
            _say(out[-1])
        nxt = res.get("nextCursorMark")
        if not nxt or nxt == cursor:
            break
        cursor = nxt
    return out


def _say(a: dict):
    raw = sum(f["kind"] == "raw" for f in a["files"])
    print(f"{a['source']:9s} {a['kind']:8s} {a['id']:14s} raw {raw:2d} data {len(a['files']) - raw:2d}  "
          f"{a['journal'][:22]:22s} {a['title'][:55]}", flush=True)


def cmd_search(a):
    old = _load(a.manifest) if os.path.exists(a.manifest) else []
    seen = {x["doi"].lower() for x in old}
    new = []
    for kind in a.kinds:
        fn = plos_search if a.source == "plos" else epmc_search
        new += fn(kind, a.per_kind, seen)
    arts = old + new
    with open(a.manifest, "w", encoding="utf-8") as f:
        json.dump({"about": "CC-BY open-access articles with raw gel/blot images (autogel test set)",
                   "updated": time.strftime("%Y-%m-%d"), "articles": arts}, f, indent=1, ensure_ascii=False)
    print(f"{len(new)} new, {len(arts)} articles in {a.manifest}")


# ---------------------------------------------------------------------------- fetch / extract
def _load(path: str) -> list[dict]:
    with open(path, encoding="utf-8") as f:
        return json.load(f)["articles"]


def _save(path: str, data: bytes):
    with open(path, "wb") as f:
        f.write(data)


def _folder(cache: str, art: dict) -> str:
    return os.path.join(cache, art["source"], art["id"])


MIME_EXT = {"application/pdf": ".pdf", "image/tiff": ".tif", "image/png": ".png", "image/jpeg": ".jpg",
            "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet": ".xlsx",
            "application/vnd.ms-excel": ".xls", "text/csv": ".csv"}


def _local_name(f: dict) -> str:
    """File name in the cache (PLOS file ids have no extension: taken from the mime type)."""
    name = os.path.basename(f["name"])
    return name if os.path.splitext(name)[1] else name + MIME_EXT.get(f.get("mimetype", ""), "")


def _download_capped(url: str, limit: int) -> bytes | None:
    buf = io.BytesIO()
    with urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=300) as r:
        while chunk := r.read(1 << 20):
            buf.write(chunk)
            if buf.tell() > limit:
                return None
    return buf.getvalue()


def cmd_fetch(a):
    total = 0
    for art in _load(a.manifest):
        d = _folder(a.cache, art)
        os.makedirs(d, exist_ok=True)
        wanted = [f for f in art["files"] if not os.path.exists(os.path.join(d, _local_name(f)))]
        if not wanted:
            continue
        if art["source"] == "plos":
            for f in wanted:
                data = _download_capped(f["url"], a.max_mb * 1_000_000)
                if data is None:
                    print(f"  skip {f['name']}: larger than {a.max_mb} MB")
                    continue
                _save(os.path.join(d, _local_name(f)), data)
                total += len(data)
                time.sleep(1.0)
        else:
            _last["epmc"] = time.time()
            data = _download_capped(f"{EPMC}/{art['id']}/supplementaryFiles", a.max_zip_mb * 1_000_000)
            if data is None:
                print(f"  skip {art['id']}: supplements larger than {a.max_zip_mb} MB")
                continue
            total += len(data)
            with zipfile.ZipFile(io.BytesIO(data)) as z:
                members = {os.path.basename(m): m for m in z.namelist()}
                for f in wanted:
                    m = members.get(os.path.basename(f["name"]))
                    if m:
                        _save(os.path.join(d, _local_name(f)), z.read(m))
            time.sleep(1.0)
        print(f"{art['source']:9s} {art['id']:14s} {len(os.listdir(d))} files", flush=True)
    print(f"downloaded {total / 1e6:.1f} MB into {a.cache}")


def _unzip_images(path: str, out: str) -> int:
    n = 0
    with zipfile.ZipFile(path) as z:
        for m in z.namelist():
            if "__MACOSX" in m or m.endswith("/"):
                continue
            ext = os.path.splitext(m)[1].lower()
            if ext in IMG_EXT or ext == ".pdf":
                dst = os.path.join(out, os.path.basename(path)[:-4] + "__" + os.path.basename(m))
                if not os.path.exists(dst):
                    _save(dst, z.read(m))
                    n += 1
    return n


def _pdf_images(path: str, out: str) -> int:
    from pypdf import PdfReader
    n = 0
    for p, page in enumerate(PdfReader(path).pages, 1):
        for k, img in enumerate(page.images, 1):
            ext = os.path.splitext(img.name)[1] or ".png"
            dst = os.path.join(out, f"{os.path.basename(path)[:-4]}_p{p}_{k}{ext}")
            if not os.path.exists(dst):
                _save(dst, img.data)
                n += 1
    return n


def cmd_extract(a):
    """Pictures out of PDFs and (nested) ZIPs -> <article>/images/; plain images are copied there."""
    n = 0
    for art in _load(a.manifest):
        d = _folder(a.cache, art)
        if not os.path.isdir(d):
            continue
        out = os.path.join(d, "images")
        os.makedirs(out, exist_ok=True)
        for fn in sorted(os.listdir(d)):
            p = os.path.join(d, fn)
            if os.path.isdir(p):
                continue
            ext = os.path.splitext(fn)[1].lower()
            try:
                if ext == ".zip":
                    n += _unzip_images(p, out)
                elif ext in IMG_EXT and not os.path.exists(os.path.join(out, fn)):
                    shutil.copy(p, out)
                    n += 1
            except Exception as e:
                print(f"  {fn}: {type(e).__name__}: {e}")
        # PDFs, also those that came out of ZIPs
        for src in [os.path.join(d, f) for f in os.listdir(d)] + [os.path.join(out, f) for f in os.listdir(out)]:
            if src.lower().endswith(".pdf"):
                try:
                    n += _pdf_images(src, out)
                except Exception as e:
                    print(f"  {os.path.basename(src)}: {type(e).__name__}: {e}")
    print(f"{n} images extracted")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("search")
    s.add_argument("--source", choices=["plos", "europepmc"], default="europepmc")
    s.add_argument("--kinds", nargs="+", default=["western", "rna", "dna"], choices=["western", "rna", "dna"])
    s.add_argument("--per-kind", type=int, default=20)
    f = sub.add_parser("fetch")
    f.add_argument("--max-mb", type=int, default=40, help="per file (PLOS)")
    f.add_argument("--max-zip-mb", type=int, default=200, help="per article ZIP (Europe PMC)")
    for p in (s, f, sub.add_parser("extract")):
        p.add_argument("--manifest", default=DEFAULT_MANIFEST)
        p.add_argument("--cache", default=DEFAULT_CACHE)
    a = ap.parse_args(argv)
    {"search": cmd_search, "fetch": cmd_fetch, "extract": cmd_extract}[a.cmd](a)


if __name__ == "__main__":
    main()
