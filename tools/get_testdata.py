"""Rebuild autogel's test data on any computer with one command (nothing is stored in git).

    python tools/get_testdata.py                 # everything (about 1.5 GB download)
    python tools/get_testdata.py --only zenodo   # hand-labelled gels + Western-blot sets
    python tools/get_testdata.py --only oa       # open-access papers (tools/oa_testset.json)
    python tools/get_testdata.py --record        # write checksums of the local open-access files
    python tools/get_testdata.py --verify        # only check what is there

Zenodo files are checked against Zenodo's own MD5 checksums; open-access supplements against the
SHA-256 values in tools/testdata_checksums.json (a changed file is reported, not used silently).
Everything goes to benchmark_data/ (ignored by git). Sources and licences:
  * GelGenie hand-labelled gels, Zenodo 10.5281/zenodo.14641949, CC-BY-4.0
  * Western blots + the authors' densitometry, Zenodo 10.5281/zenodo.18338155, CC-BY-4.0
  * Western-blot figure source data, Zenodo 10.5281/zenodo.14291079, CC-BY-4.0
  * CC-BY open-access papers (PLOS, Europe PMC) listed in tools/oa_testset.json
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
import urllib.request
import zipfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(os.path.dirname(HERE), "benchmark_data")
CHECKSUMS = os.path.join(HERE, "testdata_checksums.json")
UA = {"User-Agent": "autogel-testset (open-source gel analysis benchmark)"}

# record -> (folder in benchmark_data, files to fetch or None for all)
ZENODO = {
    "14641949": ("zenodo_14641949", ["readme.txt", "quantitation_ladder_gels.zip", "nathan_gels.zip",
                                     "external_gels.zip", "stella_gels_for_finetuning.zip", "matthew_gels_2.zip"]),
    "18338155": ("wb_18338155", ["Western blots.zip", "Raw data densitometry Westernblots.xlsx"]),
    "14291079": ("wb_14291079", None),
}
UNPACK = {"14291079": False}          # print-resolution figure pages: 4.9 GB unpacked, kept zipped


def _md5(path: str) -> str:
    h = hashlib.md5()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def _sha256(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def zenodo(verify_only: bool) -> bool:
    ok = True
    for rec, (folder, want) in ZENODO.items():
        meta = json.load(urllib.request.urlopen(urllib.request.Request(
            f"https://zenodo.org/api/records/{rec}", headers=UA), timeout=60))
        d = os.path.join(ROOT, folder)
        os.makedirs(d, exist_ok=True)
        for f in meta["files"]:
            if want and f["key"] not in want:
                continue
            dst = os.path.join(d, f["key"])
            md5 = f.get("checksum", "").replace("md5:", "")
            if os.path.exists(dst) and os.path.getsize(dst) == f["size"] and _md5(dst) == md5:
                print(f"ok        {folder}/{f['key']}")
            elif verify_only:
                print(f"MISSING   {folder}/{f['key']}")
                ok = False
                continue
            else:
                print(f"download  {folder}/{f['key']} ({f['size'] / 1e6:.0f} MB)", flush=True)
                urllib.request.urlretrieve(f["links"]["self"], dst)
                if _md5(dst) != md5:
                    print(f"CHECKSUM MISMATCH {dst}")
                    ok = False
                    continue
            target = os.path.splitext(dst)[0]
            if dst.endswith(".zip") and UNPACK.get(rec, True) and not os.path.isdir(target) and not verify_only:
                with zipfile.ZipFile(dst) as z:
                    z.extractall(target)
                print(f"unpacked  {folder}/{os.path.basename(target)}")
    return ok


def oa_files() -> list[str]:
    out = []
    base = os.path.join(ROOT, "oa")
    for dp, dns, fs in os.walk(base):
        dns[:] = [x for x in dns if x != "images"]          # extracted pictures are derived
        out += [os.path.relpath(os.path.join(dp, f), ROOT).replace("\\", "/") for f in fs]
    return sorted(out)


def oa(verify_only: bool) -> bool:
    if not verify_only:
        sys.path.insert(0, HERE)
        import oa_testset
        oa_testset.main(["fetch"])
        oa_testset.main(["extract"])
    if not os.path.exists(CHECKSUMS):
        print("no checksum file yet (run with --record once)")
        return True
    with open(CHECKSUMS, encoding="utf-8") as f:
        known = json.load(f)["sha256"]
    ok = True
    for rel, h in known.items():
        p = os.path.join(ROOT, rel)
        if not os.path.exists(p):
            print(f"MISSING   {rel}")
            ok = False
        elif _sha256(p) != h:
            print(f"CHANGED   {rel} (the publisher's file differs from the recorded one)")
            ok = False
    print(f"open-access files: {len(known)} recorded, {'all present and unchanged' if ok else 'see above'}")
    return ok


def record():
    files = oa_files()
    data = {"about": "SHA-256 of the open-access test files in benchmark_data/ (see get_testdata.py)",
            "sha256": {rel: _sha256(os.path.join(ROOT, rel)) for rel in files}}
    with open(CHECKSUMS, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=1, sort_keys=True)
    print(f"recorded {len(files)} files -> {CHECKSUMS}")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--only", choices=["zenodo", "oa"])
    ap.add_argument("--verify", action="store_true", help="check only, download nothing")
    ap.add_argument("--record", action="store_true", help="write checksums of the local open-access files")
    a = ap.parse_args(argv)
    if a.record:
        record()
        return
    ok = True
    if a.only in (None, "zenodo"):
        ok &= zenodo(a.verify)
    if a.only in (None, "oa"):
        ok &= oa(a.verify)
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
