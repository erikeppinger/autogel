#!/bin/bash
# autogel installer for macOS - double-click this file (it opens in Terminal).
#
#  1. finds Python >= 3.10 (if there is none: opens the python.org download page)
#  2. creates a private environment in ~/Library/Application Support/autogel and installs
#     autogel + its viewer (napari)
#  3. creates the app "autogel" in ~/Applications (Launchpad / Spotlight; drag it to the Dock)
#  4. writes an uninstaller next to the environment
# No administrator rights, git or command-line knowledge needed.
#
# Options (for IT / testing):  --prefix <folder>   --no-launcher   --yes (no questions, no pause)
#                              --with-ai (also the optional AI second opinion, onnxruntime ~17 MB)
set -u

SRC="$(cd "$(dirname "$0")" && pwd)"
PREFIX="$HOME/Library/Application Support/autogel"
APPDIR="$HOME/Applications"
LAUNCHER=1
YES=0
EXTRAS="gui"
while [ $# -gt 0 ]; do
    case "$1" in
        --prefix) PREFIX="$2"; shift ;;
        --no-launcher) LAUNCHER=0 ;;
        --yes) YES=1 ;;
        --with-ai) EXTRAS="gui,ai" ;;
        *) echo "unknown option: $1"; exit 2 ;;
    esac
    shift
done

say()  { printf '%s\n' "$*"; }
fail() { printf '\n*** %s\n' "$*"; [ "$YES" = 1 ] || read -r -p "Press Return to close." _; exit 1; }

say "=== autogel installer ==="
say "Source:  $SRC"
say "Target:  $PREFIX"
[ -f "$SRC/pyproject.toml" ] || fail "pyproject.toml not found next to this installer. Unzip the whole autogel folder first."

# ------------------------------------------------------------------ 1. Python
ok_python() { "$1" -c 'import sys; sys.exit(0 if sys.version_info >= (3, 10) else 1)' >/dev/null 2>&1; }
PY=""
# 3.12 first: the viewer's binary packages (Qt) are most reliably available for it
candidates=(python3.12 python3.13 python3.11 python3.10 python3.14
            /Library/Frameworks/Python.framework/Versions/3.12/bin/python3
            /Library/Frameworks/Python.framework/Versions/3.1[0-9]/bin/python3
            /opt/homebrew/bin/python3 /usr/local/bin/python3)
# Apple's /usr/bin/python3 is only a stub that pops up a "install developer tools" dialog when
# the Command Line Tools are missing - only ask it when they are installed.
if [ "$(uname)" != "Darwin" ] || xcode-select -p >/dev/null 2>&1; then candidates+=(python3); fi
for c in "${candidates[@]}"; do
    p="$(command -v "$c" 2>/dev/null || true)"
    [ -n "$p" ] && [ -x "$p" ] && ok_python "$p" && { PY="$p"; break; }
done
if [ -z "$PY" ]; then
    say ""
    say "No Python 3.10 or newer found."
    say "Install Python 3.12 from python.org (the 'macOS 64-bit universal2 installer'),"
    say "then double-click this installer again."
    [ "$(uname)" = "Darwin" ] && [ "$YES" = 0 ] && open "https://www.python.org/downloads/macos/"
    fail "Python missing."
fi
say "Python:  $PY ($("$PY" -c 'import platform; print(platform.python_version())'))"

# ------------------------------------------------------------------ 2. environment + autogel
VENV="$PREFIX/venv"
mkdir -p "$PREFIX" || fail "Cannot create $PREFIX"
if [ ! -x "$VENV/bin/python" ]; then
    say "Creating environment ..."
    "$PY" -m venv "$VENV" || { rm -rf "$VENV"; fail "Could not create the environment (Linux: install the package python3-venv)."; }
fi
VPY="$VENV/bin/python"
say "Installing autogel and its viewer (napari) - this downloads ~300 MB and takes a few minutes ..."
"$VPY" -m pip install --upgrade pip --quiet --disable-pip-version-check
"$VPY" -m pip install --quiet --disable-pip-version-check "$SRC[$EXTRAS]" \
    || fail "Installation failed (see messages above). Check the internet connection / proxy."
"$VPY" -c "import autogel, napari; print('autogel', autogel.__version__, '| napari', napari.__version__)" \
    || fail "autogel was installed but cannot be imported."

# ------------------------------------------------------------------ 3. app launcher
if [ "$LAUNCHER" = 1 ] && [ "$(uname)" = "Darwin" ]; then
    APP="$APPDIR/autogel.app"
    rm -rf "$APP"
    mkdir -p "$APP/Contents/MacOS" || fail "Cannot create $APP"
    cat > "$APP/Contents/Info.plist" <<'PLIST'
<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE plist PUBLIC "-//Apple//DTD PLIST 1.0//EN" "http://www.apple.com/DTDs/PropertyList-1.0.dtd">
<plist version="1.0">
<dict>
    <key>CFBundleName</key><string>autogel</string>
    <key>CFBundleDisplayName</key><string>autogel</string>
    <key>CFBundleIdentifier</key><string>io.github.erikeppinger.autogel</string>
    <key>CFBundleExecutable</key><string>autogel</string>
    <key>CFBundlePackageType</key><string>APPL</string>
    <key>CFBundleShortVersionString</key><string>1.0</string>
    <key>NSHighResolutionCapable</key><true/>
</dict>
</plist>
PLIST
    # the launcher only starts the program in the private environment
    printf '#!/bin/bash\ncd "$HOME"\nexec "%s" -m autogel.cli gui "$@"\n' "$VPY" > "$APP/Contents/MacOS/autogel"
    chmod +x "$APP/Contents/MacOS/autogel"
    say "App created: $APP (find it in Launchpad / Spotlight, or drag it to the Dock)."
fi

# ------------------------------------------------------------------ 4. uninstaller
UN="$PREFIX/Uninstall autogel.command"
cat > "$UN" <<UNEOF
#!/bin/bash
echo "Removing autogel ..."
rm -rf "$APPDIR/autogel.app"
rm -rf "$PREFIX"
echo "Done. Your own ladders in ~/.autogel were kept."
UNEOF
chmod +x "$UN"

say ""
say "autogel is installed."
say "Start it with the app 'autogel' (Launchpad / Spotlight), or run: \"$VPY\" -m autogel.cli gui"
say "To update: download the new version and run install_mac.command again."
say "To remove: double-click \"$UN\""
[ "$YES" = 1 ] || read -r -p "Press Return to close this window." _
exit 0
