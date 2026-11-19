# autogel installer for Windows (started by install_windows.bat).
#
#  1. finds Python >= 3.10 (offers to install Python 3.12 with winget if there is none)
#  2. creates a private environment in %LOCALAPPDATA%\autogel\venv and installs autogel + GUI
#  3. creates "autogel" shortcuts on the desktop and in the Start menu
#  4. writes an uninstaller (Uninstall autogel.bat) into the install folder
#
# Options (for IT / testing):  -InstallDir <folder>   -NoShortcut   -Yes (no questions)
#                              -WithAI (also the optional AI second opinion, onnxruntime ~17 MB)
param(
    [string]$InstallDir = (Join-Path $env:LOCALAPPDATA "autogel"),
    [switch]$NoShortcut,
    [switch]$Yes,
    [switch]$WithAI
)
$ErrorActionPreference = "Stop"
$Source = Split-Path -Parent $PSScriptRoot          # the autogel folder (contains pyproject.toml)

function Say($msg, $color = "Gray") { Write-Host $msg -ForegroundColor $color }

function Find-Python {
    $candidates = @()
    if (Get-Command py -ErrorAction SilentlyContinue) {
        foreach ($v in "3.13", "3.12", "3.11", "3.10") { $candidates += ,@("py", "-$v") }
    }
    if (Get-Command python -ErrorAction SilentlyContinue) { $candidates += ,@("python") }
    foreach ($c in $candidates) {
        try {
            $exe = $c[0]; $args_ = @($c | Select-Object -Skip 1)
            $out = & $exe @args_ -c "import sys; print(sys.executable); print(sys.version_info >= (3, 10))" 2>$null
            if ($LASTEXITCODE -eq 0 -and $out.Count -ge 2 -and $out[1] -eq "True") { return $out[0] }
        } catch { }
    }
    return $null
}

Say "=== autogel installer ===" Cyan
Say "Source:  $Source"
Say "Target:  $InstallDir"
if (-not (Test-Path (Join-Path $Source "pyproject.toml"))) {
    Say "pyproject.toml not found next to this installer. Unzip the whole autogel folder first." Red
    exit 1
}

if ($InstallDir.Length -gt 70) {
    # PyQt6 contains ~130-character file paths; Windows' classic limit is 260 characters
    Say "The install folder path is very long ($($InstallDir.Length) characters). Windows may refuse" Red
    Say "the deep file names inside the Qt library. Choose a shorter folder with -InstallDir." Red
    exit 1
}

# ---------------------------------------------------------------- 1. Python
$py = Find-Python
if (-not $py) {
    Say "No Python 3.10 or newer found." Yellow
    if (Get-Command winget -ErrorAction SilentlyContinue) {
        $ok = $Yes -or ((Read-Host "Install Python 3.12 now with winget (official Windows package manager)? [Y/n]") -notmatch "^[nN]")
        if ($ok) {
            winget install --id Python.Python.3.12 -e --scope user --accept-package-agreements --accept-source-agreements
            $env:Path = [Environment]::GetEnvironmentVariable("Path", "User") + ";" + $env:Path
            $py = Find-Python
        }
    }
    if (-not $py) {
        Say "Please install Python 3.12 from https://www.python.org/downloads/ (tick 'Add python.exe to PATH')," Red
        Say "then run this installer again." Red
        exit 1
    }
}
Say "Python:  $py" Green

# ---------------------------------------------------------------- 2. environment + autogel
$venv = Join-Path $InstallDir "venv"
New-Item -ItemType Directory -Force -Path $InstallDir | Out-Null
if (-not (Test-Path (Join-Path $venv "Scripts\python.exe"))) {
    Say "Creating environment ..." Cyan
    & $py -m venv $venv
    if ($LASTEXITCODE -ne 0) { Say "Could not create the environment." Red; exit 1 }
}
$vpy = Join-Path $venv "Scripts\python.exe"
Say "Installing autogel and its viewer (napari) - this downloads ~300 MB and takes a few minutes ..." Cyan
& $vpy -m pip install --upgrade pip --quiet --disable-pip-version-check
$extras = if ($WithAI) { "gui,ai" } else { "gui" }
& $vpy -m pip install --quiet --disable-pip-version-check "$Source[$extras]"
if ($LASTEXITCODE -ne 0) { Say "Installation failed (see messages above). Check the internet connection / proxy." Red; exit 1 }
& $vpy -c "import autogel, napari; print('autogel', autogel.__version__, '| napari', napari.__version__)"
if ($LASTEXITCODE -ne 0) { Say "autogel was installed but cannot be imported." Red; exit 1 }

# ---------------------------------------------------------------- 3. shortcuts
$pyw = Join-Path $venv "Scripts\pythonw.exe"
$cmd = Join-Path $InstallDir "autogel.bat"
Set-Content -Encoding ASCII -Path $cmd -Value "@echo off`r`nstart `"`" `"$pyw`" -m autogel.cli gui %*"
if (-not $NoShortcut) {
    $shell = New-Object -ComObject WScript.Shell
    $targets = @([Environment]::GetFolderPath("Desktop"),
                 (Join-Path ([Environment]::GetFolderPath("Programs")) ""))
    foreach ($dir in $targets) {
        $lnk = $shell.CreateShortcut((Join-Path $dir "autogel.lnk"))
        $lnk.TargetPath = $pyw
        $lnk.Arguments = "-m autogel.cli gui"
        $lnk.WorkingDirectory = [Environment]::GetFolderPath("MyDocuments")
        $lnk.Description = "autogel - gel lane and band analysis"
        $lnk.Save()
    }
    Say "Shortcuts created on the desktop and in the Start menu." Green
}

# ---------------------------------------------------------------- 4. uninstaller
$un = Join-Path $InstallDir "Uninstall autogel.bat"
$unText = @"
@echo off
echo Removing autogel from $InstallDir ...
del "%USERPROFILE%\Desktop\autogel.lnk" 2>nul
del "%APPDATA%\Microsoft\Windows\Start Menu\Programs\autogel.lnk" 2>nul
cd /d "%TEMP%"
rmdir /s /q "$InstallDir"
echo Done. Your own ladders in %USERPROFILE%\.autogel were kept.
pause
"@
Set-Content -Encoding ASCII -Path $un -Value $unText

Say ""
Say "autogel is installed." Green
Say "Start it from the desktop icon 'autogel' (or run $cmd)."
Say "To update: download the new version and run install_windows.bat again."
Say "To remove: $un"
