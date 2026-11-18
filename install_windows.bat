@echo off
rem autogel installer for Windows - double-click this file.
rem Installs autogel into %LOCALAPPDATA%\autogel and creates desktop / Start-menu shortcuts.
rem No administrator rights, git or command-line knowledge needed.
powershell.exe -NoProfile -ExecutionPolicy Bypass -File "%~dp0tools\install_windows.ps1" %*
echo.
pause
