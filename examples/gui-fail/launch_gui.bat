@echo off
REM GridGeneration GUI Launcher for Windows
REM Double-click this file to launch the GUI

echo ========================================
echo   GridGeneration.jl GUI Launcher
echo ========================================
echo.

REM Check if Julia is in PATH
where julia >nul 2>nul
if %ERRORLEVEL% NEQ 0 (
    echo ERROR: Julia not found in PATH
    echo.
    echo Please install Julia or add it to your PATH:
    echo https://julialang.org/downloads/
    echo.
    pause
    exit /b 1
)

echo Julia found!
echo.

REM Navigate to GUI directory
cd /d "%~dp0"

echo Starting GUI launcher...
echo.

julia launcher.jl

pause
