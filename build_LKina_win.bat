@echo off
:: ============================================================
::  LKina — Windows build script (MSYS2 / MinGW-w64)
::
::  LKina is a derivative of AutoDock Vina
::  Original work: Copyright (c) 2006-2010, The Scripps Research Institute
::  Licensed under the Apache License, Version 2.0
::  https://github.com/ccsb-scripps/AutoDock-Vina
::
::  Modifications: metal-aware inline AD4 map generation (LKina extensions)
::  See LICENSE file for full Apache-2.0 license text.
::
::  Requirements:
::    MSYS2  : https://www.msys2.org
::    Inside MSYS2 UCRT64 shell, run once:
::      pacman -S mingw-w64-ucrt-x86_64-gcc ^
::                mingw-w64-ucrt-x86_64-boost ^
::                make
::
::  Then run THIS SCRIPT from a standard CMD / PowerShell:
::    build_LKina_win.bat
::
::  Or run directly inside the MSYS2 shell:
::    bash build_LKina_linux.sh   (same script, auto-detects UCRT64)
::
::  Output: build\linux\release\LKina.exe
::
::  Alternatively — MinGW fully-static build via MSYS2 UCRT64 shell:
::    cd /path/to/LKina/build/linux/release
::    make -j4 LKina \
::         BASE=/ucrt64 \
::         BOOST_INCLUDE=/ucrt64/include \
::         BOOST_STATIC=y \
::         GPP=g++ \
::         C_PLATFORM="" \
::         C_OPTIONS="-O3 -DNDEBUG -std=c++14 -fPIC" \
::         BOOST_LIB_VERSION="" \
::         EXTRA_LDFLAGS="-static"
::  Result: LKina.exe depends only on KERNEL32.dll + Windows UCRT (Win10+ built-in)
:: ============================================================

:: Locate MSYS2
set MSYS2_ROOT=C:\msys64
if not exist "%MSYS2_ROOT%\usr\bin\bash.exe" (
    set MSYS2_ROOT=C:\msys2
)
if not exist "%MSYS2_ROOT%\usr\bin\bash.exe" (
    echo ERROR: MSYS2 not found at C:\msys64 or C:\msys2
    echo Install MSYS2 from https://www.msys2.org then run:
    echo   pacman -S mingw-w64-ucrt-x86_64-gcc mingw-w64-ucrt-x86_64-boost make
    exit /b 1
)

:: Get repo root (directory of this bat file)
set REPO_ROOT=%~dp0
:: Remove trailing backslash
if "%REPO_ROOT:~-1%"=="\" set REPO_ROOT=%REPO_ROOT:~0,-1%

:: Convert to MSYS2 path
set REPO_UNIX=%REPO_ROOT:\=/%
set REPO_UNIX=%REPO_UNIX:C:=/c%

echo Building LKina (MSYS2 UCRT64)...
set MSYSTEM=UCRT64
"%MSYS2_ROOT%\usr\bin\bash.exe" --login -c ^
    "cd '%REPO_UNIX%/build/linux/release' && touch dependencies 2>/dev/null; make -j%NUMBER_OF_PROCESSORS% LKina BASE=/ucrt64 BOOST_INCLUDE=/ucrt64/include BOOST_STATIC=y GPP=g++ C_PLATFORM='' C_OPTIONS='-O3 -DNDEBUG -std=c++14 -fPIC' BOOST_LIB_VERSION='' EXTRA_LDFLAGS='-static' 2>&1 && echo 'Build OK: LKina.exe' || echo 'Build FAILED'"

pause
