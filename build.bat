@echo off
setlocal

where g++ >nul 2>nul
if %ERRORLEVEL% EQU 0 (
    echo [BUILD] Compilando con g++...
    g++ -std=c++17 -O2 -Iinclude src/parser.cpp src/tablas_hash.cpp src/max_heap.cpp src/grafo_adn.cpp src/main.cpp -o genomic_engine.exe
    if %ERRORLEVEL% EQU 0 (
        echo [OK] Compilado exitosamente: genomic_engine.exe
        exit /b 0
    )
)

set "VCVARS=C:\Program Files (x86)\Microsoft Visual Studio\2019\BuildTools\VC\Auxiliary\Build\vcvars64.bat"
if exist "%VCVARS%" (
    echo [BUILD] Compilando con MSVC cl.exe...
    call "%VCVARS%" >nul
    cl /EHsc /O2 /Iinclude src\parser.cpp src\tablas_hash.cpp src\max_heap.cpp src\grafo_adn.cpp src\main.cpp /Fe:genomic_engine.exe >nul
    del *.obj >nul 2>nul
    if exist genomic_engine.exe (
        echo [OK] Compilado exitosamente: genomic_engine.exe
        exit /b 0
    )
)

echo [ERROR] No se encontro compilador de C++ (g++ o MSVC).
exit /b 1
