@echo off
cd /d "%~dp0"
if not exist genomic_engine.exe (
    call build.bat
)
python gui\app.py
