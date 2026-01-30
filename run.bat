@echo off
cd /d "%~dp0"
title PALI 2 Server

IF NOT EXIST "venv" (
    echo [INFO] Creating Virtual Environment...
    python -m venv venv
)

echo [INFO] Installing Dependencies...
venv\Scripts\python -m pip install -r requirements.txt

echo [INFO] Starting Server...
start "" "http://localhost:7777"
venv\Scripts\python -m uvicorn backend.main:app --host 127.0.0.1 --port 7777 --reload
pause
