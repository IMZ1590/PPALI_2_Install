#!/bin/bash
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd "$DIR"

if [ ! -d "venv" ]; then
    echo "[INFO] Creating Virtual Environment..."
    python3 -m venv venv
fi

echo "[INFO] Installing Dependencies..."
./venv/bin/pip install -r requirements.txt

echo "[INFO] Starting Server..."
# Open browser (Linux/Mac compatible)
(sleep 2 && (xdg-open "http://localhost:7777" || open "http://localhost:7777")) &

./venv/bin/python3 -m uvicorn backend.main:app --host 127.0.0.1 --port 7777 --reload
