from fastapi import FastAPI, Form, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import FileResponse
from fastapi.staticfiles import StaticFiles
from typing import List
import os
import re
import numpy as np
import sys

# Ensure backend directory is in path for imports
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from analyzer import load_nmr_data, run_pca_and_fit, parse_bruker_param, run_3way_analysis, fit_with_regime_constraint

app = FastAPI(title="PPALI 2")

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_methods=["*"],
    allow_headers=["*"],
)

from fastapi.responses import JSONResponse
from fastapi.requests import Request

@app.exception_handler(Exception)
async def global_exception_handler(request: Request, exc: Exception):
    import traceback
    error_msg = f"{str(exc)}\n{traceback.format_exc()}"
    print(f"CRITICAL ERROR: {error_msg}")
    return JSONResponse(
        status_code=500,
        content={"detail": str(exc), "trace": traceback.format_exc()},
    )

# Robust cross-platform root detection
ALLOWED_ROOT = os.environ.get("ALLOWED_ROOT", os.path.abspath(os.sep))

def ensure_safe_path(path: str):
    abs_path = os.path.abspath(path)
    if not abs_path.startswith(os.path.abspath(ALLOWED_ROOT)):
        raise HTTPException(status_code=403, detail=f"Access denied: Paths must be within {ALLOWED_ROOT}")
    if not os.path.exists(abs_path):
        raise HTTPException(status_code=404, detail=f"Path not found: {path}")
    return abs_path

# Logo Endpoint
@app.get("/logo")
async def get_logo():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    logo_path = os.path.join(script_dir, "KBSI_Logo.png")
    if os.path.exists(logo_path):
        return FileResponse(logo_path)
    return FileResponse("KBSI_Logo.png")

import nmrglue as ng

@app.post("/scan")
async def scan_directory(path: str = Form(...), dim: str = Form("2D")):
    safe_path = ensure_safe_path(path)
    found_experiments = []
    target_file = '2rr' if dim == '2D' else '1r'
    
    for root, dirs, files in os.walk(safe_path):
        if target_file in files:
            pdata_path = root
            try:
                exp_path = os.path.dirname(os.path.dirname(pdata_path))
                acqus_path = os.path.join(exp_path, "acqus")
                procs_path = os.path.join(pdata_path, "procs")
                
                def safe_parse(path, key, default=0):
                    val = parse_bruker_param(path, key)
                    return val if val is not None else default
                
                nc_proc = safe_parse(procs_path, "NC_proc")
                sw = safe_parse(procs_path, "SW_p") 
                o1p = safe_parse(procs_path, "OFFSET")
                td = safe_parse(procs_path, "SI")
                ns = safe_parse(acqus_path, "NS")
                rg = safe_parse(acqus_path, "RG")
                
                try:
                    exp_no = os.path.basename(exp_path)
                    exp_name = os.path.basename(os.path.dirname(exp_path))
                    display_name = f"{exp_name}/{exp_no}"
                except:
                    display_name = pdata_path

                found_experiments.append({
                    "path": pdata_path, 
                    "name": display_name, 
                    "nc_proc": nc_proc,
                    "sw": sw, "o1p": o1p, "td": td, "ns": ns, "rg": rg, "dim": dim
                })
            except Exception as e:
                continue

    found_experiments.sort(key=lambda x: x['name'])
    return {"experiments": found_experiments}

import shutil
import tempfile
from fastapi import UploadFile, File

@app.post("/upload-scan")
async def upload_scan_directory(files: List[UploadFile] = File(...), dim: str = Form("2D")):
    temp_dir = tempfile.mkdtemp(prefix="ppali_upload_")
    
    for file in files:
        if not file.filename: continue
        clean_name = file.filename.replace("..", "")
        dest_path = os.path.join(temp_dir, clean_name)
        os.makedirs(os.path.dirname(dest_path), exist_ok=True)
        with open(dest_path, "wb") as buffer:
            shutil.copyfileobj(file.file, buffer)
    
    found_experiments = []
    target_file = '2rr' if dim == '2D' else '1r'
    
    for root, dirs, filenames in os.walk(temp_dir):
        if target_file in filenames:
            pdata_path = root
            try:
                exp_path = os.path.dirname(os.path.dirname(pdata_path))
                acqus_path = os.path.join(exp_path, "acqus")
                procs_path = os.path.join(pdata_path, "procs")
                
                def safe_parse(path, key, default=0):
                    val = parse_bruker_param(path, key)
                    return val if val is not None else default
                
                nc_proc = safe_parse(procs_path, "NC_proc")
                sw = safe_parse(procs_path, "SW_p") 
                o1p = safe_parse(procs_path, "OFFSET")
                td = safe_parse(procs_path, "SI")
                ns = safe_parse(acqus_path, "NS")
                rg = safe_parse(acqus_path, "RG")
                
                try:
                    rel_path = os.path.relpath(exp_path, temp_dir)
                    display_name = rel_path if rel_path != "." else "Uploaded Data"
                except:
                    display_name = "Uploaded Experiment"

                found_experiments.append({
                    "path": pdata_path, 
                    "name": display_name, 
                    "nc_proc": nc_proc,
                    "sw": sw, "o1p": o1p, "td": td, "ns": ns, "rg": rg, "dim": dim
                })
            except Exception as e:
                continue
                
    found_experiments.sort(key=lambda x: x['name'])
    return {"experiments": found_experiments, "temp_id": temp_dir}

@app.post("/analyze")
async def analyze_data(
    paths: List[str] = Form(...), concs: List[float] = Form(...),
    protein_conc: float = Form(50.0), regime: str = Form("Intermediate"),
    ppm_start: float = Form(None), ppm_end: float = Form(None),
    no_fitting: bool = Form(False)
):
    if len(paths) != len(concs):
        raise HTTPException(status_code=400, detail="Paths and concentrations mismatch")
    try:
        sorted_indices = sorted(range(len(concs)), key=lambda k: concs[k])
        sorted_paths = [paths[i] for i in sorted_indices]
        sorted_concs = [concs[i] for i in sorted_indices]
        for path in sorted_paths:
            ensure_safe_path(path)
        ppm_range = (ppm_start, ppm_end) if ppm_start is not None and ppm_end is not None else None
        analysis_res = run_3way_analysis(sorted_paths, sorted_concs, protein_conc, regime, ppm_range=ppm_range, no_fitting=no_fitting)
        return analysis_res
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

@app.post("/refit_regime")
async def refit_regime(
    pc_scores: List[float] = Form(...), concentrations: List[float] = Form(...),
    regime: str = Form(...), protein_conc: float = Form(50.0)
):
    try:
        return fit_with_regime_constraint(np.array(concentrations), np.array(pc_scores), protein_conc, regime=regime)
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

# Mount Frontend for Static Serving (Must be LAST to avoid shadowing API)
frontend_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "frontend")
if os.path.exists(frontend_dir):
    app.mount("/", StaticFiles(directory=frontend_dir, html=True), name="frontend")
else:
    print(f"Warning: Frontend directory not found at {frontend_dir}")

# Startup Event for Auto-Open
@app.on_event("startup")
async def startup_event():
    import webbrowser
    import threading
    import time
    
    port = int(os.environ.get("PORT", 7777))
    url = f"http://127.0.0.1:{port}"
    
    def open_browser():
        # Give server a moment to bind
        time.sleep(1.0)
        print(f"[INFO] Opening Browser: {url}")
        try:
            # Force NEW TAB (new=2)
            webbrowser.open(url, new=2)
        except Exception as e:
            print(f"[WARN] Failed to open browser: {e}")

    threading.Thread(target=open_browser, daemon=True).start()

if __name__ == "__main__":
    import uvicorn
    port = int(os.environ.get("PORT", 7777))
    # Listen on localhost (127.0.0.1) explicitly
    # Pass 'app' object to avoid re-import execution loop
    uvicorn.run(app, host="127.0.0.1", port=port)