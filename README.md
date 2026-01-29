# PPALI 2 (Peak-based PCA Analysis for Ligand Interactions 2)

**PPALI 2** is a web-based tool for analyzing NMR data using Principal Component Analysis (PCA). It specifically targets ligand-protein interaction studies, providing robust analysis for both 1D and 2D NMR titration series.

## Features

- **2D HSQC Analysis**: Automates PCA on 2D Peak lists (chemical shift perturbations).
- **1D Projection Analysis**: Projects 2D data into ^1^H and ^15^N dimensions for simplified visualization.
- **PCA Visualization**: Interactive 2D and 3D plots of scores and loadings.
- **Outlier Detection**: Statistical identification of significant chemical shift changes (Mahalanobis distance & P-value).
- **Binding Affinity (Kd) Fitting**:
  - **Traditional**: Standard independent fitting per residue.
  - **Universal (Relax)**: Global fitting shared across residues for robust estimation.

## Prerequisites

- **Python 3.10** or higher
- **Modern Web Browser** (Chrome, Firefox, Edge)

## Installation

### 1. Clone the Repository
```bash
git clone https://github.com/YourUsername/PPALI_2.git
cd PPALI_2
```

### 2. Setup Python Environment
It is recommended to use a virtual environment.

**Linux / macOS:**
```bash
python3 -m venv venv
source venv/bin/activate
pip install -r requirements.txt
```

**Windows:**
```powershell
python -m venv venv
venv\Scripts\activate
pip install -r requirements.txt
```

## Running the Application

### Using Convenience Scripts
- **Windows**: Double-click `run.bat`
- **Linux / macOS**: Run `./run.sh`

### Manual Start
```bash
# Ensure venv is active
uvicorn backend.main:app --host 127.0.0.1 --port 7777 --reload
```
Once running, open your browser to: [http://localhost:7777](http://localhost:7777)

## Project Structure

- `backend/`: FastAPI server and analysis logic (`analyzer.py`).
- `frontend/`: HTML/JS user interface.

## Contact
**Min June Yang**  
Email: minjune1590@kbsi.re.kr  
LinkedIn: [Min June Yang](https://www.linkedin.com/in/bionmr/)
