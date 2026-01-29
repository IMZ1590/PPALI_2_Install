# PPALI 2 (Peak-based PCA Analysis for Ligand Interactions 2)

**PPALI 2** is a web-based tool for analyzing NMR data using Principal Component Analysis (PCA). It specifically targets ligand-protein interaction studies, providing robust analysis for both 1D and 2D NMR titration series.

## Features

<ul>
  <li>
    <strong>2D HSQC Analysis</strong>: Automates multi-parametric PCA on titration data extracted from <strong>2rr files</strong>, specifically analyzing chemical shift perturbations and peak intensity changes.
  </li>
  <li>
    <strong>1D Projection Analysis</strong>: Projects complex 2D spectra into <strong><sup>1</sup>H and <sup>15</sup>N dimensions</strong> to provide simplified visualization of residue-specific shifts.
  </li>
  <li>
    <strong>PCA Visualization</strong>: Provides interactive <strong>2D plots</strong> of scores and loadings for intuitive identification of binding patterns (3D plotting is excluded for clearer data interpretation).
  </li>
  <li>
    <strong>Outlier Detection</strong>: Statistical identification of significant binding-induced changes using <strong>Mahalanobis distance and P-value</strong> calculations.
  </li>
  <li>
    <strong>Binding Affinity (K<sub>d</sub>) Fitting</strong>:
    <ul>
      <li><strong>Traditional</strong>: Standard non-linear regression fitting using 2D or 1D spectral data for individual residues.</li>
      <li><strong>Relax</strong>: <strong>Robust global fitting</strong> shared across multiple residues to ensure stable estimation when individual peak data is ambiguous or noisy.</li>
    </ul>
  </li>
</ul>

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
