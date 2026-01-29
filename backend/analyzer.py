import numpy as np
import nmrglue as ng
from sklearn.decomposition import PCA
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d
import os
import re
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import io
import base64
def parse_bruker_param(filepath, param_name):
    if not os.path.exists(filepath):
        return None
    try:
        with open(filepath, 'r', encoding='utf-8', errors='ignore') as f:
            content = f.read()
            pattern = rf'##\${param_name}=\s*([-\d.eE+]+)'
            match = re.search(pattern, content)
            if match:
                val = match.group(1)
                try:
                    return float(val)
                except ValueError:
                    return val
    except Exception as e:
        print(f"Error parsing {filepath} for {param_name}: {e}")
    return None
def load_nmr_data(pdata_path):
    dic, data = ng.bruker.read_pdata(pdata_path)
    procs_path = os.path.join(pdata_path, "procs")
    exp_path = os.path.dirname(os.path.dirname(pdata_path))
    acqus_path = os.path.join(exp_path, "acqus")
    nc_proc = parse_bruker_param(procs_path, "NC_proc")
    ns = parse_bruker_param(acqus_path, "NS")
    rg = parse_bruker_param(acqus_path, "RG")
    td_f2 = parse_bruker_param(acqus_path, "TD") # Note: this might need more careful handling for 2D
    if nc_proc is None: nc_proc = 0
    if ns is None: ns = 1
    if rg is None: rg = 1
    if '3rr' in pdata_path:
        if data.ndim == 3:
            planes = []
            delays = extract_delays(exp_path)
            if not delays or len(delays) != data.shape[0]:
                 pass
            for i in range(data.shape[0]):
                plane = data[i, :, :].astype(np.float64) * (2.0 ** nc_proc) / (np.sqrt(ns) * rg)
                planes.append(plane)
            return {
                "type": "3D",
                "data": planes,
                "delays": delays,
                "params": {"nc_proc": nc_proc, "ns": ns, "rg": rg},
                "shape": (data.shape[1], data.shape[2])
            }
    normalized_data = data.astype(np.float64) * (2.0 ** nc_proc) / (np.sqrt(ns) * rg)
    return {
        "type": "2D",
        "data": normalized_data,
        "params": {
            "nc_proc": nc_proc,
            "ns": ns,
            "rg": rg,
            "td": td_f2
        },
        "shape": data.shape
    }
def extract_delays(exp_path):
    for fname in ['vdlist', 'vclist']:
        fpath = os.path.join(exp_path, fname)
        if os.path.exists(fpath):
            try:
                with open(fpath, 'r') as f:
                    lines = [l.strip() for l in f.readlines()]
                    delays = []
                    for line in lines:
                        if not line or line.startswith('#'): continue
                        parts = line.split()
                        if parts:
                            delays.append(float(parts[0]))
                    return delays
            except:
                pass
    return []
def get_bound_fraction(L, Kd, P_tot):
    return L / (Kd + L)
def universal_exchange_model(L, A, B, C, Kd, P_tot):
    f = get_bound_fraction(L, Kd, P_tot)
    return A + B * f + C * f * (1 - f)
def classify_exchange_regime(B, C, threshold=3.0):
    ratio = abs(B) / (abs(C) + 1e-10)
    if ratio > threshold:
        return "Slow", ratio
    elif ratio < 1.0/threshold:
        return "Fast", ratio
    else:
        return "Intermediate", ratio
def traditional_binding_model(L, Kd, amp, offset, P_tot):
    f = get_bound_fraction(L, Kd, P_tot)
    return offset + amp * f
def fit_universal_model(concentrations, scores, protein_conc):
    y = np.array(scores)
    x = np.array(concentrations)
    kd_guess = max(1.0, np.median(x) if len(x) > 0 else 10.0)
    amp_guess = np.max(y) - np.min(y)
    if y[-1] < y[0]: amp_guess *= -1
    off_guess = np.mean(y)
    p0_trad = [kd_guess, amp_guess, off_guess]
    bounds_trad = ([0.01, -np.inf, -np.inf], [np.inf, np.inf, np.inf])
    def trad_wrapper(conc, Kd, amp, offset):
        return traditional_binding_model(conc, Kd, amp, offset, protein_conc)
    try:
        popt_trad, pcov_trad = curve_fit(trad_wrapper, x, y, p0=p0_trad, bounds=bounds_trad, maxfev=10000)
        Kd_trad, amp_trad, off_trad = popt_trad
        Kd_err_trad = np.sqrt(np.diag(pcov_trad))[0]
        max_c = max(x) if len(x) > 0 else 100
        fit_x_trad = np.linspace(0, max_c * 1.5, 100)
        fit_y_trad = traditional_binding_model(fit_x_trad, Kd_trad, amp_trad, off_trad, protein_conc)
        trad_result = {
            "success": True,
            "kd": Kd_trad,
            "kd_err": Kd_err_trad,
            "fit_x": fit_x_trad.tolist(),
            "fit_y": fit_y_trad.tolist()
        }
    except:
        trad_result = {"success": False, "kd": 0, "kd_err": 0}
    univ_result = fit_with_regime_constraint(x, y, protein_conc, regime="auto")
    return {
        "traditional": trad_result,
        "universal": univ_result
    }
def fit_with_regime_constraint(x, y, protein_conc, regime="Fast"):
    x = np.array(x)
    y = np.array(y)
    A_guess = np.mean(y)
    amp_total = np.max(y) - np.min(y)
    kd_guess = max(1.0, np.median(x) if len(x) > 0 else 10.0)
    C_guess = abs(amp_total) * 2.0  
    p0 = [A_guess, C_guess, kd_guess]
    bounds = ([-np.inf, -np.inf, 0.01], [np.inf, np.inf, np.inf])
    def model_fast(conc, A, C, Kd):
        return universal_exchange_model(conc, A, 0.0, C, Kd, protein_conc)
    try:
        popt, pcov = curve_fit(model_fast, x, y, p0=p0, bounds=bounds, maxfev=10000)
        A_fit, C_fit, Kd_fit = popt
        B_fit = 0.0
        errors = np.sqrt(np.diag(pcov))
        Kd_err = errors[2]
        regime_class = "Relax" 
        max_c = max(x) if len(x) > 0 else 100
        fit_x = np.linspace(0, max_c * 1.5, 100)
        fit_y = universal_exchange_model(fit_x, A_fit, B_fit, C_fit, Kd_fit, protein_conc)
        return {
            "success": True,
            "kd": Kd_fit,
            "kd_err": Kd_err,
            "A": A_fit,
            "B": B_fit,
            "C": C_fit,
            "regime": regime_class,
            "fit_x": fit_x.tolist(),
            "fit_y": fit_y.tolist()
        }
    except Exception as e:
        return {"success": False, "kd": 0, "kd_err": 0, "error": f"Relax fit failed: {str(e)}"}
def get_ppm_scale(dic, data, axis=0):
    try:
        udic = ng.bruker.guess_udic(dic, data)
        uc = ng.fileio.fileiobase.uc_from_udic(udic, dim=axis)
        return uc.ppm_scale()
    except:
        return np.linspace(0, 1, data.shape[axis])
def load_projection_data(pdata_path):
    try:
        p_2rr = os.path.join(pdata_path, '2rr')
        p_1r = os.path.join(pdata_path, '1r')
        p_1r_stash = os.path.join(pdata_path, '1r.stash')
        stashed = False
        if os.path.exists(p_2rr) and os.path.exists(p_1r):
            try:
                os.rename(p_1r, p_1r_stash)
                stashed = True
            except:
                pass # Permission error?
        try:
            dic, data = ng.bruker.read_pdata(pdata_path)
        finally:
            if stashed and os.path.exists(p_1r_stash):
                os.rename(p_1r_stash, p_1r)
        if data.ndim == 1 and data.size == 10000:
            data = data.reshape((100, 100))
        if data.ndim != 1 and data.ndim != 2: 
            return {"error": f"Invalid dimensions {data.shape} (Expected 1D or 2D)"}
        procs_path = os.path.join(pdata_path, "procs")
        exp_path = os.path.dirname(os.path.dirname(pdata_path))
        acqus_path = os.path.join(exp_path, "acqus")
        nc_proc = parse_bruker_param(procs_path, "NC_proc")
        ns = parse_bruker_param(acqus_path, "NS")
        rg = parse_bruker_param(acqus_path, "RG")
        if nc_proc is None: nc_proc = 0
        if ns is None: ns = 1
        if rg is None: rg = 1
        if rg is None: rg = 1
        sf = parse_bruker_param(procs_path, "SF")
        if sf is None: sf = parse_bruker_param(acqus_path, "SFO1")
        if sf is None: sf = parse_bruker_param(acqus_path, "BF1")
        if sf is None: sf = 600.0 # Last resort fallback
        start_ppm = None
        end_ppm = None
        offset = parse_bruker_param(procs_path, "OFFSET")
        sw_p = parse_bruker_param(procs_path, "SW_p")
        if offset is not None and sw_p is not None:
             if sw_p > 500: sw_p = sw_p / sf 
             start_ppm = offset
             end_ppm = offset - sw_p
        if start_ppm is None:
             sw_acq = parse_bruker_param(acqus_path, "SW")
             o1p = parse_bruker_param(acqus_path, "O1P")
             if o1p is not None and sw_acq is not None:
                 sw_ppm = sw_acq
                 if sw_acq > 100: sw_ppm = sw_acq / sf
                 start_ppm = o1p + (sw_ppm / 2.0)
                 end_ppm = o1p - (sw_ppm / 2.0)
        if start_ppm is None:
             start_ppm = 10.0
             end_ppm = -2.0
        factor = (2.0 ** nc_proc) / (np.sqrt(ns) * rg)
        if data.ndim == 1:
            norm_data = data.astype(np.float64) * factor
            size = norm_data.shape[0]
            ppm_axis = np.linspace(start_ppm, end_ppm, size)
            return {
                "full_data": norm_data, 
                "expanded_n": None,     
                "expanded_h": norm_data, 
                "ppm_axis": ppm_axis,
                "params": {"nc_proc": nc_proc, "ns": ns, "rg": rg},
                "is_1d": True
            }
        norm_data = data.astype(np.float64) * factor
        size_f1, size_f2 = norm_data.shape
        ppm_axis = np.linspace(start_ppm, end_ppm, size_f2)
        ppm_axis_15n = None
        proc2s_path = os.path.join(pdata_path, "proc2s")
        acqu2s_path = os.path.join(exp_path, "acqu2s")
        sf_2 = parse_bruker_param(proc2s_path, "SF") # try proc first
        if sf_2 is None: sf_2 = parse_bruker_param(acqu2s_path, "SFO2") 
        if sf_2 is None: sf_2 = parse_bruker_param(acqu2s_path, "BF2")
        if sf_2 is None: sf_2 = 60.0 
        if os.path.exists(proc2s_path):
            offset_2 = parse_bruker_param(proc2s_path, "OFFSET")
            sw_p_2 = parse_bruker_param(proc2s_path, "SW_p")
            if offset_2 is not None and sw_p_2 is not None:
                if sw_p_2 > 500: sw_p_2 = sw_p_2 / sf_2 # Hz check
                ppm_axis_15n = np.linspace(offset_2, offset_2 - sw_p_2, size_f1)
        if ppm_axis_15n is None:
            sw_2_acq = parse_bruker_param(acqu2s_path, "SW")
            o1p_2 = parse_bruker_param(acqu2s_path, "O1P")
            if o1p_2 is not None and sw_2_acq is not None:
                 sw_2_ppm = sw_2_acq
                 if sw_2_acq > 500: sw_2_ppm = sw_2_acq / sf_2
                 start_2 = o1p_2 + (sw_2_ppm / 2.0)
                 end_2 = o1p_2 - (sw_2_ppm / 2.0)
                 ppm_axis_15n = np.linspace(start_2, end_2, size_f1)
        if ppm_axis_15n is None:
            acqu2s_path = os.path.join(exp_path, "acqu2s")
            o1p_2 = parse_bruker_param(acqu2s_path, "O1P")
            sw_2 = parse_bruker_param(acqu2s_path, "SW")
            sf_2 = parse_bruker_param(acqu2s_path, "SFO2") 
            if sf_2 is None: sf_2 = parse_bruker_param(acqu2s_path, "BF2") # Spectrometer Freq
            if sf_2 is None: sf_2 = 60.0 # Approximate 15N freq on 600MHz
            if sw_2 is not None:
                 if sw_2 > 500: sw_2_ppm = sw_2 / sf_2
                 else: sw_2_ppm = sw_2
            else: sw_2_ppm = 40.0
            if o1p_2 is not None:
                 start_2 = o1p_2 + (sw_2_ppm / 2.0)
                 end_2 = o1p_2 - (sw_2_ppm / 2.0)
                 ppm_axis_15n = np.linspace(start_2, end_2, size_f1)
        proj_n_mean = np.mean(norm_data, axis=1) 
        proj_h_mean = np.mean(norm_data, axis=0) 
        return {
            "full_data": norm_data,
            "expanded_n": proj_n_mean, 
            "expanded_h": proj_h_mean,
            "ppm_axis": ppm_axis,
            "ppm_axis_15n": ppm_axis_15n,
            "params": {"nc_proc": nc_proc, "ns": ns, "rg": rg},
            "is_1d": False
        }
    except Exception as e:
        import traceback
        return {"error": f"Error loading {os.path.basename(pdata_path)}: {str(e)}"}
def align_spectra(data_list, ppm_axes_list):
    if not ppm_axes_list or len(ppm_axes_list) != len(data_list):
        return data_list, None # Cannot align
    global_min_ppm = -np.inf
    global_max_ppm = np.inf
    for ppm in ppm_axes_list:
        p_min, p_max = np.min(ppm), np.max(ppm)
        if p_min > global_min_ppm: global_min_ppm = p_min
        if p_max < global_max_ppm: global_max_ppm = p_max
    if global_max_ppm <= global_min_ppm:
        return data_list, None # No overlap
    max_points = 0
    for ppm in ppm_axes_list:
        mask = (ppm >= global_min_ppm) & (ppm <= global_max_ppm)
        count = np.sum(mask)
        if count > max_points: max_points = count
    if max_points < 2: max_points = 100 # Fallback
    common_axis = np.linspace(global_max_ppm, global_min_ppm, max_points)
    aligned_data = []
    for i, data in enumerate(data_list):
        ppm = ppm_axes_list[i]
        if data.ndim == 1:
            try:
                f = interp1d(ppm, data, kind='linear', bounds_error=False, fill_value=0.0)
                aligned_spec = f(common_axis)
                aligned_data.append(aligned_spec)
            except Exception as e:
                print(f"Alignment failed for index {i}: {e}")
                aligned_data.append(data) # Fallback
        elif data.ndim == 2:
             new_shape = (data.shape[0], len(common_axis))
             aligned_mat = np.zeros(new_shape)
             for r in range(data.shape[0]):
                 row = data[r, :]
                 f = interp1d(ppm, row, kind='linear', bounds_error=False, fill_value=0.0)
                 aligned_mat[r, :] = f(common_axis)
             aligned_data.append(aligned_mat)
    return aligned_data, common_axis
def run_pca_and_fit(data_list, concentrations, protein_conc=50.0, regime="Intermediate", ppm_axes=None, ppm_range=None, exclude_ranges=None, no_fitting=False, ppm_15n=None, shape_2d=None):
    if not data_list or len(data_list) < 2:
        return {"success": False, "error": "Insufficient data"}
    common_axis = None
    if ppm_axes and len(ppm_axes) == len(data_list):
        try:
             data_list, common_axis = align_spectra(data_list, ppm_axes)
        except Exception as e:
             return {"success": False, "error": f"Alignment failed: {str(e)}"}
    elif ppm_axes and len(ppm_axes) > 0:
        common_axis = ppm_axes[0]
    if ppm_range and common_axis is not None:
        try:
            p_start, p_end = ppm_range
            p_min = min(p_start, p_end)
            p_max = max(p_start, p_end)
            mask = (common_axis >= p_min) & (common_axis <= p_max)
            if not np.any(mask):
                return {"success": False, "error": f"PPM Range {ppm_range} out of data bounds."}
            sliced_list = []
            for d in data_list:
                if d.ndim == 1:
                    sliced_list.append(d[mask])
                elif d.ndim == 2:
                    sliced_list.append(d[:, mask])
            data_list = sliced_list
            common_axis = common_axis[mask]
            if shape_2d is not None and len(sliced_list) > 0:
                 new_cols = sliced_list[0].shape[1] if sliced_list[0].ndim == 2 else sliced_list[0].shape[0]
                 shape_2d = (shape_2d[0], new_cols)
            elif sliced_list and len(sliced_list) > 0:
                 new_shape = sliced_list[0].shape
                 if len(new_shape) == 2:
                      shape_2d = new_shape
        except Exception as e:
             return {"success": False, "error": f"Slicing failed: {str(e)}"}
    if exclude_ranges and common_axis is not None:
        try:
            mask = np.ones(len(common_axis), dtype=bool)
            for (p_start, p_end) in exclude_ranges:
                p_min = min(p_start, p_end)
                p_max = max(p_start, p_end)
                bad_indices = (common_axis >= p_min) & (common_axis <= p_max)
                mask[bad_indices] = False
            sliced_list = []
            for d in data_list:
                if d.ndim == 1:
                    sliced_list.append(d[mask])
                elif d.ndim == 2:
                    sliced_list.append(d[:, mask])
            data_list = sliced_list
            common_axis = common_axis[mask]
            if sliced_list and len(sliced_list) > 0:
                 new_shape = sliced_list[0].shape
                 if len(new_shape) == 2:
                      shape_2d = new_shape
        except Exception as e:
             return {"success": False, "error": f"Exclusion failed: {str(e)}"}
    try:
        concs_arr = np.array(concentrations, dtype=float)
        first_shape = data_list[0].shape
        if shape_2d is None and len(first_shape) == 2:
            shape_2d = first_shape
        for i, d in enumerate(data_list):
            if d.shape != first_shape:
                return {"success": False, "error": f"Data dimensions mismatch. File {i+1} has shape {d.shape}, expected {first_shape}. All spectra must have the same dimensions."}
        sample = data_list[0]
        if sample.ndim > 1:
            X = np.array([d.flatten() for d in data_list])
        else:
            X = np.vstack(data_list)
        n_comp = 2
        pca = PCA(n_components=n_comp)
        scores = pca.fit_transform(X) # (N_samples, n_comp)
        var = pca.explained_variance_ratio_.tolist()
        pc1 = scores[:, 0]
        pc2 = scores[:, 1] if scores.shape[1] > 1 else np.zeros_like(pc1)
        mean_spec = np.mean(X, axis=0)
        final_mean = mean_spec
        final_ppm = common_axis
        final_ppm_15n = ppm_15n
        final_shape = shape_2d
        final_all = X.tolist() # Default for 1D
        if shape_2d:
            final_all = [] 
            rows, cols = shape_2d
            MAX_DIM = 256 # Constraint for performance of overlays
            def downsample_2d(flat_arr, r, c, r_step, c_step, nr, nc):
                try:
                    mat = flat_arr.reshape(r, c)
                    mat = mat[:nr*r_step, :nc*c_step]
                    return mat.reshape(nr, r_step, nc, c_step).max(axis=(1, 3)).flatten()
                except:
                    return flat_arr # Fallback
            if rows > MAX_DIM or cols > MAX_DIM:
                try:
                    r_step = max(1, int(np.ceil(rows / MAX_DIM)))
                    c_step = max(1, int(np.ceil(cols / MAX_DIM)))
                    new_rows = rows // r_step
                    new_cols = cols // c_step
                    final_mean = downsample_2d(mean_spec, rows, cols, r_step, c_step, new_rows, new_cols)
                    final_shape = (new_rows, new_cols)
                    final_all = []
                    for spec in X:
                        ds_spec = downsample_2d(spec, rows, cols, r_step, c_step, new_rows, new_cols)
                        final_all.append(ds_spec.tolist())
                    if final_ppm is not None:
                        final_ppm = final_ppm[:new_cols*c_step:c_step]
                    if final_ppm_15n is not None:
                        final_ppm_15n = final_ppm_15n[:new_rows*r_step:r_step]
                except Exception as e:
                    print(f"Downsampling error: {e}")
                    pass # Fallback
            else:
                 final_all = X.tolist()
        hsqc_img_b64 = None
        spectra_payload = {
            "ppm": final_ppm.tolist() if final_ppm is not None else [],
            "ppm_15n": final_ppm_15n.tolist() if final_ppm_15n is not None else [],
            "shape_2d": final_shape if final_shape else [],
            "mean_spectrum": final_mean.tolist(),
            "all_spectra": final_all,
            "hsqc_image": hsqc_img_b64
        }
        while len(var) < 3: var.append(0.0)
        if no_fitting:
             dummy_fit = {"success": False, "error": "Fitting skipped by user"}
             return {
                "success": True,
                "scores": pc1.tolist(),
                "pc2_scores": pc2.tolist(),
                "pc3_scores": pc3.tolist(),
                "variance": var,
                "pc1_fit": {"traditional": dummy_fit, "universal": dummy_fit},
                "pc2_fit": {"traditional": dummy_fit, "universal": dummy_fit},
                "pc3_fit": {"traditional": dummy_fit, "universal": dummy_fit},
                "pc3_fit": {"traditional": dummy_fit, "universal": dummy_fit},
                "spectra_check": spectra_payload
            }
        def fit_trad(scores_arr):
            y = np.array(scores_arr)
            x = concs_arr
            kd_guess = max(1.0, np.median(x) if len(x) > 0 else 10.0)
            amp_guess = np.max(y) - np.min(y)
            if len(y) > 1 and y[-1] < y[0]: amp_guess *= -1
            off_guess = np.mean(y)
            p0 = [kd_guess, amp_guess, off_guess]
            bounds = ([0.01, -np.inf, -np.inf], [np.inf, np.inf, np.inf])
            def wrapper(conc, Kd, amp, offset):
                return traditional_binding_model(conc, Kd, amp, offset, protein_conc)
            try:
                popt, pcov = curve_fit(wrapper, x, y, p0=p0, bounds=bounds, maxfev=10000)
                Kd_fit, amp_fit, off_fit = popt
                Kd_err = np.sqrt(np.diag(pcov))[0]
                max_c = max(x) if len(x) > 0 else 100
                fit_x = np.linspace(0, max_c * 1.5, 100)
                fit_y = traditional_binding_model(fit_x, Kd_fit, amp_fit, off_fit, protein_conc)
                return {"success": True, "kd": Kd_fit, "kd_err": Kd_err, "fit_x": fit_x.tolist(), "fit_y": fit_y.tolist()}
            except:
                return {"success": False, "kd": 0, "kd_err": 0}
        return {
            "success": True,
            "scores": pc1.tolist(),
            "pc2_scores": pc2.tolist(),
            "variance": var,
            "pc1_fit": {"traditional": fit_trad(pc1), "universal": fit_with_regime_constraint(concs_arr, pc1, protein_conc, regime)},
            "pc2_fit": {"traditional": fit_trad(pc2), "universal": fit_with_regime_constraint(concs_arr, pc2, protein_conc, regime)},
            "spectra_check": spectra_payload
        }
    except Exception as e:
        return {"success": False, "error": str(e)}
def run_3way_analysis(path_list, concentrations, protein_conc, regime="Intermediate", ppm_range=None, no_fitting=False):
    list_2d = []
    list_h = []
    list_n = []
    ppm_axes_list = []
    has_ppm = True
    for path in path_list:
        res = load_projection_data(path)
        if res and 'ppm_axis' in res:
            ppm_axes_list.append(res['ppm_axis'])
        else:
            has_ppm = False
    list_2d = []
    list_h = []
    list_n = []
    ppm_axes = [] # Aligned order
    errors = []
    valid_indices = []
    final_list_2d = []
    final_list_h = []
    final_list_n = []
    final_ppm_axes = []     # For 1H
    final_ppm_axes_15n = [] # For 15N
    ppm_15n_ref = None
    shape_2d_ref = None
    for i, path in enumerate(path_list):
        res = load_projection_data(path)
        if res and 'full_data' in res:
             valid_indices.append(i)
             if ppm_range and not res.get('is_1d'):
                 data_2d = res['full_data']
                 ppm_1h = res['ppm_axis']
                 p_min, p_max = min(ppm_range), max(ppm_range)
                 mask = (ppm_1h >= p_min) & (ppm_1h <= p_max)
                 if np.any(mask):
                     res['expanded_n'] = np.mean(data_2d[:, mask], axis=1)
             if res.get('is_1d'):
                 final_list_h.append(res['expanded_h'])
             else:
                 final_list_2d.append(res['full_data'])
                 final_list_h.append(res['expanded_h'])
                 final_list_n.append(res['expanded_n'])
                 if 'ppm_axis_15n' in res:
                     final_ppm_axes_15n.append(res['ppm_axis_15n'])
                 if shape_2d_ref is None:
                     if 'ppm_axis_15n' in res:
                         ppm_15n_ref = res['ppm_axis_15n']
                     if 'full_data' in res:
                         shape_2d_ref = res['full_data'].shape
             if 'ppm_axis' in res:
                 final_ppm_axes.append(res['ppm_axis'])
    if not valid_indices:
          return {"error": "No valid data loaded."}
    valid_concs = [concentrations[i] for i in valid_indices]
    list_2d = final_list_2d
    list_h = final_list_h
    list_n = final_list_n
    ppm_axes = final_ppm_axes
    ppm_axes_15n = final_ppm_axes_15n
    concentrations = valid_concs # Use filtered concs for analysis
    if not list_2d and not list_h:
         return {"error": "No data found"}
    if list_2d and shape_2d_ref is None and len(list_2d) > 0:
        if hasattr(list_2d[0], 'shape'):
             shape_2d_ref = list_2d[0].shape
    res_2d = run_pca_and_fit(
        list_2d, 
        concentrations, 
        protein_conc, 
        regime, 
        ppm_axes=ppm_axes if list_2d else None, 
        ppm_range=ppm_range, 
        exclude_ranges=[(4.6, 5.0)],
        no_fitting=no_fitting,
        ppm_15n=ppm_15n_ref, 
        shape_2d=shape_2d_ref
    ) if list_2d else None
    res_h = run_pca_and_fit(
        list_h, 
        concentrations, 
        protein_conc, 
        regime, 
        ppm_axes=ppm_axes, 
        ppm_range=ppm_range, 
        exclude_ranges=[(4.6, 5.0)],
        no_fitting=no_fitting
    ) if list_h else None
    res_n = run_pca_and_fit(
        list_n, 
        concentrations, 
        protein_conc, 
        regime, 
        ppm_axes=ppm_axes_15n if list_n else None,
        no_fitting=no_fitting
    ) if list_n else None
    return {
        "result_2d": res_2d,
        "result_h": res_h,
        "result_n": res_n,
        "concentrations": concentrations
    }
def run_residue_pca(residue_data, feature_names=None):
    data = np.array(residue_data)
    if data.ndim != 2 or data.shape[1] < 2:
        return {"error": "Insufficient data for PCA. Need at least 2 columns (ID + Feature)"}
    residue_nos = data[:, 0]
    features = data[:, 1:].astype(float) # All but first column
    n_features = features.shape[1]
    if not feature_names or len(feature_names) != n_features:
        feature_names = [f"Col {i+1}" for i in range(n_features)]
    mean = np.mean(features, axis=0)
    std = np.std(features, axis=0)
    std[std == 0] = 1.0
    features_scaled = (features - mean) / std
    n_components = min(3, n_features, features.shape[0])
    pca = PCA(n_components=n_components)
    scores = pca.fit_transform(features_scaled)
    loadings = pca.components_ 
    variance_ratio = pca.explained_variance_ratio_.tolist()
    results = []
    for i in range(n_components):
        results.append({
            "pc_index": i + 1,
            "scores": scores[:, i].tolist(),
            "loadings": {name: float(val) for name, val in zip(feature_names, loadings[i])},
            "success": False # "success" triggers Titration fit in frontend; FALSE triggers Residue view
        })
    return {
        "results": results,
        "variance_ratio": variance_ratio,
        "residue_nos": residue_nos.tolist(),
        "feature_names": feature_names
    }
def generate_hsqc_image(spectra_list, ppm_h, ppm_n, concentrations=None):
    try:
        if not spectra_list or not ppm_h or not ppm_n:
            return None
        fig, ax = plt.subplots(figsize=(10, 8)) # Larger figure
        if concentrations:
            pairs = sorted(zip(concentrations, spectra_list), key=lambda x: x[0])
            sorted_specs = [p[1] for p in pairs]
            sorted_concs = [p[0] for p in pairs]
        else:
            sorted_specs = spectra_list
            sorted_concs = list(range(len(sorted_specs)))
        n_spectra = len(sorted_specs)
        cmap = cm.nipy_spectral
        colors = cm.rainbow(np.linspace(0, 1, n_spectra)) 
        X_grid, Y_grid = np.meshgrid(ppm_h, ppm_n)
        for i, spec in enumerate(sorted_specs):
            s_max = np.max(spec)
            if s_max <= 0: continue
            base_lvl = s_max * 0.05 
            factor = 1.6 # 1.4^n spacing roughly
            levels = [base_lvl * (factor ** k) for k in range(8)]
            levels = [l for l in levels if l < s_max]
            if not levels: levels = [s_max * 0.5, s_max * 0.9] # Fallback
            ax.contour(X_grid, Y_grid, spec, levels=levels, colors=[colors[i]], linewidths=0.7)
        ax.set_xlabel('1H (ppm)', fontsize=12)
        ax.set_ylabel('15N (ppm)', fontsize=12)
        ax.invert_xaxis()
        ax.invert_yaxis()
        if sorted_concs:
            from matplotlib.lines import Line2D
            custom_lines = [Line2D([0], [0], color=colors[0], lw=2),
                            Line2D([0], [0], color=colors[-1], lw=2)]
            ax.legend(custom_lines, [f'{sorted_concs[0]} uM', f'{sorted_concs[-1]} uM'], loc='upper right')
        plt.tight_layout()
        buf = io.BytesIO()
        plt.savefig(buf, format='png', dpi=150, transparent=True) # Transparent background
        plt.close(fig)
        buf.seek(0)
        b64_str = base64.b64encode(buf.read()).decode('utf-8')
        return f"data:image/png;base64,{b64_str}"
    except Exception as e:
        print(f"Plotting Error: {e}")
        return None