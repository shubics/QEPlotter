"""
DOS and PDOS plotting functions.
Extracted verbatim from qep.py (plot_dos, plot_pdos_dir).
"""
import os
import glob
import re
import matplotlib.pyplot as plt
import numpy as np


def plot_dos(dos_file, fermi_level=None, shift_fermi=False, y_range=None, x_range=None, dpi=None,
        save_dir="saved", savefig=None, vertical=False):
    """
    Plot the total Density of States (DOS) from a QE DOS file.
    """
    try:
        data = np.loadtxt(dos_file, comments='#')
    except Exception as e:
        try:
             data = np.loadtxt(dos_file, skiprows=1)
        except Exception:
             raise ValueError(f"Could not read DOS file {dos_file}. Check format. Error: {e}")

    if data.ndim != 2 or data.shape[1] < 2:
        raise ValueError(f"Unexpected DOS file format: {dos_file}")
    E = data[:, 0]
    DOS = data[:, 1]
    if shift_fermi and fermi_level is not None:
        E = E - fermi_level
    if dpi is not None:
        plt.figure(figsize=(6,6), dpi=dpi)
    else:
        plt.figure(figsize=(6,6))

    if vertical:
        # DOS on X, Energy on Y
        plt.plot(DOS, E, 'k-', lw=1, label='Total DOS')
        if fermi_level is not None:
            y0 = 0.0 if shift_fermi else fermi_level
            plt.axhline(y0, color='r', ls='--', lw=1.2, label=f'Fermi = {fermi_level:.2f} eV')
        
        ylabel = 'E - E_F (eV)' if (shift_fermi and fermi_level is not None) else 'Energy (eV)'
        plt.ylabel(ylabel)
        plt.xlabel('DOS')
        if y_range:
            plt.ylim(y_range)
        if x_range:
            plt.xlim(x_range)
    else:
        # Energy on X, DOS on Y
        plt.plot(E, DOS, 'k-', lw=1, label='Total DOS')
        if fermi_level is not None:
            x0 = 0.0 if shift_fermi else fermi_level
            plt.axvline(x0, color='r', ls='--', lw=1.2, label=f'Fermi = {fermi_level:.2f} eV')
        
        xlabel = 'E - E_F (eV)' if (shift_fermi and fermi_level is not None) else 'Energy (eV)'
        plt.xlabel(xlabel)
        plt.ylabel('DOS')
        if y_range:
            plt.ylim(y_range)
        if x_range:
            plt.xlim(x_range)
    plt.title('Total DOS')
    plt.grid(True, ls='--', alpha=0.4)
    if fermi_level is not None:
        plt.legend()
    plt.tight_layout()
    if savefig:
        if not os.path.exists(save_dir):
            os.makedirs(save_dir)
        out = os.path.join(save_dir, os.path.basename(savefig))
        plt.savefig(out, dpi=dpi or plt.rcParams['figure.dpi'])
        print(f"Saved figure to {out}")

    plt.show()

def plot_pdos_dir(pdos_dir, fermi_level=None,
                  shift_fermi=False, y_range=None, dpi=None,pdos_mode='atomic',
                  save_dir="saved", savefig=None):
    """
    Plot projected Density of States (PDOS) from a set of QE projwfc/pdos files.
    """
    pat = re.compile(r'atm#\d+\(([A-Za-z]+)\)_wfc#\d+\(([spdfgpxyz]+)(?:_j[0-9.]+)?\)')
    fallback_pat = re.compile(r'atm#\d+\(([A-Za-z]+)\)_wfc#\d+\(([spdfgpxyz]+)\)')

    files = glob.glob(os.path.join(pdos_dir, '*pdos*'))
    if not files:
        raise FileNotFoundError(f"No PDOS files found in {pdos_dir}")
    grouped = {}
    E = None
    for fn in files:
        base = os.path.basename(fn)
        m = pat.search(base)
        if not m:
            m = fallback_pat.search(base)
        
        if not m:
            continue
            
        elem, orb = m.group(1), m.group(2)
        
        if pdos_mode == 'atomic':
            key = elem
        elif pdos_mode == 'orbital':
            key = orb
        elif pdos_mode == 'element_orbital':
            key = f"{elem}-{orb}"
        else:
            raise ValueError(f"Unknown pdos_mode: {pdos_mode}")
            
        data = np.loadtxt(fn, comments='#')
        if data.ndim != 2 or data.shape[1] < 2:
            continue
            
        if np.all(data[:,0] == np.arange(data.shape[0])):
             print(f"Warning: File {base} looks like a fatband file (Col 0 is index), but we expect Energy. Skipping to avoid bad plot.")
             continue

        if E is None:
            E = data[:, 0].copy()
            if shift_fermi and fermi_level is not None:
                E = E - fermi_level
                
        if data.shape[1] >= 2:
             pd = data[:, 1]
        else:
             pd = data[:, -1]
             
        grouped.setdefault(key, np.zeros_like(pd))
        grouped[key] += pd
        
    if not grouped:
        raise RuntimeError("No PDOS channels matched; check filenames and pdos_mode")
    if dpi is not None:
        plt.figure(figsize=(6,6), dpi=dpi)
    else:
        plt.figure(figsize=(6,6))
    for k, pd in sorted(grouped.items()):
        plt.plot(E, pd, lw=1, label=k)
    if fermi_level is not None:
        x0 = 0.0 if shift_fermi else fermi_level
        plt.axvline(x0, color='r', ls='--', lw=1.2, label=f'Fermi = {fermi_level:.2f} eV')
    
    xlabel = 'E - E_F (eV)' if (shift_fermi and fermi_level is not None) else 'Energy (eV)'
    plt.xlabel(xlabel)

    plt.ylabel('Projected DOS')
    if y_range:
        plt.ylim(y_range)
    plt.title(f'Projected DOS ({pdos_mode})')
    plt.grid(True, ls='--', alpha=0.4)
    plt.legend(fontsize='small', ncol=2)
    plt.tight_layout()
    if savefig:
        if not os.path.exists(save_dir):
            os.makedirs(save_dir)
        out = os.path.join(save_dir, os.path.basename(savefig))
        plt.savefig(out, dpi=dpi or plt.rcParams['figure.dpi'])
        print(f"Saved figure to {out}")

    plt.show()
