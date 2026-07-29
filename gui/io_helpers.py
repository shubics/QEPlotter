"""Shared helpers for the GUI: temp-file uploads, Fermi auto-detect, channel scan."""
import glob
import os
import re
import tempfile

import streamlit as st


def ensure_temp_dir():
    """Create (once per session) and return a temp dir for uploaded files."""
    if "temp_dir" not in st.session_state:
        st.session_state.temp_dir = tempfile.mkdtemp()
    return st.session_state.temp_dir


def save_file(uploaded_file, subdir=None):
    """Persist a Streamlit UploadedFile to the session temp dir; return its path."""
    if uploaded_file is None:
        return None
    target_dir = ensure_temp_dir()
    if subdir:
        target_dir = os.path.join(target_dir, subdir)
        os.makedirs(target_dir, exist_ok=True)
    # Browser-provided names are display metadata, not trusted paths.
    safe_name = os.path.basename(str(uploaded_file.name).replace("\\", "/"))
    if not safe_name:
        raise ValueError("The uploaded file has no usable filename.")
    file_path = os.path.join(target_dir, safe_name)
    with open(file_path, "wb") as f:
        f.write(uploaded_file.getbuffer())
    return file_path


@st.cache_data
def get_fermi_from_scf(scf_path):
    """Parse scf.out to find the Fermi energy (metals, insulators, semiconductors)."""
    try:
        with open(scf_path, "r", errors="ignore") as f:
            content = f.read()

        # Metal: "the Fermi energy is     1.2345 eV"
        m = re.search(r"the Fermi energy is\s+([-+]?\d*\.?\d+)\s+eV", content)
        if m:
            return float(m.group(1))

        # Insulator/semiconductor: "highest occupied, lowest unoccupied level (ev): HOMO LUMO"
        m2 = re.search(
            r"highest occupied, lowest unoccupied level \(ev\):\s+"
            r"([-+]?\d*\.?\d+)\s+([-+]?\d*\.?\d+)",
            content,
        )
        if m2:
            homo = float(m2.group(1))
            lumo = float(m2.group(2))
            return 0.5 * (homo + lumo)  # mid-gap estimate

        # "highest occupied level (ev):     1.2345"
        m3 = re.search(r"highest occupied level\s*\(ev\):\s+([-+]?\d*\.?\d+)", content)
        if m3:
            return float(m3.group(1))

        return None
    except Exception:
        return None


def get_available_channels(pdos_dir):
    """Scan a PDOS dir for available atoms / elements / orbitals / element-orbitals."""
    if not pdos_dir or not os.path.exists(pdos_dir):
        return [], [], [], []
    files = glob.glob(os.path.join(pdos_dir, "*pdos*"))

    atoms, elements, orbitals, elem_orbs = set(), set(), set(), set()
    pattern = re.compile(r"atm#(\d+)\(([A-Za-z]+)\)_wfc#\d+\(([a-zA-Z0-9_.]+)\)")
    for f in files:
        m = pattern.search(os.path.basename(f))
        if m:
            num, elem, orb = m.groups()
            base_orb = orb.split("_")[0]
            atoms.add(f"{elem}{num}")
            elements.add(elem)
            orbitals.add(base_orb)
            elem_orbs.add(f"{elem}-{base_orb}")

    return (
        sorted(atoms),
        sorted(elements),
        sorted(orbitals),
        sorted(elem_orbs),
    )
