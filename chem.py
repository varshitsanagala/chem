import streamlit as st
from rdkit import Chem
from rdkit.Chem import Draw, AllChem, rdMolDescriptors
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit import Geometry
import base64
from io import BytesIO
from PIL import Image
import pandas as pd

import sys
import streamlit as st

st.write("Python path:", sys.executable)

try:
    from rdkit import Chem
    st.write("RDKit imported successfully ✅")
    test = Chem.MolFromSmiles("CCO")
    st.write("Test molecule:", test)
except Exception as e:
    st.error(f"RDKit failed ❌: {e}")


# ── Page config ──────────────────────────────────────────────────────────────
st.set_page_config(
    page_title="Bleomycin · Stereochemistry Explorer",
    page_icon="🧬",
    layout="wide",
)

# ── Custom CSS ────────────────────────────────────────────────────────────────
st.markdown("""
<style>
@import url('https://fonts.googleapis.com/css2?family=Space+Mono:ital,wght@0,400;0,700;1,400&family=Syne:wght@400;600;800&display=swap');

:root {
    --bg:        #0a0d14;
    --surface:   #111520;
    --border:    #1e2535;
    --accent:    #00e5c0;
    --accent2:   #ff6b6b;
    --accent3:   #ffd93d;
    --text:      #e8eaf6;
    --muted:     #6b7280;
    --r-col:     #ff6b6b;
    --s-col:     #00e5c0;
}

html, body, [class*="css"] {
    background-color: var(--bg) !important;
    color: var(--text);
    font-family: 'Syne', sans-serif;
}

/* Hide Streamlit chrome */
#MainMenu, footer, header { visibility: hidden; }
.block-container { padding: 2rem 3rem; max-width: 1400px; }

/* Hero header */
.hero {
    border-bottom: 1px solid var(--border);
    padding-bottom: 1.5rem;
    margin-bottom: 2rem;
}
.hero h1 {
    font-family: 'Syne', sans-serif;
    font-weight: 800;
    font-size: 3rem;
    letter-spacing: -0.02em;
    background: linear-gradient(135deg, var(--accent) 0%, #7c3aed 100%);
    -webkit-background-clip: text;
    -webkit-text-fill-color: transparent;
    background-clip: text;
    margin: 0;
}
.hero .sub {
    font-family: 'Space Mono', monospace;
    font-size: 0.78rem;
    color: var(--muted);
    letter-spacing: 0.12em;
    text-transform: uppercase;
    margin-top: 0.4rem;
}

/* Cards */
.card {
    background: var(--surface);
    border: 1px solid var(--border);
    border-radius: 12px;
    padding: 1.5rem;
    margin-bottom: 1.2rem;
}
.card-title {
    font-family: 'Space Mono', monospace;
    font-size: 0.7rem;
    letter-spacing: 0.15em;
    text-transform: uppercase;
    color: var(--accent);
    margin-bottom: 0.8rem;
}

/* Stat pills */
.stat-row { display: flex; gap: 1rem; flex-wrap: wrap; margin-bottom: 1rem; }
.stat-pill {
    background: #0d1117;
    border: 1px solid var(--border);
    border-radius: 8px;
    padding: 0.6rem 1rem;
    font-family: 'Space Mono', monospace;
}
.stat-pill .val {
    font-size: 1.6rem;
    font-weight: 700;
    line-height: 1;
}
.stat-pill .lbl {
    font-size: 0.65rem;
    color: var(--muted);
    text-transform: uppercase;
    letter-spacing: 0.1em;
}
.pill-r .val { color: var(--r-col); }
.pill-s .val { color: var(--s-col); }
.pill-t .val { color: var(--accent3); }

/* Table */
.chiral-table { width: 100%; border-collapse: collapse; font-family: 'Space Mono', monospace; font-size: 0.8rem; }
.chiral-table th {
    text-align: left; padding: 0.5rem 0.8rem;
    border-bottom: 1px solid var(--border);
    color: var(--muted); text-transform: uppercase; letter-spacing: 0.1em; font-size: 0.65rem;
}
.chiral-table td { padding: 0.45rem 0.8rem; border-bottom: 1px solid #1a1f2e; }
.chiral-table tr:last-child td { border-bottom: none; }
.badge {
    display: inline-block; border-radius: 4px; padding: 2px 8px;
    font-weight: 700; font-size: 0.75rem; letter-spacing: 0.05em;
}
.badge-R { background: rgba(255,107,107,0.15); color: var(--r-col); border: 1px solid rgba(255,107,107,0.3); }
.badge-S { background: rgba(0,229,192,0.12); color: var(--s-col); border: 1px solid rgba(0,229,192,0.3); }

/* Legend */
.legend { display: flex; gap: 1.5rem; align-items: center; font-family: 'Space Mono', monospace; font-size: 0.75rem; }
.leg-dot { width: 12px; height: 12px; border-radius: 50%; display: inline-block; margin-right: 6px; }

/* Info box */
.info-box {
    background: rgba(0,229,192,0.04);
    border-left: 3px solid var(--accent);
    padding: 0.8rem 1rem;
    border-radius: 0 8px 8px 0;
    font-size: 0.85rem;
    line-height: 1.6;
    color: #b0bec5;
}

/* Molecule container */
.mol-wrap {
    background: #f5f5f0;
    border-radius: 10px;
    overflow: hidden;
    border: 1px solid var(--border);
}

div[data-testid="stSelectbox"] label,
div[data-testid="stCheckbox"] label {
    font-family: 'Space Mono', monospace !important;
    font-size: 0.75rem !important;
    text-transform: uppercase;
    letter-spacing: 0.08em;
    color: var(--muted) !important;
}
</style>
""", unsafe_allow_html=True)

# ── Bleomycin A2 SMILES (canonical) ──────────────────────────────────────────
BLEOMYCIN_SMILES = (
    "[NH3+][C@@H](CCc1c[nH]cn1)C(=O)N[C@@H](C[C@H](O)[C@@H](NC(=O)"
    "[C@@H](CC(N)=O)NC(=O)[C@H](C)NC(=O)[C@@H](NC(=O)/C(=C/C)NC(=O)"
    "[C@@H](NC(=O)[C@@H]1CCCN1C(=O)[C@@H](NC(=O)[C@H](O)[C@@H](NC(=O)"
    "c1nc2c(s1)cc(cc2)C(=O)N[C@@H](CO)[C@H](O)c1ccc(OC)cc1)"
    "c1ccc(O)cc1)C[S+](C)CC[C@@H](NC(=O)CC)C(N)=O)CC(=O)N)"
    "C(N)=O)C(=O)N[C@@H]1C[C@H](O)[C@@H](n2cnc3c(N)ncnc23)O1"
)

# Fallback to a simpler known SMILES if parsing fails
FALLBACK_SMILES = (
    "CC1=C(N=C(N=C1NC2=CC=CC=C2)NC3=CC=CC=C3)NC4=CC=CC=C4"  # just for fallback safety
)

@st.cache_data
def load_molecule():
    mol = Chem.MolFromSmiles(BLEOMYCIN_SMILES, sanitize=False)

    if mol is None:
        st.warning("⚠️ Using fallback molecule (Bleomycin too complex for strict parsing)")
        mol = Chem.MolFromSmiles(FALLBACK_SMILES)
    else:
        try:
            Chem.SanitizeMol(mol)
        except Exception as e:
            st.warning(f"Sanitization warning: {e}")

    AllChem.Compute2DCoords(mol)
    return mol  

@st.cache_data
def get_chiral_centers(mol_smiles):
    mol = Chem.MolFromSmiles(mol_smiles, sanitize=False)

    if mol is None:
        return []

    try:
        Chem.SanitizeMol(mol)
    except:
        pass

    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
    return centers

def draw_molecule(mol, highlight_atoms=None, highlight_colors=None,
                  width=900, height=600, show_atom_idx=False):
    drawer = rdMolDraw2D.MolDraw2DSVG(width, height)
    opts = drawer.drawOptions()
    opts.addStereoAnnotation = True
    opts.addAtomIndices = show_atom_idx
    opts.bondLineWidth = 2.0
    opts.addStereoAnnotation = True
    opts.addAtomIndices = show_atom_idx
    opts.bondLineWidth = 2.0
    opts.padding = 0.08

    # White background
    drawer.SetBackgroundColour((1, 1, 1, 1))

    if highlight_atoms and highlight_colors:
        atom_cols = {a: highlight_colors.get(a, (0.9, 0.9, 0.9)) for a in highlight_atoms}
        bond_cols = {}
        drawer.DrawMolecule(mol,
                            highlightAtoms=highlight_atoms,
                            highlightAtomColors=atom_cols,
                            highlightBonds=[],
                            highlightBondColors=bond_cols)
    else:
        drawer.DrawMolecule(mol)

    drawer.FinishDrawing()
    svg = drawer.GetDrawingText()
    return svg

def svg_to_img_tag(svg: str) -> str:
    b64 = base64.b64encode(svg.encode()).decode()
    return f'<img src="data:image/svg+xml;base64,{b64}" style="width:100%; border-radius:10px;" />'

# ── Load data ─────────────────────────────────────────────────────────────────
mol = load_molecule()
chiral_centers = get_chiral_centers(BLEOMYCIN_SMILES) if mol else []

r_atoms = [idx for idx, cfg in chiral_centers if cfg == 'R']
s_atoms = [idx for idx, cfg in chiral_centers if cfg == 'S']
unassigned = [idx for idx, cfg in chiral_centers if cfg == '?']

# ── Hero ──────────────────────────────────────────────────────────────────────
st.markdown("""
<div class="hero">
  <h1>Bleomycin A₂</h1>
  <div class="sub">Stereochemistry Explorer &nbsp;·&nbsp; Chiral Centre Mapper &nbsp;·&nbsp; R/S Configuration</div>
</div>
""", unsafe_allow_html=True)

# ── Layout: 2 columns ─────────────────────────────────────────────────────────
col_mol, col_info = st.columns([3, 2], gap="large")

with col_info:
    # Stats
    st.markdown(f"""
    <div class="card">
      <div class="card-title">Stereochemistry Summary</div>
      <div class="stat-row">
        <div class="stat-pill pill-t">
          <div class="val">{len(chiral_centers)}</div>
          <div class="lbl">Total Chiral Centres</div>
        </div>
        <div class="stat-pill pill-r">
          <div class="val">{len(r_atoms)}</div>
          <div class="lbl">R Configuration</div>
        </div>
        <div class="stat-pill pill-s">
          <div class="val">{len(s_atoms)}</div>
          <div class="lbl">S Configuration</div>
        </div>
      </div>
      <div class="legend">
        <span><span class="leg-dot" style="background:#ff6b6b"></span>R carbon</span>
        <span><span class="leg-dot" style="background:#00e5c0"></span>S carbon</span>
      </div>
    </div>
    """, unsafe_allow_html=True)

    # Controls
    st.markdown('<div class="card"><div class="card-title">Display Options</div>', unsafe_allow_html=True)
    highlight_mode = st.selectbox(
        "Highlight",
        ["All chiral centres", "R centres only", "S centres only", "None"],
        index=0,
    )
    show_idx = st.checkbox("Show atom indices", value=False)
    st.markdown('</div>', unsafe_allow_html=True)

    # Chiral table
    if chiral_centers:
        rows_html = ""
        for idx, cfg in sorted(chiral_centers, key=lambda x: x[0]):
            badge = f'<span class="badge badge-{cfg}">{cfg}</span>' if cfg in ('R','S') else f'<span class="badge" style="background:#2a2a2a;color:#888">{cfg}</span>'
            atom = mol.GetAtomWithIdx(idx)
            sym = atom.GetSymbol()
            rows_html += f"<tr><td style='color:#6b7280'>{idx}</td><td>{sym}</td><td>{badge}</td></tr>"

        st.markdown(f"""
        <div class="card">
          <div class="card-title">Chiral Centres — All {len(chiral_centers)}</div>
          <table class="chiral-table">
            <thead><tr><th>Atom #</th><th>Element</th><th>Config</th></tr></thead>
            <tbody>{rows_html}</tbody>
          </table>
        </div>
        """, unsafe_allow_html=True)

    # Molecule info
    if mol:
        mw = rdMolDescriptors.CalcExactMolWt(mol)
        formula = rdMolDescriptors.CalcMolFormula(mol)
        hbd = rdMolDescriptors.CalcNumHBD(mol)
        hba = rdMolDescriptors.CalcNumHBA(mol)
        rings = rdMolDescriptors.CalcNumRings(mol)

        st.markdown(f"""
        <div class="card">
          <div class="card-title">Molecular Properties</div>
          <table class="chiral-table">
            <tbody>
              <tr><td style="color:var(--muted)">Formula</td><td style="font-family:'Space Mono',monospace">{formula}</td></tr>
              <tr><td style="color:var(--muted)">Exact MW</td><td>{mw:.2f} Da</td></tr>
              <tr><td style="color:var(--muted)">H-bond donors</td><td>{hbd}</td></tr>
              <tr><td style="color:var(--muted)">H-bond acceptors</td><td>{hba}</td></tr>
              <tr><td style="color:var(--muted)">Rings</td><td>{rings}</td></tr>
            </tbody>
          </table>
        </div>
        """, unsafe_allow_html=True)

    # Info
    st.markdown("""
    <div class="info-box">
        <strong style="color:#e8eaf6">Bleomycin A₂</strong> is a glycopeptide antibiotic and antineoplastic agent
        produced by <em>Streptomyces verticillus</em>. Its multiple stereocentres are critical
        for DNA strand scission activity and metal-ion chelation.
        <br><br>Chiral centres are highlighted using CIP (Cahn–Ingold–Prelog) priority rules
        as implemented in RDKit.
    </div>
    """, unsafe_allow_html=True)

with col_mol:
    st.markdown('<div class="card"><div class="card-title">2D Structure Viewer</div>', unsafe_allow_html=True)

    if mol:
        # Build highlight sets based on selection
        if highlight_mode == "All chiral centres":
            h_atoms = r_atoms + s_atoms + unassigned
            h_colors = {a: (1.0, 0.42, 0.42) for a in r_atoms}
            h_colors.update({a: (0.0, 0.9, 0.75) for a in s_atoms})
            h_colors.update({a: (1.0, 0.85, 0.24) for a in unassigned})
        elif highlight_mode == "R centres only":
            h_atoms = r_atoms
            h_colors = {a: (1.0, 0.42, 0.42) for a in r_atoms}
        elif highlight_mode == "S centres only":
            h_atoms = s_atoms
            h_colors = {a: (0.0, 0.9, 0.75) for a in s_atoms}
        else:
            h_atoms = None
            h_colors = None

        svg = draw_molecule(mol, highlight_atoms=h_atoms, highlight_colors=h_colors,
                            width=900, height=580, show_atom_idx=show_idx)
        st.markdown(
            f'<div class="mol-wrap">{svg_to_img_tag(svg)}</div>',
            unsafe_allow_html=True
        )
    else:
        st.error("⚠️ Could not parse Bleomycin SMILES. Please check RDKit installation.")

    st.markdown('</div>', unsafe_allow_html=True)

    # SMILES display
    with st.expander("📋 SMILES String"):
        st.code(BLEOMYCIN_SMILES, language=None)

# ── R/S breakdown bar ─────────────────────────────────────────────────────────
if chiral_centers:
    st.markdown("---")
    st.markdown('<div class="card-title" style="font-family:\'Space Mono\',monospace;font-size:0.7rem;letter-spacing:0.15em;text-transform:uppercase;color:#00e5c0;margin-bottom:0.6rem">Configuration Distribution</div>', unsafe_allow_html=True)

    df = pd.DataFrame(chiral_centers, columns=["Atom Index", "Configuration"])
    cfg_counts = df["Configuration"].value_counts().reset_index()
    cfg_counts.columns = ["Configuration", "Count"]

    c1, c2, c3 = st.columns(3)
    total = len(chiral_centers)

    for col, cfg, color in zip([c1, c2, c3],
                                ["R", "S", "?"],
                                ["#ff6b6b", "#00e5c0", "#ffd93d"]):
        n = len([x for x in chiral_centers if x[1] == cfg])
        pct = (n / total * 100) if total else 0
        with col:
            st.markdown(f"""
            <div style="background:#111520;border:1px solid #1e2535;border-radius:10px;padding:1rem;text-align:center">
              <div style="font-size:2.2rem;font-weight:800;color:{color};font-family:'Syne',sans-serif">{n}</div>
              <div style="font-family:'Space Mono',monospace;font-size:0.65rem;color:#6b7280;text-transform:uppercase;letter-spacing:0.1em">{cfg} Configuration</div>
              <div style="margin-top:0.6rem;background:#0a0d14;border-radius:4px;height:6px;overflow:hidden">
                <div style="background:{color};width:{pct:.0f}%;height:100%;border-radius:4px;transition:width 0.4s ease"></div>
              </div>
              <div style="font-size:0.7rem;color:#6b7280;margin-top:0.3rem;font-family:'Space Mono',monospace">{pct:.0f}%</div>
            </div>
            """, unsafe_allow_html=True)