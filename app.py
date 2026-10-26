import streamlit as st
import pandas as pd
from io import BytesIO

from rdkit import Chem
from rdkit.Chem import Draw, AllChem, Descriptors, Crippen, Lipinski, rdMolDescriptors
import py3Dmol
import streamlit.components.v1 as components
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.enums import TA_CENTER
from reportlab.lib.colors import HexColor
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle
from reportlab.lib import colors

st.set_page_config(page_title="Techie Lab | Cheminformatics", page_icon="🧪", layout="wide")

st.markdown("""
<style>
.block-container {padding-top:1.4rem; max-width:1400px}
.hero {padding:1.4rem 1.6rem;border-radius:16px;margin-bottom:1rem;
background:linear-gradient(120deg,#102a43,#17365d 55%,#167d8d);color:white}
.hero h1 {color:white;margin:0 0 .3rem 0}
.hero p {color:#e5f4f7;margin:0}
</style>
""", unsafe_allow_html=True)

st.markdown('<div class="hero"><h1>🧪 Techie Lab</h1><p>Cheminformatics Virtual Laboratory — explore molecular structures and properties.</p></div>', unsafe_allow_html=True)

EXAMPLES = {
    "Ethanol": "CCO",
    "Benzene": "c1ccccc1",
    "Phenol": "Oc1ccccc1",
    "Acetic acid": "CC(=O)O",
    "Aspirin": "CC(=O)Oc1ccccc1C(=O)O",
    "Caffeine": "Cn1c(=O)c2c(ncn2C)n(C)c1=O",
    "Paracetamol": "CC(=O)NC1=CC=C(O)C=C1",
    "Ibuprofen": "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O",
    "Glucose": "C([C@@H]1[C@H]([C@@H]([C@H](C(O1)O)O)O)O)O",
}

def properties(mol):
    return {
        "Molecular Formula": rdMolDescriptors.CalcMolFormula(mol),
        "Molecular Weight": round(Descriptors.MolWt(mol), 3),
        "Exact Molecular Weight": round(Descriptors.ExactMolWt(mol), 4),
        "LogP": round(Crippen.MolLogP(mol), 3),
        "TPSA (Å²)": round(rdMolDescriptors.CalcTPSA(mol), 3),
        "Hydrogen Bond Donors": Lipinski.NumHDonors(mol),
        "Hydrogen Bond Acceptors": Lipinski.NumHAcceptors(mol),
        "Rotatable Bonds": Lipinski.NumRotatableBonds(mol),
        "Ring Count": Lipinski.RingCount(mol),
        "Aromatic Rings": rdMolDescriptors.CalcNumAromaticRings(mol),
        "Aliphatic Rings": rdMolDescriptors.CalcNumAliphaticRings(mol),
        "Heavy Atoms": mol.GetNumHeavyAtoms(),
        "Bonds": mol.GetNumBonds(),
        "Fraction Csp3": round(rdMolDescriptors.CalcFractionCSP3(mol), 3),
        "Molar Refractivity": round(Crippen.MolMR(mol), 3),
    }

def make_pdf(smiles, props):
    buf = BytesIO()
    doc = SimpleDocTemplate(buf, pagesize=A4)
    styles = getSampleStyleSheet()
    title = ParagraphStyle("TitleCustom", parent=styles["Title"], textColor=HexColor("#17365D"),
                           alignment=TA_CENTER, spaceAfter=14)
    story = [Paragraph("Techie Lab — Molecular Report", title),
             Paragraph(f"<b>SMILES:</b> {smiles}", styles["BodyText"]), Spacer(1, 12)]
    rows = [["Property", "Value"]] + [[str(k), str(v)] for k, v in props.items()]
    table = Table(rows, colWidths=[250, 230])
    table.setStyle(TableStyle([
        ("BACKGROUND", (0,0), (-1,0), HexColor("#17365D")),
        ("TEXTCOLOR", (0,0), (-1,0), colors.white),
        ("GRID", (0,0), (-1,-1), 0.5, colors.lightgrey),
        ("ROWBACKGROUNDS", (0,1), (-1,-1), [colors.white, HexColor("#F2F6FA")]),
        ("PADDING", (0,0), (-1,-1), 7),
    ]))
    story.append(table)
    doc.build(story)
    buf.seek(0)
    return buf.getvalue()

left, right = st.columns([1, 1])
with left:
    chosen = st.selectbox("Choose an example molecule", list(EXAMPLES.keys()))
    default_smiles = EXAMPLES[chosen]
    smiles = st.text_input("Or enter a SMILES string", value=default_smiles,
                           help="SMILES is a text representation of a molecule.")
    show_h = st.checkbox("Include hydrogens in 3D view", value=True)
with right:
    st.info("Enter a valid SMILES string to view a 2D structure, calculate descriptors, inspect a 3D model, and export a report.")

mol = Chem.MolFromSmiles(smiles.strip()) if smiles.strip() else None
if mol is None:
    st.error("This SMILES string could not be parsed. Please check it and try again.")
    st.stop()

tab1, tab2, tab3 = st.tabs(["🧬 2D Structure & Properties", "🌐 3D Viewer", "📄 Export"])
with tab1:
    c1, c2 = st.columns([1, 1])
    with c1:
        st.subheader("2D molecular structure")
        try:
            st.image(Draw.MolToImage(mol, size=(650, 450)), use_container_width=True)
        except Exception as e:
            st.warning(f"2D drawing could not be generated: {e}")
    props = properties(mol)
    with c2:
        st.subheader("Molecular descriptors")
        st.dataframe(pd.DataFrame(list(props.items()), columns=["Property", "Value"]),
                     use_container_width=True, hide_index=True)
with tab2:
    st.subheader("Interactive 3D molecular viewer")
    mol3d = Chem.AddHs(Chem.Mol(mol)) if show_h else Chem.Mol(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = 42
    status = AllChem.EmbedMolecule(mol3d, params)
    if status == 0:
        try:
            AllChem.MMFFOptimizeMolecule(mol3d, maxIters=200)
        except Exception:
            pass
        viewer = py3Dmol.view(width=900, height=550)
        viewer.addModel(Chem.MolToMolBlock(mol3d), "mol")
        viewer.setStyle({"stick": {"radius": 0.16}, "sphere": {"scale": 0.25}})
        viewer.setBackgroundColor("#fbfdff")
        viewer.zoomTo()
        components.html(viewer._make_html(), height=570, scrolling=False)
    else:
        st.warning("RDKit could not generate a 3D conformation for this molecule.")
with tab3:
    st.subheader("Download your results")
    csv = pd.DataFrame([{"SMILES": smiles, **props}]).to_csv(index=False).encode("utf-8")
    st.download_button("⬇️ Download descriptors (CSV)", csv, "molecular_descriptors.csv", "text/csv")
    pdf = make_pdf(smiles, props)
    st.download_button("⬇️ Download molecular report (PDF)", pdf, "techie_lab_molecular_report.pdf", "application/pdf")

st.caption("For educational use only. Molecular descriptors are computational estimates and should be independently verified for research or clinical use.")
