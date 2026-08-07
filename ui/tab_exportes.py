import streamlit as st
from chemio_tools import exportar_bloco_3d


def renderizar_aba_exportacao(smiles_input: str):
    """Renderiza a Aba 3: Download de arquivos 3D (.sdf e .pdb)."""
    st.subheader("📁 Download de Coordenadas 3D para Editores")
    st.caption("Baixe a molécula otimizada para abrir em softwares como Avogadro, PyMOL ou ChemSketch.")

    c_sdf, c_pdb, _ = st.columns([1.2, 1.2, 2])

    with c_sdf:
        sdf_data = exportar_bloco_3d(smiles_input, formato="SDF")
        if sdf_data:
            st.download_button(
                label="⬇️ Baixar Bloco SDF (.sdf)",
                data=sdf_data,
                file_name="molecula_3d.sdf",
                mime="chemical/x-mdl-sdfile"
            )

    with c_pdb:
        pdb_data = exportar_bloco_3d(smiles_input, formato="PDB")
        if pdb_data:
            st.download_button(
                label="⬇️ Baixar Arquivo PDB (.pdb)",
                data=pdb_data,
                file_name="molecula_3d.pdb",
                mime="chemical/x-pdb"
            )