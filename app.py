import streamlit as st
from chemio_tools import validar_smiles
from ui import (
    renderizar_secao_entradas,
    renderizar_aba_visualizacao,
    renderizar_aba_propriedades,
    renderizar_aba_exportacao
)

# Configuração da página
st.set_page_config(page_title="ChemioChat", layout="wide")
st.title("🧪 ChemioChat - Laboratório de Testes")

# 1. Seção Superior (Entrada de SMILES, Histórico e Upload)
smiles_input = renderizar_secao_entradas()

if not smiles_input.strip():
    st.stop()

if not validar_smiles(smiles_input):
    st.error("SMILES inválido! Verifique a sintaxe química digitada.")
    st.stop()

st.success("SMILES válido! A estrutura foi reconhecida pelo RDKit.")
st.write(f"SMILES capturado: `{smiles_input}`")

# Reset de estados ao alterar o SMILES
if "ultimo_smiles" not in st.session_state or st.session_state.ultimo_smiles != smiles_input:
    st.session_state.ultimo_smiles = smiles_input
    st.session_state.exibir_estereo = False
    st.session_state.exibir_isomeros = False
    st.session_state.grupo_ativo = "Nenhum"

# 2. Estrutura de Abas
tab_vis, tab_props, tab_export = st.tabs([
    "👁️ Visualização & Análise Interativa",
    "📊 Propriedades & Composição",
    "💾 Exportar Estrutura"
])

# 3. Chamada dos Módulos da Pasta UI
with tab_vis:
    renderizar_aba_visualizacao(smiles_input)

with tab_props:
    renderizar_aba_propriedades(smiles_input)

with tab_export:
    renderizar_aba_exportacao(smiles_input)