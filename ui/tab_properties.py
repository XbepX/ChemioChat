import streamlit as st
from chemio_tools import obter_propriedades_basicas


def renderizar_aba_propriedades(smiles_input: str):
    """Renderiza a Aba 2: Propriedades físico-químicas e contagem atômica."""
    dados = obter_propriedades_basicas(smiles_input)

    if dados:
        st.subheader("Métricas Físico-Químicas")
        col1, col2, col3, col4 = st.columns(4)
        col1.metric("Fórmula Molecular", dados["formula"])
        col2.metric("Massa Molecular", f"{dados['massa_molecular']} g/mol")
        col3.metric("LogP (Lipofilicidade)", dados["logp"])
        col4.metric("TPSA", f"{dados['tpsa']} Å²")

        st.divider()
        st.subheader("Contagem de Átomos (Com H Explícitos)")

        elementos = dados["contagem_elementos"]
        cols_elem = st.columns(len(elementos))
        for col, (simbolo, quantidade) in zip(cols_elem, elementos.items()):
            col.metric(label=f"Átomos de {simbolo}", value=quantidade)