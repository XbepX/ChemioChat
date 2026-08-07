import streamlit as st
from chemio_tools.io_handlers import carregar_mol_de_texto


def renderizar_secao_entradas() -> str:
    """Renderiza a barra de entrada SMILES, histórico e o carregador de arquivos."""
    SMILES_INICIAL = "O=C[C@H](O)[C@@H](O)[C@H](O)[C@H](O)CO"

    if "historico_smiles" not in st.session_state:
        st.session_state.historico_smiles = [
            SMILES_INICIAL,
            "CC(=O)O",
            "CCO",
            "c1ccccc1"
        ]

    col_input, col_hist = st.columns([3, 1])

    with col_hist:
        smiles_selecionado = st.selectbox(
            "📜 Histórico / Exemplos:",
            options=st.session_state.historico_smiles,
            index=0
        )

    with col_input:
        smiles_input = st.text_input(
            "Digite o SMILES da molécula:",
            value=smiles_selecionado
        )

    with st.expander("📂 Ou carregar arquivo de estrutura (.mol, .sdf, .pdb)"):
        arquivo_enviado = st.file_uploader(
            "Arraste o arquivo exportado (ChemSketch, Avogadro, etc.):",
            type=["mol", "sdf", "pdb"]
        )
        if arquivo_enviado is not None:
            conteudo_texto = arquivo_enviado.getvalue().decode("utf-8")
            extensao = arquivo_enviado.name.split(".")[-1].upper()

            smiles_extraido = carregar_mol_de_texto(conteudo_texto, formato=extensao)
            if smiles_extraido:
                smiles_input = smiles_extraido
                st.success(f"Estrutura importada do arquivo: `{arquivo_enviado.name}`")
            else:
                st.error("Não foi possível ler a estrutura deste arquivo.")

    if smiles_input and smiles_input not in st.session_state.historico_smiles:
        st.session_state.historico_smiles.insert(0, smiles_input)

    return smiles_input