import streamlit as st
import streamlit.components.v1 as components
from chemio_tools import (
    gerar_imagem_2d,
    gerar_html_3d,
    identificar_grupos_funcionais,
    identificar_centros_quirais,
    gerar_estereoisomeros
)


def renderizar_aba_visualizacao(smiles_input: str):
    """Renderiza a Aba 1: Painel visual 2D/3D, seleção de grupos e estereoquímica."""
    grupos_funcionais = identificar_grupos_funcionais(smiles_input)

    # Controle de sessão dos painéis inferiores
    if "exibir_estereo" not in st.session_state:
        st.session_state.exibir_estereo = False
    if "exibir_isomeros" not in st.session_state:
        st.session_state.exibir_isomeros = False
    if "grupo_ativo" not in st.session_state:
        st.session_state.grupo_ativo = "Nenhum"

    col_render, col_controles = st.columns([3, 1])

    # PAINEL LATERAL DE CONTROLES
    with col_controles:
        st.subheader("⚙️ Opções")
        modo_vis = st.radio("Dimensão:", ("2D", "3D"), horizontal=True)
        st.divider()

        indices_grupo_selecionado = []

        if modo_vis == "2D":
            with st.expander("🏷️ Átomos e Rótulos", expanded=True):
                opt_mostrar_h = st.checkbox("Hidrogênios explícitos", value=False)
                opt_destacar_q = st.checkbox("Destacar quirais (* e R/S)", value=False)
                opt_mostrar_idx = st.checkbox("Índices dos carbonos (C1, C2...)", value=True)
                opt_metilas = st.checkbox("Exibir metilas (-CH₃)", value=False)

            with st.expander("🧪 Grupos Funcionais Mapeados", expanded=True):
                if grupos_funcionais:
                    st.caption("Clique no grupo para destacar no desenho:")
                    
                    if st.button("⚪ Limpar Destaque", key="btn_limpar_grupo"):
                        st.session_state.grupo_ativo = "Nenhum"

                    for g in grupos_funcionais:
                        nome_g = g["grupo"]
                        rotulo_btn = f"📌 {nome_g}" if st.session_state.grupo_ativo == nome_g else f"🧪 {nome_g}"
                        
                        if st.button(rotulo_btn, key=f"btn_grupo_{nome_g}"):
                            st.session_state.grupo_ativo = nome_g

                    if st.session_state.grupo_ativo != "Nenhum":
                        for g in grupos_funcionais:
                            if g["grupo"] == st.session_state.grupo_ativo:
                                indices_grupo_selecionado = g["indices"]
                else:
                    st.caption("Nenhum grupo funcional mapeado.")

            with st.expander("🎨 Estilo e Traço", expanded=False):
                opt_espessura = st.slider("Espessura da linha", min_value=1, max_value=5, value=2)
                opt_fonte = st.slider("Tamanho do texto", min_value=8, max_value=20, value=12)
                opt_fundo_escuro_2d = st.checkbox("Fundo escuro", value=False)

        else:
            with st.expander("🧊 Geometria e Estilo", expanded=True):
                opt_estilo_3d = st.selectbox(
                    "Estilo do Modelo",
                    ["Ball & Stick", "Wireframe / Line", "Stick (Bastonetes)", "CPK / Esferas"]
                )
                opt_espessura_3d = st.slider("Espessura das ligações/linhas", 1, 10, 5)
                opt_esquema_cor_3d = st.selectbox(
                    "Esquema de Cores",
                    ["Monocromático Escuro", "Elemento (Jmol)", "Carbono Preto (CPK)", "Monocromático Claro"]
                )
                opt_mostrar_h_3d = st.checkbox("Mostrar Hidrogênios", value=True)
                opt_destacar_q_3d = st.checkbox("Destacar quirais (Malha 3D)", value=False)
                opt_raio_malha = st.slider("Raio da malha quiral (Å)", 0.8, 2.5, 1.2, step=0.1)

            with st.expander("🌐 Superfície e Efeitos", expanded=False):
                opt_superficie = st.selectbox(
                    "Superfície Eletrônica",
                    ["Nenhuma", "Van der Waals (VDW)", "Superfície Acessível ao Solvente (SAS)"]
                )
                opt_coloracao = st.selectbox(
                    "Coloração da Superfície",
                    ["Mapa Eletrostático (MEP)", "Monocromático (Azul Claro)"]
                )
                opt_opacidade = st.slider("Opacidade da superfície", 0.1, 1.0, 0.5, step=0.1)
                opt_auto_spin = st.checkbox("Giro automático (Spin)", value=False)
                opt_fundo_escuro = st.checkbox("Fundo Escuro", value=False)

    # RENDERIZAÇÃO CENTRAL
    with col_render:
        if modo_vis == "2D":
            img_2d = gerar_imagem_2d(
                smiles_input,
                mostrar_h=opt_mostrar_h,
                destacar_quirais=opt_destacar_q,
                mostrar_indices=opt_mostrar_idx,
                mostrar_metilas=opt_metilas,
                espessura_linha=opt_espessura,
                tamanho_fonte=opt_fonte,
                fundo_escuro=opt_fundo_escuro_2d,
                indices_destaque_extra=indices_grupo_selecionado
            )
            if img_2d:
                st.image(img_2d, caption="Estrutura 2D Interativa")
        else:
            html_3d = gerar_html_3d(
                smiles_input,
                estilo=opt_estilo_3d,
                mostrar_h=opt_mostrar_h_3d,
                fundo_escuro=opt_fundo_escuro,
                superficie=opt_superficie,
                coloracao_superficie=opt_coloracao,
                opacidade_superficie=opt_opacidade,
                auto_spin=opt_auto_spin,
                destacar_quirais=opt_destacar_q_3d,
                espessura_3d=opt_espessura_3d,
                esquema_cor=opt_esquema_cor_3d,
                raio_malha_quiral=opt_raio_malha
            )
            if html_3d:
                components.html(html_3d, height=430, width=540)
            else:
                st.warning("Não foi possível calcular a geometria 3D desta molécula.")

        st.divider()

        # Botões de Ação Compactos
        c_b1, c_b2, _ = st.columns([1.2, 1.5, 2])

        with c_b1:
            if st.button("🧬 Analisar Estereoquímica"):
                st.session_state.exibir_estereo = not st.session_state.exibir_estereo

        with c_b2:
            if st.button("🖼️ Verificar Isômeros Possíveis"):
                st.session_state.exibir_isomeros = not st.session_state.exibir_isomeros

        if st.session_state.exibir_estereo:
            dados_quirais = identificar_centros_quirais(smiles_input)
            st.markdown("---")
            st.subheader("Análise de Centros Quirais & Atividade Óptica")

            if dados_quirais and dados_quirais["tem_centros"]:
                status_optico = "✨ Opticamente Ativa" if dados_quirais["opticamente_ativo"] else "⚪ Inativa / Indefinida"
                st.markdown(f"**Status da Molécula:** `{status_optico}`")
                st.write(f"Encontrado(s) **{dados_quirais['total_centros']}** carbono(s) quiral(is):")

                for centro in dados_quirais["centros"]:
                    st.write(
                        f"• Carbono **{centro['rotulo_carbono']}** "
                        f"(Índice RDKit: `{centro['indice_rdkit']}`): Configuração **{centro['configuracao']}**"
                    )
            else:
                st.info("Esta molécula é aquiral (não possui centros quirais identificados).")

        if st.session_state.exibir_isomeros:
            st.markdown("---")
            st.subheader("Galeria de Estereoisômeros Possíveis ($2^n$)")

            with st.spinner("Enumerando estereoisômeros com o RDKit..."):
                isomeros_smiles = gerar_estereoisomeros(smiles_input)

            st.success(f"Total de combinações geradas ($2^n$): **{len(isomeros_smiles)}**")

            cols_por_linha = 4
            for i in range(0, len(isomeros_smiles), cols_por_linha):
                cols = st.columns(cols_por_linha)
                for j in range(cols_por_linha):
                    idx_iso = i + j
                    if idx_iso < len(isomeros_smiles):
                        smi_iso = isomeros_smiles[idx_iso]
                        img_iso = gerar_imagem_2d(smi_iso, espessura_linha=1, tamanho_fonte=10)
                        if img_iso:
                            cols[j].image(img_iso, caption=f"Isômero {idx_iso + 1}")
                            cols[j].code(smi_iso, language="text")
