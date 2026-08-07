import py3Dmol
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
from .utils import validar_smiles, obter_mapa_carbonos


def gerar_imagem_2d(
    smiles: str,
    mostrar_h: bool = False,
    destacar_quirais: bool = False,
    mostrar_indices: bool = False,
    mostrar_metilas: bool = False,
    espessura_linha: int = 2,
    tamanho_fonte: int = 12,
    fundo_escuro: bool = False,
    indices_destaque_extra: list[int] = None
):
    """Gera a imagem 2D customizada alinhada às opções e destaque de centros quirais."""
    if not validar_smiles(smiles):
        return None

    mol = Chem.MolFromSmiles(smiles.strip())

    if mostrar_h:
        mol = Chem.AddHs(mol)

    mapa_c = obter_mapa_carbonos(mol)

    if mostrar_indices:
        for atom in mol.GetAtoms():
            idx = atom.GetIdx()
            if idx in mapa_c:
                atom.SetProp("atomNote", mapa_c[idx])

    highlight_atoms = set()

    if destacar_quirais:
        # Força a verificação de estereocentros mesmo sem estereoquímica explícita
        mol_estereo = Chem.Mol(mol)
        Chem.AssignStereochemistry(mol_estereo, cleanIt=True, force=True, flagPossibleStereoCenters=True)
        centros = Chem.FindMolChiralCenters(mol_estereo, includeUnassigned=True)

        for idx, config in centros:
            highlight_atoms.add(idx)
            atomo = mol.GetAtomWithIdx(idx)
            note_existente = atomo.GetProp("atomNote") if atomo.HasProp("atomNote") else ""

            simbolo_config = f"({config})" if config != "?" else ""
            if note_existente:
                rotulo = f"{note_existente} {simbolo_config}*".strip()
            else:
                rotulo = f"{simbolo_config}*".strip()

            atomo.SetProp("atomNote", rotulo)

    if indices_destaque_extra:
        for idx in indices_destaque_extra:
            highlight_atoms.add(idx)

    opts = Draw.MolDrawOptions()
    opts.bondLineWidth = espessura_linha
    opts.minFontSize = tamanho_fonte
    opts.maxFontSize = tamanho_fonte + 6
    opts.explicitMethyl = mostrar_metilas

    if fundo_escuro:
        opts.clearBackground = False
        opts.setBackgroundColor((0.15, 0.17, 0.22, 1.0))

    return Draw.MolToImage(
        mol,
        size=(500, 420),
        options=opts,
        highlightAtoms=list(highlight_atoms) if highlight_atoms else None
    )


def gerar_html_3d(
    smiles: str,
    estilo: str = "Ball & Stick",
    mostrar_h: bool = True,
    fundo_escuro: bool = False,
    superficie: str = "Nenhuma",
    coloracao_superficie: str = "Mapa Eletrostático (MEP)",
    opacidade_superficie: float = 0.5,
    auto_spin: bool = False,
    destacar_quirais: bool = False,
    espessura_3d: int = 4,
    esquema_cor: str = "Elemento (Jmol)",
    raio_malha_quiral: float = 1.2
) -> str | None:
    """Gera a visualização 3D interativa sincronizando os centros quirais com o modelo final."""
    if not validar_smiles(smiles):
        return None

    mol = Chem.MolFromSmiles(smiles.strip())
    mol_3d = Chem.AddHs(mol) if mostrar_h else mol

    if AllChem.EmbedMolecule(mol_3d, randomSeed=42) != 0:
        return None

    AllChem.MMFFOptimizeMolecule(mol_3d)

    try:
        AllChem.ComputeGasteigerCharges(mol_3d)
    except Exception:
        pass

    mol_block = Chem.MolToMolBlock(mol_3d)

    view = py3Dmol.view(width=520, height=420)
    view.addModel(mol_block, "sdf")

    bg_color = "#1e1e1e" if fundo_escuro else "white"
    view.setBackgroundColor(bg_color)

    spec_cor = {}
    if esquema_cor == "Monocromático Escuro":
        spec_cor = {"color": "#111111"}
    elif esquema_cor == "Monocromático Claro":
        spec_cor = {"color": "#f0f0f0"}
    elif esquema_cor == "Carbono Preto (CPK)":
        spec_cor = {"colorscheme": "blackCarbon"}
    else:
        spec_cor = {"colorscheme": "Jmol"}

    if estilo == "Wireframe / Line":
        raio_wireframe = 0.02 * espessura_3d
        dict_stick = {"radius": raio_wireframe}
        dict_stick.update(spec_cor)
        view.setStyle({"stick": dict_stick})

    elif estilo == "Stick (Bastonetes)":
        dict_stick = {"radius": espessura_3d * 0.04}
        dict_stick.update(spec_cor)
        view.setStyle({"stick": dict_stick})

    elif estilo == "CPK / Esferas":
        dict_sphere = {"scale": 0.8}
        dict_sphere.update(spec_cor)
        view.setStyle({"sphere": dict_sphere})

    else:  # Ball & Stick
        dict_stick = {"radius": espessura_3d * 0.03}
        dict_stick.update(spec_cor)
        dict_sphere = {"scale": 0.25}
        dict_sphere.update(spec_cor)
        view.setStyle({"stick": dict_stick, "sphere": dict_sphere})

    # Destaque dos Centros Quirais no modelo 3D
    if destacar_quirais:
        Chem.AssignStereochemistry(mol_3d, cleanIt=True, force=True, flagPossibleStereoCenters=True)
        centros = Chem.FindMolChiralCenters(mol_3d, includeUnassigned=True)
        conf = mol_3d.GetConformer()

        for idx, _ in centros:
            pos = conf.GetAtomPosition(idx)
            view.addSphere({
                "center": {"x": pos.x, "y": pos.y, "z": pos.z},
                "radius": raio_malha_quiral,
                "color": "#00aaff",
                "wireframe": True,
                "opacity": 0.45
            })

    # Superfície Eletrônica
    if superficie != "Nenhuma":
        surf_type = py3Dmol.VDW
        if superficie == "Superfície Acessível ao Solvente (SAS)":
            surf_type = py3Dmol.SAS

        if coloracao_superficie == "Mapa Eletrostático (MEP)":
            surf_spec = {"opacity": opacidade_superficie, "colorscheme": "electrostatic"}
        else:
            surf_spec = {"opacity": opacidade_superficie, "color": "lightblue"}

        view.addSurface(surf_type, surf_spec)

    if auto_spin:
        view.spin(True)

    view.zoomTo()
    return view._make_html()