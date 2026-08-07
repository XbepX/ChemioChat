import math
import py3Dmol
from rdkit import Chem
from rdkit.Chem import AllChem, Draw, Descriptors, rdMolDescriptors
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers, StereoEnumerationOptions

# Mapeamento de padrões SMARTS para grupos funcionais comuns
GRUPOS_SMARTS = {
    "Ácido Carboxílico": "[CX3](=O)[OX2H1]",
    "Aldeído": "[CX3H1](=O)",
    "Cetona": "[#6][CX3](=O)[#6]",
    "Álcool": "[OX2H1][CX4]",
    "Amina": "[NX3;!$(NC=O)]",
    "Éster": "[CX3](=O)[OX2H0][#6]",
    "Éter": "[OD2]([#6])[#6]",
    "Amida": "[CX3](=O)[NX3]"
}


def validar_smiles(smiles: str) -> bool:
    """Valida se a string SMILES é reconhecida pelo RDKit."""
    if not smiles or not smiles.strip():
        return False
    mol = Chem.MolFromSmiles(smiles.strip())
    return mol is not None


def obter_mapa_carbonos(mol: Chem.Mol) -> dict[int, str]:
    """Mapeia os índices dos átomos do RDKit para a numeração sequencial de Carbonos (C1, C2...)."""
    mapa = {}
    contador_c = 1
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == "C":
            mapa[atom.GetIdx()] = f"C{contador_c}"
            contador_c += 1
    return mapa


def obter_propriedades_basicas(smiles: str) -> dict | None:
    """Calcula massa molecular, LogP, TPSA e contagem de elementos."""
    if not validar_smiles(smiles):
        return None

    mol = Chem.MolFromSmiles(smiles.strip())
    mol_com_h = Chem.AddHs(mol)

    contagem_elementos = {}
    for atom in mol_com_h.GetAtoms():
        simbolo = atom.GetSymbol()
        contagem_elementos[simbolo] = contagem_elementos.get(simbolo, 0) + 1

    return {
        "formula": rdMolDescriptors.CalcMolFormula(mol),
        "massa_molecular": round(Descriptors.MolWt(mol), 2),
        "logp": round(Descriptors.MolLogP(mol), 2),
        "tpsa": round(Descriptors.TPSA(mol), 2),
        "contagem_elementos": contagem_elementos
    }


def identificar_centros_quirais(smiles: str) -> dict | None:
    """Identifica centros quirais alinhados à numeração C1, C2 e avalia atividade óptica."""
    if not validar_smiles(smiles):
        return None

    mol = Chem.MolFromSmiles(smiles.strip())
    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    centros = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
    mapa_c = obter_mapa_carbonos(mol)

    resultado_centros = []
    for idx, config in centros:
        atomo = mol.GetAtomWithIdx(idx)
        rotulo_c = mapa_c.get(idx, f"Átomo {idx}")
        resultado_centros.append({
            "indice_rdkit": idx,
            "rotulo_carbono": rotulo_c,
            "elemento": atomo.GetSymbol(),
            "configuracao": config if config != "?" else "Não especificada (?)"
        })

    is_opticamente_ativo = len(resultado_centros) > 0 and any(
        c["configuracao"] in ["R", "S"] for c in resultado_centros
    )

    return {
        "centros": resultado_centros,
        "total_centros": len(resultado_centros),
        "opticamente_ativo": is_opticamente_ativo,
        "tem_centros": len(resultado_centros) > 0
    }


def identificar_grupos_funcionais(smiles: str) -> list[dict]:
    """Varre a estrutura em busca de grupos funcionais via SMARTS."""
    if not validar_smiles(smiles):
        return []

    mol = Chem.MolFromSmiles(smiles.strip())
    grupos_encontrados = []

    for nome_grupo, smarts in GRUPOS_SMARTS.items():
        patt = Chem.MolFromSmarts(smarts)
        if patt and mol.HasSubstructMatch(patt):
            matches = mol.GetSubstructMatches(patt)
            indices_unicos = list(set([idx for tuple_match in matches for idx in tuple_match]))
            grupos_encontrados.append({
                "grupo": nome_grupo,
                "indices": indices_unicos
            })

    return grupos_encontrados


def gerar_estereoisomeros(smiles: str) -> list[str]:
    """Gera todas as combinações de estereoisômeros possíveis (2^n SMILES)."""
    if not validar_smiles(smiles):
        return []

    mol = Chem.MolFromSmiles(smiles.strip())
    opts = StereoEnumerationOptions(onlyUnassigned=False, unique=True)
    isomeros = list(EnumerateStereoisomers(mol, options=opts))

    return [Chem.MolToSmiles(iso, isomericSmiles=True) for iso in isomeros]


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
    """Gera a imagem 2D customizada adicionando asterisco (*) nos carbonos quirais."""
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

    # Adiciona o destaque e a marcação tradicional com asterisco (*)
    if destacar_quirais:
        Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
        centros = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
        for idx, config in centros:
            highlight_atoms.add(idx)
            atomo = mol.GetAtomWithIdx(idx)
            note_existente = atomo.GetProp("atomNote") if atomo.HasProp("atomNote") else ""

            if note_existente:
                rotulo = f"{note_existente} ({config})*" if config != "?" else f"{note_existente}*"
            else:
                rotulo = f"({config})*" if config != "?" else "*"

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


# função para gerar a visualização 3D interativa
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
    esquema_cor: str = "Elemento (Jmol)"
) -> str | None:
    """Gera visualização 3D interativa com malha quiral de diâmetro customizável."""
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

    else:
        dict_stick = {"radius": espessura_3d * 0.03}
        dict_stick.update(spec_cor)
        dict_sphere = {"scale": 0.25}
        dict_sphere.update(spec_cor)
        view.setStyle({"stick": dict_stick, "sphere": dict_sphere})

    # Destaque dos Centros Quirais com Malha Azul de Diâmetro Personalizado
    if destacar_quirais:
        Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
        centros = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
        conf = mol_3d.GetConformer()

        for idx, _ in centros:
            pos = conf.GetAtomPosition(idx)
            view.addSphere({
                "center": {"x": pos.x, "y": pos.y, "z": pos.z},
                "radius": 1.2,  # <--- Raio em Ångströms (Diâmetro = 2 * radius)
                "color": "#00aaff",
                "wireframe": True,
                "opacity": 0.45
            })

    # Superfície Eletrônica Geral
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