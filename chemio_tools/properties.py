from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors
from .utils import validar_smiles

# Dicionário SMARTS de grupos funcionais de interesse
GRUPOS_SMARTS = {
    "Ácido Carboxílico": "[CX3](=O)[OX2H1]",
    "Aldeído": "[CX3H1](=O)",
    "Cetona": "[#6][CX3](=O)[#6]",
    "Fenol": "[OX2H1][c]",
    "Álcool": "[OX2H1][CX4]",
    "Amina": "[NX3;!$(NC=O)]",
    "Éster": "[CX3](=O)[OX2H0][#6]",
    "Éter": "[OD2]([#6])[#6]",
    "Amida": "[CX3](=O)[NX3]"
}


def obter_propriedades_basicas(smiles: str) -> dict | None:
    """Calcula massa molecular, LogP, TPSA e contagem atômica elementar."""
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


def identificar_grupos_funcionais(smiles: str) -> list[dict]:
    """Varre a estrutura em busca de grupos funcionais via SMARTS e retorna os índices associados."""
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
