from rdkit import Chem


def validar_smiles(smiles: str) -> bool:
    """Valida se a string SMILES é reconhecida pelo RDKit."""
    if not smiles or not smiles.strip():
        return False
    mol = Chem.MolFromSmiles(smiles.strip())
    return mol is not None


def obter_mapa_carbonos(mol: Chem.Mol) -> dict[int, str]:
    """
    Mapeia os índices dos átomos do RDKit para a numeração sequencial de Carbonos (C1, C2...).
    Utiliza a topologia do grafo de leitura para garantir estabilidade e previsibilidade.
    """
    mapa = {}
    contador_c = 1
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == "C":
            mapa[atom.GetIdx()] = f"C{contador_c}"
            contador_c += 1
    return mapa