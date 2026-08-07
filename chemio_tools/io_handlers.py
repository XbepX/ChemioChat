import io
from rdkit import Chem
from rdkit.Chem import AllChem
from .utils import validar_smiles


def carregar_mol_de_texto(conteudo_texto: str, formato: str = "MOL") -> str | None:
    """
    Converte arquivos de texto (.mol, .sdf, .pdb) em string SMILES.
    Base para carregar estruturas salvas em editores de moléculas.
    """
    if not conteudo_texto:
        return None

    mol = None
    formato_upper = formato.upper()

    try:
        if formato_upper in ["MOL", "SDF"]:
            mol = Chem.MolFromMolBlock(conteudo_texto)
        elif formato_upper == "PDB":
            mol = Chem.MolFromPDBBlock(conteudo_texto)

        if mol:
            return Chem.MolToSmiles(mol)
    except Exception:
        return None

    return None


def exportar_bloco_3d(smiles: str, formato: str = "SDF") -> str | None:
    """Exporta a estrutura 3D minimizada em formato SDF ou PDB para uso em editores externos."""
    if not validar_smiles(smiles):
        return None

    mol = Chem.MolFromSmiles(smiles.strip())
    mol_3d = Chem.AddHs(mol)

    if AllChem.EmbedMolecule(mol_3d, randomSeed=42) != 0:
        return None

    AllChem.MMFFOptimizeMolecule(mol_3d)

    formato_upper = formato.upper()
    if formato_upper == "SDF":
        return Chem.MolToMolBlock(mol_3d)
    elif formato_upper == "PDB":
        return Chem.MolToPDBBlock(mol_3d)

    return None
