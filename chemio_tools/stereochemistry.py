from rdkit import Chem
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers, StereoEnumerationOptions
from .utils import validar_smiles, obter_mapa_carbonos


def identificar_centros_quirais(smiles: str) -> dict | None:
    """Identifica centros quirais com configuração R/S e mapeamento C1, C2..."""
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


def gerar_estereoisomeros(smiles: str) -> list[str]:
    """Gera todas as combinações de estereoisômeros possíveis (2^n SMILES)."""
    if not validar_smiles(smiles):
        return []

    mol = Chem.MolFromSmiles(smiles.strip())
    opts = StereoEnumerationOptions(onlyUnassigned=False, unique=True)
    isomeros = list(EnumerateStereoisomers(mol, options=opts))

    return [Chem.MolToSmiles(iso, isomericSmiles=True) for iso in isomeros]
