import tempfile

import psi4
from rdkit import Chem
from rdkit.Chem import rdchem, rdDetermineBonds


def psi4mol_to_rdmol(psi4mol: psi4.core.Molecule) -> rdchem.Mol:
    with tempfile.TemporaryDirectory() as tmpdir:
        xyz_file_path = f"{tmpdir}/temp.xyz"
        with open(xyz_file_path, "w") as f:
            f.write(psi4mol.save_string_xyz_file())
        mol = Chem.MolFromXYZFile(xyz_file_path)
        rdDetermineBonds.DetermineBonds(mol)
    return mol


def rdmol_to_psi4input(rdmol: rdchem.Mol) -> str:
    elements = [atom.GetSymbol() for atom in rdmol.GetAtoms()]
    positions = rdmol.GetConformer(0).GetPositions()
    psi4_geom_input = []
    for elem, pos in zip(elements, positions):
        psi4_geom_input.append(f"{elem} {pos[0]} {pos[1]} {pos[2]}")
    psi4_geom_input_join = "\n".join(psi4_geom_input)
    return psi4_geom_input_join
