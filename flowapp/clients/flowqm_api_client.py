import json
from typing import Any, Dict, Optional

import grpc
from flowqm import qm_pb2, qm_pb2_grpc
from rdkit import Chem
from rdkit.Chem import AllChem, rdchem


class QMClient:
    def __init__(self, host: str = "flowqm", port: int = 50051):
        """Flowqm client."""
        self.channel = grpc.insecure_channel(f"{host}:{port}")
        self.stub = qm_pb2_grpc.QMServiceStub(self.channel)

    def compute_func(
        self,
        rdkit_mol: rdchem.Mol,
        calc_type: str = "opt",
        theory_level: str = "hf",
        basis_set: str = "sto-3g",
        additional_params: Optional[Dict[Any, Any]] = None,
    ):
        if additional_params is None:
            additional_params = {}

        elements, coords = self._convert_rdkit_mol_to_atoms(rdkit_mol)

        charge = Chem.rdmolops.GetFormalCharge(rdkit_mol)
        multiplicity = 1

        request = qm_pb2.QMRequest(
            atoms=qm_pb2.Atoms(
                elements=elements,
                coordinates=coords,
                charge=charge,
                multiplicity=multiplicity,
            ),
            params=qm_pb2.CalculationParams(
                calc_type=calc_type,
                theory_level=theory_level,
                basis_set=basis_set,
                additional_params=additional_params,
            ),
        )

        response = self.stub.CalculateQM(request)

        if not response.success:
            raise RuntimeError(f"QM calculation failed: {response.error_message}")

        optimized_mol = self._update_rdkit_mol_coords(rdkit_mol, response.atoms)

        energy = response.energy
        additional_data = response.additional_data
        additional_data = json.loads(additional_data)
        return energy, optimized_mol, additional_data

    def _convert_rdkit_mol_to_atoms(self, rdkit_mol: rdchem.Mol):
        elements = []
        coords = []
        conf = rdkit_mol.GetConformer()
        for atom in rdkit_mol.GetAtoms():
            elements.append(atom.GetSymbol())
            pos = conf.GetAtomPosition(atom.GetIdx())
            coords.extend([pos.x, pos.y, pos.z])
        return elements, coords

    def _update_rdkit_mol_coords(self, rdkit_mol: rdchem.Mol, atoms_msg: qm_pb2.Atoms):
        new_mol = Chem.RWMol(rdkit_mol)
        conf = new_mol.GetConformer()

        new_coords = atoms_msg.coordinates
        for i, atom in enumerate(new_mol.GetAtoms()):
            x = new_coords[3 * i]
            y = new_coords[3 * i + 1]
            z = new_coords[3 * i + 2]
            conf.SetAtomPosition(i, (x, y, z))

        new_mol.UpdatePropertyCache()
        return Chem.Mol(new_mol)


def main():
    client = QMClient(host="localhost", port=50051)
    mol = Chem.MolFromSmiles("O")
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol)

    energy, optimized_mol = client.compute_func(
        rdkit_mol=mol,
        calc_type="opt",
        theory_level="hf",
        basis_set="sto-3g",
        additional_params={"key1": "value1"},
    )

    print("Energy =", energy)
    print("Optimized Mol:", Chem.MolToMolBlock(optimized_mol))


if __name__ == "__main__":
    main()
