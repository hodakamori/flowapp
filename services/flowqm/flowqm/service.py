import json
import os
import tempfile
from concurrent import futures
from typing import Any, Dict, Tuple

import grpc
import numpy as np
import psi4

from flowqm import qm_pb2, qm_pb2_grpc


class QMService(qm_pb2_grpc.QMServiceServicer):
    def __init__(self) -> None:
        """Initialize the QM service."""
        psi4.core.be_quiet()
        self.calculation_methods = {
            "opt": self._run_optimization,
        }

    def CalculateQM(
        self, request: qm_pb2.QMRequest, context: grpc.ServicerContext
    ) -> qm_pb2.QMResponse:
        try:
            mol_string = self._create_molecule_string(request.atoms)

            calc_type = request.params.calc_type
            theory = f"{request.params.theory_level}/{request.params.basis_set}"

            with tempfile.TemporaryDirectory() as tmpdir:
                psi4.set_output_file(f"{tmpdir}/psi4_output.dat", False)

                molecule = psi4.geometry(mol_string)

                calc_method = self.calculation_methods.get(calc_type)
                if calc_method is None:
                    raise ValueError(f"Unsupported calculation type: {calc_type}")

                result_atoms, energy, additional_data = calc_method(
                    molecule, theory, request.params.additional_params
                )

                if os.path.exists("timer.dat"):
                    os.remove("timer.dat")

                return qm_pb2.QMResponse(
                    success=True,
                    error_message="",
                    atoms=result_atoms,
                    energy=energy,
                    additional_data=json.dumps(additional_data),
                )

        except Exception as e:
            return qm_pb2.QMResponse(
                success=False,
                error_message=str(e),
                atoms=request.atoms,
                energy=0.0,
                additional_data="{}",
            )

    def _create_molecule_string(self, atoms: qm_pb2.Atoms) -> str:
        mol_lines = []

        mol_lines.append(f"{atoms.charge} {atoms.multiplicity}")

        coords = np.array(atoms.coordinates).reshape(-1, 3)
        for element, (x, y, z) in zip(atoms.elements, coords):
            mol_lines.append(f"{element} {x:.10f} {y:.10f} {z:.10f}")

        # if atoms.pbc and atoms.HasField("cell"):
        #     cell = np.array(atoms.cell.values).reshape(3, 3)
        #     mol_lines.append("units angstrom")
        #     mol_lines.append("no_reorient")
        #     mol_lines.append("no_com")

        return "\n".join(mol_lines)

    def _create_atoms_message(self, molecule: psi4.core.Molecule) -> qm_pb2.Atoms:
        elements = []
        coordinates = []

        for i in range(molecule.natom()):
            elements.append(molecule.symbol(i))
            coord = molecule.geometry().np[i]
            coordinates.extend(coord)

        return qm_pb2.Atoms(
            elements=elements,
            coordinates=coordinates,
            charge=molecule.molecular_charge(),
            multiplicity=molecule.multiplicity(),
        )

    def _run_optimization(
        self, molecule: psi4.core.Molecule, theory: str, options: Dict[str, str]
    ) -> Tuple[qm_pb2.Atoms, float, Dict[str, Any]]:
        energy, wfn, history = psi4.optimize(
            theory, molecule=molecule, return_wfn=True, return_history=True
        )
        additional_data = {
            "converged": True,
            "n_steps": len(history["energy"]),
            "energy_history": history["energy"],
        }

        return self._create_atoms_message(wfn.molecule()), energy, additional_data


def serve() -> None:
    server = grpc.server(futures.ThreadPoolExecutor(max_workers=4))
    qm_pb2_grpc.add_QMServiceServicer_to_server(QMService(), server)  # type: ignore[no-untyped-call]
    server.add_insecure_port("[::]:50051")
    server.start()
    print("QM Server started on port 50051")
    server.wait_for_termination()


if __name__ == "__main__":
    serve()
