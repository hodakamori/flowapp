from typing import Any

import plotly.express as px

# import psi4
from barfi import Block

from flowapp.clients.flowqm_api_client import QMClient
from flowapp.utils.logger import log_exceptions


def opt() -> Block:
    block = Block(name="Opt")
    block.add_input(name="In(mol)")
    block.add_output(name="Out(mol)")
    block.add_output(name="Out(energy)")
    block.add_option(
        name="Basis set",
        type="select",
        items=["6-31g"],
        value="6-31g",
    )
    block.add_option(
        name="Theory",
        type="select",
        items=["b3lyp"],
        value="b3lyp",
    )

    @log_exceptions(block._name)
    def compute_func(self: Any) -> None:
        mol = self.get_interface(name="In(mol)")
        basis_set = self.get_option(name="Basis set")
        theory = self.get_option(name="Theory")
        client = QMClient(host="flowqm", port=50051)
        energy, optimized_mol, additional_data = client.compute_func(
            rdkit_mol=mol,
            calc_type="opt",
            theory_level=theory,
            basis_set=basis_set,
            additional_params={"key1": "value1"},
        )
        energy_history = additional_data["energy_history"]
        fig = px.line(
            x=range(len(energy_history)),
            y=energy_history,
            labels={"x": "Step", "y": "Energy"},
        )
        self.set_interface(name="Out(energy)", value=fig)
        self.set_interface(name="Out(mol)", value=optimized_mol)

    block.add_compute(compute_func)
    return block
