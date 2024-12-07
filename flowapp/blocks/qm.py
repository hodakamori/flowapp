import os
import tempfile
from typing import Any

import plotly.express as px
import psi4
from barfi import Block

from flowapp.utils.converter import psi4mol_to_rdmol, rdmol_to_psi4input
from flowapp.utils.logger import log_exceptions


def opt() -> Block:
    block = Block(name="Opt")
    block.add_input(name="In(mol)")
    block.add_output(name="Out(mol)")
    block.add_output(name="Out(energy)")
    block.add_output(name="Out(wfn)")
    block.add_option(
        name="select calc level",
        type="select",
        items=["b3lyp/6-31g"],
        value="b3lyp/6-31g",
    )

    @log_exceptions(block._name)
    def compute_func(self: Any) -> None:
        mol = self.get_interface(name="In(mol)")
        calc_level = self.get_option(name="select calc level")
        psi4_geom_input = rdmol_to_psi4input(mol)
        with tempfile.TemporaryDirectory() as tmpdir:
            psi4.set_output_file(f"{tmpdir}/psi4_output.dat", False)
            _ = psi4.geometry(psi4_geom_input)
            _, opt_wfn, history = psi4.optimize(
                calc_level, return_wfn=True, return_history=True
            )
            os.remove("timer.dat")
        rdkit_mol = psi4mol_to_rdmol(opt_wfn.molecule())
        energy_values = history["energy"]
        fig = px.line(
            x=range(len(energy_values)),
            y=energy_values,
            labels={"x": "Step", "y": "Energy"},
        )
        self.set_interface(name="Out(energy)", value=fig)
        self.set_interface(name="Out(wfn)", value=opt_wfn)
        self.set_interface(name="Out(mol)", value=rdkit_mol)

    block.add_compute(compute_func)
    return block
