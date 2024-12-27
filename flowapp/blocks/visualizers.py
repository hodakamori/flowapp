from typing import Any

import py3Dmol
from barfi import Block
from rdkit import Chem

from flowapp.components.plots import column_selectable_scatter, create_parity_plot
from flowapp.utils.logger import log_exceptions


def scatter() -> Block:
    block = Block(name="Scatter")
    block.add_input(name="In(df)")
    block.add_output(name="Figure")
    block.add_output(name="Columns")

    @log_exceptions(block._name)
    def compute_func(self: Any) -> None:
        df = self.get_interface(name="In(df)")

        fig = column_selectable_scatter(df)

        self.set_interface(name="Figure", value=fig)
        self.set_interface(name="Columns", value=df.columns)

    block.add_compute(compute_func)
    return block


def parity_plot() -> Block:
    block = Block(name="Parity Plot")
    block.add_input(name="In(y_train_true)")
    block.add_input(name="In(y_train_pred)")
    block.add_input(name="In(y_test_true)")
    block.add_input(name="In(y_test_pred)")
    block.add_output(name="Figure")

    @log_exceptions(block._name)
    def compute_func(self: Any) -> None:
        y_train_true = self.get_interface(name="In(y_train_true)")
        y_train_pred = self.get_interface(name="In(y_train_pred)")
        y_test_true = self.get_interface(name="In(y_test_true)")
        y_test_pred = self.get_interface(name="In(y_test_pred)")
        fig = create_parity_plot(y_train_true, y_train_pred, y_test_true, y_test_pred)

        self.set_interface(name="Figure", value=fig)

    block.add_compute(compute_func)
    return block


def view_mol3d() -> Block:
    block = Block(name="3D Molecule Viewer")
    block.add_input(name="In(mol)")
    block.add_output(name="Figure")

    @log_exceptions(block._name)
    def compute_func(self: Any) -> None:
        mol = self.get_interface(name="In(mol)")
        sdf2 = Chem.MolToMolBlock(mol)
        view = py3Dmol.view(data=sdf2)
        view.setStyle({"stick": {}, "sphere": {"scale": 0.3}})
        self.set_interface(name="Figure", value=view)

    block.add_compute(compute_func)
    return block
