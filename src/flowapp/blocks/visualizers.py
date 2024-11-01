from typing import Any

from barfi import Block

from flowapp.components.plots import column_selectable_scatter


def scatter() -> Block:
    block = Block(name="Scatter")
    block.add_input(name="In(df)")
    block.add_output(name="Figure")
    block.add_output(name="Columns")

    def compute_func(self: Any) -> None:
        df = self.get_interface(name="In(df)")

        fig = column_selectable_scatter(df)

        self.set_interface(name="Figure", value=fig)
        self.set_interface(name="Columns", value=df.columns)

    block.add_compute(compute_func)
    return block
