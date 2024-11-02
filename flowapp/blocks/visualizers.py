from typing import Any

from barfi import Block

from flowapp.components.plots import column_selectable_scatter, create_parity_plot


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


def parity_plot() -> Block:
    block = Block(name="Parity Plot")
    block.add_input(name="In(y_train_true)")
    block.add_input(name="In(y_train_pred)")
    block.add_input(name="In(y_test_true)")
    block.add_input(name="In(y_test_pred)")
    block.add_output(name="Figure")

    def compute_func(self: Any) -> None:
        y_train_true = self.get_interface(name="In(y_train_true)")
        y_train_pred = self.get_interface(name="In(y_train_pred)")
        y_test_true = self.get_interface(name="In(y_test_true)")
        y_test_pred = self.get_interface(name="In(y_test_pred)")
        fig = create_parity_plot(y_train_true, y_train_pred, y_test_true, y_test_pred)

        self.set_interface(name="Figure", value=fig)

    block.add_compute(compute_func)
    return block
