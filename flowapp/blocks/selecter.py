from typing import Any

from barfi import Block


def xy_selecter() -> Block:
    block = Block(name="XY Selecter")

    block.add_input(name="In(df)")
    block.add_output(name="Out(X, df)")
    block.add_output(name="Out(y, df)")
    block.add_option(
        name="select y",
        type="input",
    )

    def compute_func(self: Any) -> None:
        df = self.get_interface(name="In(df)")
        y_col = self.get_option(name="select y")
        y = df[y_col]
        X = df.loc[:, df.columns != y_col]
        self.set_interface(name="Out(X, df)", value=X)
        self.set_interface(name="Out(y, df)", value=y)

    block.add_compute(compute_func)

    return block
