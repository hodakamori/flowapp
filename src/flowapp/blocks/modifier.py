from typing import Any

from barfi import Block


def scaler() -> Block:
    block = Block(name="Scaling")
    block.add_option(
        name="select scaling type",
        type="select",
        items=["Normalization", "Standardization"],
    )
    block.add_input(name="In(df)")
    block.add_output(name="Out(df)")

    def compute_func(self: Any) -> None:
        df = self.get_interface(name="Input 1")
        scaling_type = self.get_option(name="select scaling type")
        if scaling_type == "Normalization":
            df = (df - df.min()) / (df.max() - df.min())
        elif scaling_type == "Standardization":
            df = (df - df.mean()) / df.std()
        self.set_interface(name="Out(df)", value=df)

    block.add_compute(compute_func)
    return block


def dropna() -> Block:
    block = Block(name="DropNaN")
    block.add_input(name="In(df)")
    block.add_output(name="Out(df)")

    def compute_func(self: Any) -> None:
        df = self.get_interface(name="In(df)")
        self.set_interface(name="Out(df)", value=df.dropna())

    block.add_compute(compute_func)
    return block
