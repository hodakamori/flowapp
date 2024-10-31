from typing import Any

from barfi import Block
from sklearn.linear_model import LinearRegression


def train_regression() -> Block:
    block = Block(name="Train Regression")
    block.add_option(
        name="select model",
        type="select",
        items=["Linear"],
    )

    block.add_input(name="In(X, df)")
    block.add_input(name="In(y, df)")
    block.add_output(name="Out(model)")

    def compute_func(self: Any) -> None:
        X = self.get_interface(name="In(X, df)")
        y = self.get_interface(name="In(y, df)")
        model = self.get_option(name="select model")
        if model == "Linear":
            model = LinearRegression()
            model.fit(X, y)
        self.set_interface(name="Out(model)", value=model)

    block.add_compute(compute_func)
    return block
