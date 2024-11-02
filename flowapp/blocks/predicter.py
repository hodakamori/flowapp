from typing import Any

from barfi import Block


def predict() -> Block:
    block = Block(name="Predict")
    block.add_input(name="In(X, df)")
    block.add_input(name="In(model)")
    block.add_output(name="Out(y, df)")

    def compute_func(self: Any) -> None:
        X = self.get_interface(name="In(X, df)")
        model = self.get_interface(name="In(model)")
        y_pred = model.predict(X)
        self.set_interface(name="Out(y, df)", value=y_pred)

    block.add_compute(compute_func)
    return block
