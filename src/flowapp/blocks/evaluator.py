from typing import Any

from barfi import Block
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score


def regression_score() -> Block:
    block = Block(name="Evaluator")

    block.add_input(name="In(y_true, df)")
    block.add_input(name="Out(y_pred, df)")
    block.add_output(name="Score")

    block.add_option(
        name="select metric",
        type="select",
        items=["MAE", "R2", "RMSE"],
    )

    def compute_func(self: Any) -> None:
        metric = self.get_option(name="select metric")
        y_true = self.get_interface(name="In(y_true, df)")
        y_pred = self.get_interface(name="Out(y_pred, df)")

        if metric == "MAE":
            score = mean_absolute_error(y_true, y_pred)
        elif metric == "R2":
            score = r2_score(y_true, y_pred)
        elif metric == "RMSE":
            score = mean_squared_error(y_true, y_pred, squared=False)
        self.set_interface(name="Score", value=score)

    block.add_compute(compute_func)
    return block
