from typing import Any

from barfi import Block

from flowapp.utils.logger import log_exceptions
from flowapp.utils.regressor import RegressionTrainer, available_models


def train_regression() -> Block:
    block = Block(name="Regression")
    block.add_option(
        name="select model",
        type="select",
        items=["Linear", "Lasso", "Ridge", "SVR", "PLS", "RandomForest"],
    )

    block.add_input(name="In(X, df)")
    block.add_input(name="In(y, df)")
    block.add_output(name="Out(model)")

    @log_exceptions(block._name)
    def compute_func(self: Any) -> None:
        X = self.get_interface(name="In(X, df)")
        y = self.get_interface(name="In(y, df)")
        model_type = self.get_option(name="select model")

        trainer = RegressionTrainer(
            model_type=available_models[model_type],
        )
        trained_model = trainer.train(X, y)

        self.set_interface(name="Out(model)", value=trained_model)

    block.add_compute(compute_func)
    return block
