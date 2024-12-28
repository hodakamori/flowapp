from typing import Any, Dict, Optional

import pandas as pd
from sklearn.cross_decomposition import PLSRegression
from sklearn.ensemble import RandomForestRegressor
from sklearn.linear_model import Lasso, LinearRegression, Ridge
from sklearn.svm import SVR

available_models = {
    "Linear": "linear",
    "Lasso": "lasso",
    "Ridge": "ridge",
    "SVR": "svr",
    "PLS": "pls",
    "RandomForest": "random_forest",
}


class RegressionTrainer:
    def __init__(
        self, model_type: str = "linear", hyperparams: Optional[Dict[str, Any]] = None
    ):
        """Regression model trainer."""
        self.model_type = model_type.lower()
        self.hyperparams = hyperparams if hyperparams is not None else {}

    def train(self, X: pd.DataFrame, y: pd.DataFrame):
        if self.model_type == "linear":
            model = LinearRegression(**self.hyperparams)
        elif self.model_type == "lasso":
            model = Lasso(**self.hyperparams)
        elif self.model_type == "ridge":
            model = Ridge(**self.hyperparams)
        elif self.model_type == "pls":
            model = PLSRegression(**self.hyperparams)
        elif self.model_type == "svr":
            model = SVR(**self.hyperparams)
        elif self.model_type == "random_forest":
            model = RandomForestRegressor(**self.hyperparams)
        else:
            raise ValueError(f"Unknown model type: {self.model_type}")

        model.fit(X, y)
        return model
