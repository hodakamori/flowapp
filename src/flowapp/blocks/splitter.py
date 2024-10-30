from typing import Any

import pandas as pd
from barfi import Block
from sklearn.model_selection import train_test_split


def train_test_splitter() -> Block:
    block = Block(name="Train/Test split")
    block.add_input(name="In(df)")
    block.add_output(name="Train(df)")
    block.add_output(name="Test(df)")

    def compute_func(self: Any) -> None:
        df = self.get_interface(name="In(df)")
        if df is not None and isinstance(df, pd.DataFrame):
            train_df, test_df = train_test_split(df, test_size=0.2, random_state=42)
            self.set_interface(name="Train(df)", value=train_df)
            self.set_interface(name="Test(df)", value=test_df)
        else:
            self.set_interface(name="Train(df)", value=None)
            self.set_interface(name="Test(df)", value=None)

    block.add_compute(compute_func)
    return block
