import glob
from typing import Any

import pandas as pd
from barfi import Block


def csv_loader() -> Block:
    block = Block(name="CSV Loader")
    block.add_option(
        name="select data to load", type="select", items=glob.glob("data/*.csv")
    )
    block.add_output(name="out(df)")

    def compute_func(self: Any) -> None:
        path = self.get_option(name="select data to load")
        df = pd.read_csv(path)
        self.set_interface(name="out(df)", value=df)

    block.add_compute(compute_func)

    return block
