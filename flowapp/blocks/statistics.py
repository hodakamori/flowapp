import json
from typing import Any

from barfi.flow import Block
from ydata_profiling import ProfileReport

from flowapp.utils.logger import log_exceptions


def profiling() -> Block:
    block = Block(name="Profiling")
    block.add_input(name="In(df)")
    block.add_output(name="Out(profile)")

    @log_exceptions(block.name)
    def compute_func(self: Any) -> None:
        df = self.get_interface(name="In(df)")
        profile = ProfileReport(df, title="Profiling Report")
        data = json.loads(profile.to_json())
        table_items = [
            "n",
            "n_var",
            "n_cells_missing",
            "n_vars_with_missing",
            "n_vars_all_missing",
            "types",
            "n_duplicates",
        ]
        table_view = {}

        for key in data["table"].keys():
            if key in table_items:
                table_view[key] = data["table"][key]

        self.set_interface(name="Out(profile)", value=table_view)

    block.add_compute(compute_func)
    return block
