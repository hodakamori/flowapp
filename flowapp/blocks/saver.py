from typing import Any

from barfi.flow import Block


def save_as_csv() -> Block:
    block = Block(name="save")
    block.add_input(name="In(df)")
    block.add_option(
        name="explain", type="display", value="enter a file name without extension"
    )
    block.add_option(name="result", type="input", items="output")

    def save_func(self: Any) -> None:
        df = self.get_interface(name="In(df)")
        file_name = self.get_option(name="file name")
        df.to_csv(f"output/{file_name}.csv", index=False)

    block.add_compute(save_func)

    return block
