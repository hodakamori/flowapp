import glob
import os
from typing import Any

import pandas as pd
import streamlit as st
from barfi import Block
from dotenv import load_dotenv

load_dotenv()


def csv_loader() -> Block:
    block = Block(name="CSV Loader")
    demo_files = glob.glob("demo/*.csv")

    UPLOAD_ROOT_DIR = os.getenv("UPLOAD_ROOT_DIR", "uploads")
    UPLOAD_DIR = os.path.join(UPLOAD_ROOT_DIR, st.session_state.session_id)
    uploaded_files = glob.glob(f"{UPLOAD_DIR}/*.csv")
    block.add_option(
        name="select data to load", type="select", items=demo_files + uploaded_files
    )
    block.add_output(name="out(df)")

    def compute_func(self: Any) -> None:
        path = self.get_option(name="select data to load")
        df = pd.read_csv(path)
        self.set_interface(name="out(df)", value=df)

    block.add_compute(compute_func)

    return block
