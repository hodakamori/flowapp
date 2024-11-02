from typing import Any, Dict

import py3Dmol
import streamlit as st
from stmol import showmol


def visualize_block_details(barfi_result: Dict[Any, Any]) -> None:
    selected_block = st.selectbox("Select a block:", barfi_result.keys())
    input_keys = barfi_result[selected_block]["block"]._inputs.keys()
    output_keys = barfi_result[selected_block]["block"]._outputs.keys()
    block_keys = list(output_keys) + list(input_keys)
    if block_keys:
        selected_key = st.selectbox("Select:", block_keys)
        if selected_key:
            value = barfi_result[selected_block]["block"].get_interface(selected_key)
            if isinstance(value, py3Dmol.view):
                showmol(value)
            else:
                st.write(value)
