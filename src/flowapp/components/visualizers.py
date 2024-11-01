from typing import Any, Dict

import streamlit as st


def visualize_block_details(barfi_result: Dict[Any, Any]) -> None:
    selected_block = st.selectbox("Select a block:", barfi_result.keys())
    input_keys = barfi_result[selected_block]["block"]._inputs.keys()
    output_keys = barfi_result[selected_block]["block"]._outputs.keys()
    block_keys = list(input_keys) + list(output_keys)
    if block_keys:
        selected_key = st.selectbox("Select:", block_keys)
        if selected_key:
            value = barfi_result[selected_block]["block"].get_interface(selected_key)
            st.write(value)
