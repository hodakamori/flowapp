from typing import Any, Dict

import py3Dmol
import streamlit as st
from stmol import showmol


def visualize_block_details(barfi_result: Dict[Any, Any]) -> None:
    tabs = st.tabs(list(barfi_result.keys()))

    for tab, block_name in zip(tabs, barfi_result.keys()):
        with tab:
            input_keys = barfi_result[block_name]["block"]._inputs.keys()
            output_keys = barfi_result[block_name]["block"]._outputs.keys()
            block_keys = list(output_keys) + list(input_keys)

            if block_keys:
                selected_key = st.selectbox(
                    "Select:",
                    block_keys,
                    key=f"select_key_{block_name}",
                )

                if selected_key:
                    value = barfi_result[block_name]["block"].get_interface(
                        selected_key
                    )

                    if isinstance(value, py3Dmol.view):
                        showmol(value)
                    else:
                        st.write(value)

            else:
                st.info("No input/output keys available for this block.")
