import py3Dmol
import streamlit as st
from barfi.flow.streamlit.types import StreamlitFlowResponse
from stmol import showmol


def visualize_block_details(barfi_result: StreamlitFlowResponse) -> None:
    if len(barfi_result.editor_schema.nodes) > 0:
        labels = [node.label for node in barfi_result.editor_schema.nodes]
        tabs = st.tabs(labels)
        for tab, block, label in zip(tabs, barfi_result.editor_schema.nodes, labels):
            with tab:
                block = barfi_result.editor_schema.block(node_label=label)
                input_keys = list(block._inputs.keys())
                output_keys = list(block._outputs.keys())
                block_keys = output_keys + input_keys

                if block_keys:
                    selected_key = st.selectbox(
                        "Select:",
                        block_keys,
                        key=f"select_key_{label}",
                    )

                    if selected_key:
                        with st.container():
                            value = block.get_interface(selected_key)
                            st.write("\n" * 2)
                            if isinstance(value, py3Dmol.view):
                                showmol(value)
                            else:
                                st.write(value)

                else:
                    st.info("No input/output keys available for this block.")
