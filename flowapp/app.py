import os
import uuid
from pathlib import Path

import streamlit as st
from barfi.flow import ComputeEngine, SchemaManager
from barfi.flow.streamlit import st_flow
from dotenv import load_dotenv

from flowapp.blocks.dataloader import csv_loader, pdb_loader, smiles_loader
from flowapp.blocks.descriptors import onehot_encoding, smi2fp
from flowapp.blocks.evaluator import regression_score
from flowapp.blocks.modifier import dropcolumn, dropna, scaler
from flowapp.blocks.predicter import predict
from flowapp.blocks.qm import opt
from flowapp.blocks.regression import train_regression
from flowapp.blocks.saver import save_as_csv
from flowapp.blocks.selecter import xy_selecter
from flowapp.blocks.splitter import train_test_splitter
from flowapp.blocks.statistics import profiling
from flowapp.blocks.visualizers import parity_plot, scatter, view_mol3d
from flowapp.components.block_visualizer import visualize_block_details
from flowapp.components.file_handler import file_uploader

st.set_page_config(layout="wide")
load_dotenv()

if "session_id" not in st.session_state:
    st.session_state.session_id = str(uuid.uuid4())

UPLOAD_ROOT_DIR = os.getenv("UPLOAD_ROOT_DIR", "uploads")
upload_root_path = Path(UPLOAD_ROOT_DIR)
upload_root_path.mkdir(exist_ok=True)

with st.sidebar:
    file_uploader(upload_root_path, st.session_state.session_id)

col1, col2 = st.columns([6, 4])
blocks = {
    "Data": [
        csv_loader(),
        pdb_loader(),
        smiles_loader(),
    ],
    "Statistics": [
        profiling(),
    ],
    "Modifier": [scaler(), dropna(), dropcolumn(), onehot_encoding(), smi2fp()],
    "Selecter": [
        train_test_splitter(),
        xy_selecter(),
    ],
    "QM": [
        opt(),
    ],
    "ML": [
        train_regression(),
        predict(),
        regression_score(),
    ],
    "Visualize": [
        scatter(),
        parity_plot(),
        view_mol3d(),
    ],
    "Othres": [
        save_as_csv(),
    ],
}

with col1:
    st.header("Flow")
    schema_manager = SchemaManager(filepath="./schemas")
    load_schema = st.selectbox(
        "Select a saved schema:", [""] + schema_manager.schema_names
    )
    if load_schema and load_schema != "None":
        flow_schema = schema_manager.load_schema(load_schema)
        barfi_result = st_flow(blocks=blocks, editor_schema=flow_schema)
    else:
        barfi_result = st_flow(blocks=blocks)
    compute_engine = ComputeEngine(blocks)

    @st.dialog("Save schema")
    def save_schema(barfi_result):
        placeholder = load_schema if load_schema else ""
        saved_name = st.text_input("Plese input name:", value=placeholder)
        if st.button("Save"):
            if saved_name in schema_manager.schema_names:
                schema_manager.update_schema(saved_name, barfi_result.editor_schema)
            else:
                schema_manager.save_schema(saved_name, barfi_result.editor_schema)

    if barfi_result.command == "save":
        save_schema(barfi_result)
    flow_schema = barfi_result.editor_schema
    compute_engine.execute(flow_schema)

with col2:
    st.header("Block Details")
    if barfi_result:
        visualize_block_details(barfi_result)
