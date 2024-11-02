import os
import uuid
from pathlib import Path
from typing import List

import streamlit as st
from barfi import barfi_schemas, st_barfi
from dotenv import load_dotenv

from flowapp.blocks.dataloader import csv_loader
from flowapp.blocks.descriptors import onehot_encoding, smi2fp
from flowapp.blocks.evaluator import regression_score
from flowapp.blocks.modifier import dropcolumn, dropna, scaler
from flowapp.blocks.predicter import predict
from flowapp.blocks.regression import train_regression
from flowapp.blocks.saver import save_as_csv
from flowapp.blocks.selecter import xy_selecter
from flowapp.blocks.splitter import train_test_splitter
from flowapp.blocks.visualizers import parity_plot, scatter
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

col1, col2 = st.columns([7, 3])

with col1:
    st.header("Flow")
    saved_schemas: List[str] = barfi_schemas()
    load_schema = st.selectbox("Select a saved schema:", saved_schemas)
    barfi_result = st_barfi(
        base_blocks=[
            csv_loader(),
            scaler(),
            train_test_splitter(),
            dropna(),
            dropcolumn(),
            onehot_encoding(),
            smi2fp(),
            xy_selecter(),
            train_regression(),
            predict(),
            regression_score(),
            save_as_csv(),
            scatter(),
            parity_plot(),
        ],
        compute_engine=True,
        load_schema=load_schema,
    )

with col2:
    st.header("Block Details")
    if barfi_result:
        visualize_block_details(barfi_result)
