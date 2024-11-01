from typing import List

import streamlit as st
from barfi import barfi_schemas, st_barfi

from flowapp.blocks.dataloader import csv_loader
from flowapp.blocks.descriptors import onehot_encoding
from flowapp.blocks.evaluator import regression_score
from flowapp.blocks.modifier import dropcolumn, dropna, scaler
from flowapp.blocks.predicter import predict
from flowapp.blocks.regression import train_regression
from flowapp.blocks.saver import save_as_csv
from flowapp.blocks.selecter import xy_selecter
from flowapp.blocks.splitter import train_test_splitter
from flowapp.blocks.visualizers import parity_plot, scatter
from flowapp.components.visualizers import visualize_block_details

st.set_page_config(layout="wide")

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
