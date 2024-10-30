from typing import List

import streamlit as st
from barfi import barfi_schemas, st_barfi

from flowapp.blocks.dataloader import csv_loader
from flowapp.blocks.evaluator import regression_score
from flowapp.blocks.ml import train_regression
from flowapp.blocks.modifier import dropna, scaler
from flowapp.blocks.predicter import predict
from flowapp.blocks.saver import save_as_csv
from flowapp.blocks.selecter import xy_selecter
from flowapp.blocks.splitter import train_test_splitter

st.set_page_config(layout="wide")
saved_schemas: List[str] = barfi_schemas()
load_schema = st.selectbox("Select a saved schema:", saved_schemas)

barfi_result = st_barfi(
    base_blocks=[
        csv_loader(),
        scaler(),
        train_test_splitter(),
        dropna(),
        xy_selecter(),
        train_regression(),
        predict(),
        regression_score(),
        save_as_csv(),
    ],
    compute_engine=True,
    load_schema=load_schema,
)

if barfi_result:
    st.write(barfi_result)
