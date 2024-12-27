import os
from pathlib import Path

import streamlit as st


def file_uploader(upload_root_path: Path, session_id: str) -> None:
    UPLOAD_DIR = os.path.join(upload_root_path, session_id)
    upload_path = Path(UPLOAD_DIR)
    upload_path.mkdir(exist_ok=True)
    st.header("Upload data")
    uploaded_file = st.file_uploader("Upload data", type=["csv", "pdb"])
    if uploaded_file is not None:
        file_path = os.path.join(UPLOAD_DIR, uploaded_file.name)
        with open(file_path, "wb") as f:
            f.write(uploaded_file.getbuffer())
