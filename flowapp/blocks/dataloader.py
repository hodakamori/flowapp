import glob
import os
from typing import Any

import pandas as pd
import streamlit as st
from barfi import Block
from dotenv import load_dotenv
from rdkit import Chem
from rdkit.Chem import AllChem

from flowapp.utils.logger import BlockExecutionLogger

load_dotenv()


def csv_loader() -> Block:
    block = Block(name="CSV")
    logger = BlockExecutionLogger.get_logger()

    demo_files = glob.glob("demo/*.csv")

    UPLOAD_ROOT_DIR = os.getenv("UPLOAD_ROOT_DIR", "uploads")
    UPLOAD_DIR = os.path.join(UPLOAD_ROOT_DIR, st.session_state.session_id)
    uploaded_files = glob.glob(f"{UPLOAD_DIR}/*.csv")
    block.add_option(
        name="select data to load", type="select", items=demo_files + uploaded_files
    )
    block.add_output(name="out(df)")

    def compute_func(self: Any) -> None:
        try:
            path = self.get_option(name="select data to load")
            df = pd.read_csv(path)
            self.set_interface(name="out(df)", value=df)
            logger.info(f"{self._name}: Suucessfully loaded {path}")
        except Exception as e:
            logger.error(f"{self._name}: Failed to load {path}")
            logger.error(e)

    block.add_compute(compute_func)

    return block


def pdb_loader() -> Block:
    block = Block(name="PDB")
    logger = BlockExecutionLogger.get_logger()

    demo_files = glob.glob("demo/*.pdb")

    UPLOAD_ROOT_DIR = os.getenv("UPLOAD_ROOT_DIR", "uploads")
    UPLOAD_DIR = os.path.join(UPLOAD_ROOT_DIR, st.session_state.session_id)
    uploaded_files = glob.glob(f"{UPLOAD_DIR}/*.pdb")
    block.add_option(
        name="select data to load", type="select", items=demo_files + uploaded_files
    )
    block.add_output(name="out(mol)")

    def compute_func(self: Any) -> None:
        try:
            path = self.get_option(name="select data to load")
            mol = Chem.MolFromPDBFile(path)
            mol = AllChem.AddHs(mol)
            AllChem.EmbedMolecule(mol)
            self.set_interface(name="out(mol)", value=mol)
            logger.info(f"{self._name}: Suucessfully loaded {path}")
        except Exception as e:
            logger.error(f"{self._name}: Failed to load {path}")
            logger.error(e)

    block.add_compute(compute_func)

    return block


def smiles_loader() -> Block:
    block = Block(name="SMILES")
    block.add_option(name="input SMILES", type="input")
    block.add_output(name="out(mol)")

    def compute_func(self: Any) -> None:
        smiles = self.get_option(name="input SMILES")
        mol = Chem.MolFromSmiles(smiles)
        mol = AllChem.AddHs(mol)
        AllChem.EmbedMolecule(mol)
        self.set_interface(name="out(mol)", value=mol)

    block.add_compute(compute_func)

    return block
