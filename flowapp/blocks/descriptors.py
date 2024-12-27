from typing import Any

import numpy as np
import pandas as pd
from barfi import Block
from rdkit import Chem
from sklearn.compose import ColumnTransformer
from sklearn.preprocessing import OneHotEncoder

from flowapp.utils.logger import log_exceptions


def onehot_encoding() -> Block:
    block = Block(name="One hot encoding")
    block.add_input(name="In(df)")
    block.add_output(name="Out(df)")
    block.add_option(
        name="select y",
        type="input",
    )

    @log_exceptions(block._name)
    def compute_func(self: Any) -> None:
        df = self.get_interface(name="In(df)")
        col = self.get_option(name="select y")
        categorical_columns = [col]
        transformer = ColumnTransformer(
            transformers=[
                (
                    "onehot",
                    OneHotEncoder(handle_unknown="ignore"),
                    categorical_columns,
                )
            ],
            remainder="passthrough",
        )

        encoded_array = transformer.fit_transform(df)

        encoded_feature_names = transformer.named_transformers_[
            "onehot"
        ].get_feature_names_out(categorical_columns)

        non_encoded_columns = [
            col for col in df.columns if col not in categorical_columns
        ]

        all_feature_names = np.concatenate([encoded_feature_names, non_encoded_columns])

        encoded_df = pd.DataFrame(
            encoded_array, columns=all_feature_names, index=df.index
        )
        self.set_interface(name="Out(df)", value=encoded_df)

    block.add_compute(compute_func)
    return block


def smi2fp() -> Block:
    block = Block(name="Smiles to Fingerprint")
    block.add_input(name="In(df)")
    block.add_output(name="Out(df)")
    block.add_option(
        name="select y",
        type="input",
    )

    @log_exceptions(block._name)
    def compute_func(self: Any) -> None:
        df = self.get_interface(name="In(df)")
        col = self.get_option(name="select y")
        rdkitfp_feat = []
        for smi in df[col].values:
            mol = Chem.MolFromSmiles(smi)
            fp_rdkit = Chem.RDKFingerprint(mol)
            rdkitfp_feat.append(fp_rdkit.ToList())
        rdkitfp_feat_df = pd.DataFrame(
            rdkitfp_feat,
        )
        rdkitfp_feat_df.columns = [
            f"rdkit_{i}" for i in range(rdkitfp_feat_df.shape[1])
        ]
        encoded_df = pd.concat([df, rdkitfp_feat_df], axis=1).drop(col, axis=1)
        self.set_interface(name="Out(df)", value=encoded_df)

    block.add_compute(compute_func)
    return block
