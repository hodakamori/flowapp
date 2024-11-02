# flowapp

flowapp is a tool designed to simplify data flow management and processing.
This project aims to make it easy to build, manage, and monitor data pipelines.

## Features

- **Intuitive UI**: User-friendly interface for easy design and management of data flows.
- **Extensibility**: Plugin architecture allows for easy addition of custom features.

## Installation

You can install flowapp using the following commands:

```sh
git clone https://github.com/hodakamori/flowapp.git
cd flowapp
conda-lock install -n flowapp conda-lock.yml
conda activate flowapp
pip install -e .
streamlit run flowapp/app.py
```