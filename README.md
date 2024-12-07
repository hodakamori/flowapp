# flowapp
flowapp is a flow-based programming tool designed for materials informatics and computational science applications.

### Geometry optimization by Quantum mechanics
![](assets/qmopt_example.png)

### Regression of logS with SMILES
![](assets/regression_example.png)

## Features
- Simple and intuitive user interface
- Easy to extend: Add any custom processing by writing Python code
- Built for data pipeline management and monitoring

## Installation and Usage

You can install and run flowapp either using conda or Docker.

### Option 1: Using Conda

#### Prerequisites
flowapp requires conda and uses conda-lock for environment management.

Install conda-lock with the following command:

```bash
conda install -c conda-forge conda-lock
```

#### Installation

Install flowapp using these commands:

```bash
git clone https://github.com/hodakamori/flowapp.git
cd flowapp
conda-lock install -n flowapp conda-lock.yml
conda activate flowapp
pip install -e .
streamlit run flowapp/app.py
```

### Option 2: Using Docker

#### Prerequisites
- [Docker](https://www.docker.com/get-started)
- [Docker Compose](https://docs.docker.com/compose/install/)

#### Installation and Running

1. Clone the repository:
```bash
git clone https://github.com/hodakamori/flowapp.git
cd flowapp
```

2. Start the application using Docker Compose:
```bash
docker compose up --build
```

The application will be available at `http://localhost:8501`

To stop the application, use:
```bash
docker compose down
```