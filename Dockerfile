FROM continuumio/miniconda3:latest

WORKDIR /app
COPY . .
RUN conda install -c conda-forge -y mamba && \
    mamba install -c conda-forge -y conda-lock && \
    conda-lock install -n flowapp conda-lock.yml && \
    echo "conda activate flowapp" >> ~/.bashrc && \
    /bin/bash -c "source activate flowapp && pip install -e ."

EXPOSE 8501
