FROM continuumio/miniconda3:latest
RUN conda install -y -c conda-forge \
    psi4 \
    grpcio \
    grpcio-tools \
    numpy \
    && conda clean --all -f -y
WORKDIR /app
COPY . /app
EXPOSE 50051
RUN pip install .
