FROM python:slim

WORKDIR /app
RUN apt-get update && \
    apt-get install -y --no-install-recommends curl && \
    rm -rf /var/lib/apt/lists/*
COPY . .
RUN pip install streamlit kubernetes
EXPOSE 8501