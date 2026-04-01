FROM continuumio/miniconda3

WORKDIR /app

# Copy files
COPY . /app

# Create conda environment with RDKit
RUN conda create -n chem-env python=3.10 -y && \
    conda install -n chem-env -c conda-forge rdkit streamlit pandas numpy pillow -y

# Activate environment
SHELL ["conda", "run", "-n", "chem-env", "/bin/bash", "-c"]

EXPOSE 8501

CMD ["conda", "run", "-n", "chem-env", "streamlit", "run", "chem.py", "--server.port=8501", "--server.address=0.0.0.0"]