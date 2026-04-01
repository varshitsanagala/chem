FROM python:3.10

# Install system dependencies for RDKit
RUN apt-get update && apt-get install -y \
    build-essential \
    libglib2.0-0 \
    libxrender1 \
    libxext6 \
    libsm6 \
    libboost-all-dev \
    libeigen3-dev \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

COPY . /app

RUN pip install --upgrade pip
RUN pip install --no-cache-dir -r requirements.txt

EXPOSE 8501

CMD ["streamlit", "run", "chem.py", "--server.port=8501", "--server.address=0.0.0.0"]