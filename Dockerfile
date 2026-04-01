# Base image
FROM python:3.10-slim

# System dependencies for RDKit
RUN apt-get update && apt-get install -y \
    build-essential \
    libglib2.0-0 \
    libxrender1 \
    libxext6 \
    libsm6 \
    && rm -rf /var/lib/apt/lists/*

# Set working directory
WORKDIR /app

# Copy files
COPY . /app

# Install Python dependencies
RUN pip install --no-cache-dir --upgrade pip
RUN pip install --no-cache-dir -r requirements.txt

# Expose port
EXPOSE 8501

# Run app
CMD ["streamlit", "run", "chem.py", "--server.port=8501", "--server.address=0.0.0.0"]
