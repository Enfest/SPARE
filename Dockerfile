FROM python:3.8.10

RUN apt-get update && \
    apt-get install -y --no-install-recommends vim curl bzip2 && \
    rm -rf /var/lib/apt/lists/*

WORKDIR /Spare
COPY . .

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        build-essential \
        cmake \
        ninja-build \
        python3-dev \
        libopenblas-dev \
        liblapack-dev \
        libffi-dev \
        libatlas-base-dev \
        && rm -rf /var/lib/apt/lists/*

RUN python -m venv /Spare/venv && \
    /Spare/venv/bin/pip install --upgrade pip && \
    /Spare/venv/bin/pip install --no-cache-dir -r requirements.txt

ENV CONDA_DIR=/opt/conda

RUN curl -sSL https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
        -o /tmp/miniconda.sh && \
    bash /tmp/miniconda.sh -b -p ${CONDA_DIR} && \
    rm /tmp/miniconda.sh && \
    ${CONDA_DIR}/bin/conda clean -afy
ENV PATH=${CONDA_DIR}/bin:$PATH

SHELL ["/bin/bash", "-c"]
RUN chmod +777 ./run_scripts/*
CMD ./run_scripts/run_unqomp_baselines.sh && ./run_scripts/run_all.sh && exec bash