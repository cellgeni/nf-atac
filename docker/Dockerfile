FROM python:3.11

ENV DEBIAN_FRONTEND=noninteractive
ENV VENV_PATH="/env"
ENV MALLET_PATH="/opt/Mallet"
ENV PATH="${VENV_PATH}/bin:$PATH"
ENV PATH="${MALLET_PATH}/bin:$PATH"

# Update and install system dependencies
RUN apt-get update && apt-get install -y \
    build-essential libcurl4-openssl-dev zlib1g-dev libfftw3-dev libc++-dev hdf5-tools \
    libbz2-dev liblzma-dev procps libfftw3-dev ant git

# Install uv package manager
RUN pip install uv

# Create a new environment with uv and install packages
RUN uv venv "${VENV_PATH}" && \
    uv pip install --no-cache-dir \
    setuptools wheel jupyterlab notebook papermill tqdm ipywidgets

# Install scenicplus with handling for datrie issues
RUN git clone https://github.com/aertslab/scenicplus /opt/scenicplus && \
    uv pip install --no-cache-dir --only-binary=datrie -r /opt/scenicplus/requirements.txt || \
    (sed -i '/^datrie/d' /opt/scenicplus/requirements.txt && uv pip install -r /opt/scenicplus/requirements.txt --no-cache-dir) && \
    uv pip install /opt/scenicplus

# Install Mallet
RUN mkdir -p /usr/share/man/man1 && \
    git clone --depth=1 https://github.com/mimno/Mallet.git "${MALLET_PATH}" && \
    cd "${MALLET_PATH}" && ant

COPY Dockerfile /docker/
RUN chmod -R 755 /docker
RUN chmod -R 755 /docker