FROM python:3.11

ENV DEBIAN_FRONTEND=noninteractive
ENV VENV_PATH="/env"
ENV MALLET_PATH="/opt/Mallet"
ENV PATH="${VENV_PATH}/bin:$PATH"
ENV PATH="${MALLET_PATH}/bin:$PATH"

RUN apt-get update && apt-get install -y \
    build-essential libcurl4-openssl-dev zlib1g-dev libfftw3-dev libc++-dev hdf5-tools \
    openjdk-17-jdk-headless libbz2-dev liblzma-dev procps libfftw3-dev ant

RUN python -m venv "${VENV_PATH}" && \
    . "${VENV_PATH}/bin/activate" && \
    pip install -U pip setuptools wheel && \
    pip install jupyterlab papermill

RUN git clone https://github.com/aertslab/scenicplus /opt/scenicplus && \
    . "${VENV_PATH}/bin/activate" && \
    pip install -r /opt/scenicplus/requirements.txt --no-cache-dir && \
    pip install /opt/scenicplus

RUN mkdir -p /usr/share/man/man1 && \
    git clone --depth=1 https://github.com/mimno/Mallet.git "${MALLET_PATH}" && \
    cd "${MALLET_PATH}" && ant

COPY Dockerfile /docker/
RUN chmod -R 755 /docker
RUN chmod -R 755 /docker