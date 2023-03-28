FROM ubuntu:focal
ENV TZ="America/New_York"
RUN apt update -y && apt install -y tzdata
RUN apt install -y --no-install-recommends \
    curl \
    git \
    apt-transport-https \
    ca-certificates \
    libegl1-mesa-dev


RUN curl -s -o /tmp/miniconda.sh https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
    && /bin/bash /tmp/miniconda.sh -b -p /opt/conda \
    && rm /tmp/miniconda.sh
RUN /opt/conda/bin/conda config --set always_yes yes \
    && /opt/conda/bin/conda config --add channels conda-forge \
    && /opt/conda/bin/conda install --freeze-installed \
        pip \
        setuptools \
        nomkl \
        sfepy=2022.2 \
        mayavi \
        scikit-umfpack \
        pytest \
        pyvista \
    && find /opt/conda/ -follow -type f -name '*.a' -delete \
    && find /opt/conda/ -follow -type f -name '*.pyc' -delete \
    && find /opt/conda/ -follow -type f -name '*.js.map' -delete \
    && /opt/conda/bin/conda clean -afy
RUN /opt/conda/bin/conda init bash \
    && . /root/.bashrc \
    && conda activate \
    && pip install tetgen pydantic
RUN mkdir -p /app/input
RUN mkdir -p /app/output
RUN mkdir -p /app/principalstresslines
RUN mkdir -p /app/testdata
COPY principalstresslines/ /app/principalstresslines/
COPY entrypoint.sh /app/entrypoint.sh
ENTRYPOINT [ "/bin/sh","/app/entrypoint.sh" ]