FROM debian:buster-slim AS base

ENV LANG=C.UTF-8 LC_ALL=C.UTF-8
ENV PATH /opt/conda/bin:$PATH

# Instalacao do Gfortran 7 e Miniconda
RUN set -x && \
    apt-get update --fix-missing && \
    apt-get install -y --no-install-recommends \
    gcc \
    gfortran-7 \
    ca-certificates \
    git \
    ssh \
    wget \
    rsync \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir -p /opt/conda \
    && wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /opt/conda/miniconda.sh \
    && bash /opt/conda/miniconda.sh -b -u -p /opt/conda \
    && rm -rf /opt/conda/miniconda.sh \
    && /opt/conda/bin/conda init bash \
    && find /opt/conda/ -follow -type f -name '*.a' -delete \
    && find /opt/conda/ -follow -type f -name '*.js.map' -delete \
    && /opt/conda/bin/conda clean -afy \
    && chmod =2775 /opt/conda \
    && ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh

# Create Groups ton, conda
# Grupo ton existe no ambiente do linea
RUN groupadd -r ton --gid 15010 \
    && groupadd -r conda --gid 900


# -------------- Python 3.8 Environment Stage --------------
FROM base AS py3_build
COPY ./hello_world/environment.yaml .
RUN /bin/bash --login -c "conda init bash \
    && source ~/.bashrc \
    && conda env create -f environment.yaml \
    && conda activate pipe_hello_world \
    && rm -rf environment.yaml"

# -------------- Runtime Stage --------------
FROM base

# If this is set to a non-empty string, Python won’t try
# to write .pyc files on the import of source modules
ENV PYTHONDONTWRITEBYTECODE=1

# Force the stdout and stderr streams to be unbuffered.
# This option has no effect on the stdin stream.
ENV PYTHONUNBUFFERED=1

ENV BSP_PLANETARY_NAME=de440.bsp
ENV LEAP_SECOND_NAME=naif0012.tls

# # Download da BSP planetary
# # OBS. o Download demora bastante!
# RUN wget --no-verbose --show-progress \
#     --progress=bar:force:noscroll \
#     https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/$BSP_PLANETARY_NAME

# # Download Leap Second
# RUN wget --no-verbose --show-progress \
#     --progress=bar:force:noscroll \
#     https://naif.jpl.nasa.gov/pub/naif/generic_kernels/lsk/$LEAP_SECOND_NAME

# Python 3.10 environment
COPY --chown=:conda --chmod=775 --from=py3_build /opt/conda/envs/pipe_hello_world /opt/conda/envs/pipe_hello_world

RUN set -x && \
    apt-get update --fix-missing && \
    apt-get install -y --no-install-recommends \
    bash-completion \
    nano \
    vim \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

ARG APP_HOME=/app
ARG USERNAME=app.tno
ARG USERUID=1000
ARG USERGID=1000

# Create NonRoot user and add to groups
RUN groupadd --gid ${USERGID} ${USERNAME}  \
    && useradd --uid ${USERUID} --gid ${USERGID} --shell /bin/bash --create-home ${USERNAME} \
    && usermod -a -G 15010,900 ${USERNAME}

ENV CONDA_EXE=/opt/conda/bin

COPY --chown=${USERNAME}:ton --chmod=775 . /app
# COPY --chown=${USERNAME}:ton --chmod=775 entrypoint.sh ${APP_HOME}
# COPY --chown=${USERNAME}:ton --chmod=775 check_enviroment.py ${APP_HOME}

# Copia o enviroment utilizado dentro do container para atualizar ou criar o enviroment fora do container.
# COPY --chown=${USERNAME}:ton --chmod=2775 environment.py3.yml ${APP_HOME}

# RUN chmod -R =2775 ${PIPELINE_ROOT} \
#     && chmod =2775 ${APP_HOME} \
#     && mv de440.bsp ${PIPELINE_PATH}/$BSP_PLANETARY_NAME \
#     && mv naif0012.tls ${PIPELINE_PATH}/$LEAP_SECOND_NAME

WORKDIR ${APP_HOME}


USER ${USERNAME}
RUN /bin/bash --login -c "conda init bash"
ENV PATH=${PATH}:/home/${USERNAME}/.local/bin

# ENTRYPOINT [ "./entrypoint.sh" ]
