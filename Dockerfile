FROM rockylinux/rockylinux:9.5-minimal AS base

ENV LANG=C.UTF-8 LC_ALL=C.UTF-8
ENV PATH /opt/conda/bin:$PATH
ENV SHELL /bin/zsh

# Instalacao das dependências no Rocky Linux Minimal
RUN set -x && \
    microdnf update -y && \
    microdnf install -y \
    findutils \
    tar \
    gzip \
    gcc \
    gcc-gfortran \
    ca-certificates \
    make \
    libpq-devel \
    git \
    openssh-clients \
    wget \
    rsync \
    which \
    zsh \
    curl \
    && microdnf clean all \
    && rm -rf /var/cache/dnf \
    && mkdir -p /opt/conda \
    && wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /opt/conda/miniconda.sh \
    && bash /opt/conda/miniconda.sh -b -u -p /opt/conda \
    && rm -rf /opt/conda/miniconda.sh \
    && /opt/conda/bin/conda init bash \
    && /opt/conda/bin/conda init zsh \
    && find /opt/conda/ -follow -type f -name '*.a' -delete \
    && find /opt/conda/ -follow -type f -name '*.js.map' -delete \
    && /opt/conda/bin/conda clean -afy \
    && chmod =2775 /opt/conda \
    && ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh

# Aceitar os Termos de Serviço do Conda
RUN /opt/conda/bin/conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main && \
    /opt/conda/bin/conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r

# Create Groups ton, conda
RUN groupadd -r ton --gid 15010 \
    && groupadd -r conda --gid 900

# -------------- PRAIA OCC compile Stage --------------
FROM base AS praia_occ
ADD praia_occ_src /tmp/praia_occ_src
RUN mkdir /tmp/praia_occ \
    && cd /tmp/praia_occ_src \
    && gfortran geradata.f -o geradata spicelib.a \
    && mv geradata /tmp/praia_occ \
    && gfortran elimina.f -o elimina \
    && mv elimina /tmp/praia_occ \
    && gfortran PRAIA_occ_star_search_12.f -o PRAIA_occ_star_search_12 \
    && mv PRAIA_occ_star_search_12 /tmp/praia_occ \
    && cd ~/ \
    && rm -r /tmp/praia_occ_src

# -------------- Python 3.8 Environment Stage --------------
FROM base AS py3_build

COPY ./predict_occultation/environment.yaml .

# Aceitar TOS e criar environment
RUN /bin/bash --login -c "conda init bash \
    && source ~/.bashrc \
    && conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main \
    && conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r \
    && conda env create -f environment.yaml \
    && conda activate pipe_pred_occ \
    && rm -rf environment.yaml"

# -------------- Runtime Stage --------------
FROM base

ENV PYTHONDONTWRITEBYTECODE=1
ENV PYTHONUNBUFFERED=1
ENV SHELL /bin/zsh

ARG APP_HOME=/app
ARG USERNAME=app.tno
ARG USERUID=1000
ARG USERGID=1000
ARG BSP_PLANETARY_NAME=de440.bsp
ARG LEAP_SECOND_NAME=naif0012.tls

WORKDIR ${APP_HOME}

# Download da BSP planetary e Leap Second
RUN wget --no-verbose --show-progress \
    --progress=bar:force:noscroll \
    https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/${BSP_PLANETARY_NAME} \
    -O ${APP_HOME}/${BSP_PLANETARY_NAME} \
    && wget --no-verbose --show-progress \
    --progress=bar:force:noscroll \
    https://naif.jpl.nasa.gov/pub/naif/generic_kernels/lsk/${LEAP_SECOND_NAME} \
    -O ${APP_HOME}/${LEAP_SECOND_NAME}

# Python environment
COPY --chown=:conda --chmod=775 --from=py3_build /opt/conda/envs/pipe_pred_occ /opt/conda/envs/pipe_pred_occ

# PRAIA OCC binaries
COPY --from=praia_occ /tmp/praia_occ/* /usr/local/bin

# Instalar ferramentas adicionais (usando microdnf na minimal)
# Incluindo dependências necessárias para VS Code Server
RUN microdnf install -y \
    findutils \
    tar \
    gzip \
    bash-completion \
    nano \
    vim \
    zsh \
    curl \
    procps-ng \
    lsof \               
    shadow-utils \   
    git \    
    && microdnf clean all

# Create NonRoot user and add to groups
RUN groupadd --gid ${USERGID} ${USERNAME}  \
    && useradd --uid ${USERUID} --gid ${USERGID} --shell /bin/zsh --create-home ${USERNAME} \
    && usermod -a -G 15010,900 ${USERNAME}

# Instalar Oh My Zsh para o usuário
USER ${USERNAME}

# Instalar Oh My Zsh de forma não interativa
RUN sh -c "$(curl -fsSL https://raw.githubusercontent.com/ohmyzsh/ohmyzsh/master/tools/install.sh)" "" --unattended

# Instalar plugins do Zsh
RUN git clone https://github.com/zsh-users/zsh-syntax-highlighting.git ${HOME}/.oh-my-zsh/custom/plugins/zsh-syntax-highlighting \
    && git clone https://github.com/zsh-users/zsh-autosuggestions ${HOME}/.oh-my-zsh/custom/plugins/zsh-autosuggestions

# Configurar Conda no Zsh ANTES de modificar o .zshrc
RUN /bin/zsh -c "conda init zsh"

# Usar um tema que não requer powerline (robbyrussell é o padrão e funciona bem)
RUN sed -i 's/^ZSH_THEME=.*/ZSH_THEME="robbyrussell"/' ${HOME}/.zshrc

# Configurar plugins
RUN sed -i 's/^plugins=.*/plugins=(git conda zsh-syntax-highlighting zsh-autosuggestions)/' ${HOME}/.zshrc

# Configurações adicionais para melhorar a experiência do Zsh
RUN echo 'export TERM="xterm-256color"' >> ${HOME}/.zshrc \
    && echo 'export EDITOR=nano' >> ${HOME}/.zshrc \
    && echo 'export VISUAL=nano' >> ${HOME}/.zshrc \
    && echo 'alias ll="ls -alh"' >> ${HOME}/.zshrc \
    && echo 'alias py="python"' >> ${HOME}/.zshrc \
    && echo 'alias cls="clear"' >> ${HOME}/.zshrc

# Configurar Conda para não ativar automaticamente a base
RUN echo 'export CONDA_AUTO_ACTIVATE_BASE=false' >> ${HOME}/.zshrc

USER root

ENV CONDA_EXE=/opt/conda/bin

COPY --chown=${USERNAME}:ton --chmod=775 . /app

RUN mkdir /data \
    && chown -R ${USERUID}:${USERGID} /data \
    && chmod -R g+w /data

USER ${USERNAME}

# Configurar Conda e instalar pacotes Python - usando bash para evitar problemas com zshrc
RUN /bin/bash -c "source /opt/conda/bin/activate \
    && conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main \
    && conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r \
    && conda activate base \
    && pip install pydantic PyYaml black isort"

ENV PATH=${PATH}:/home/${USERNAME}/.local/bin

# Definir Zsh como shell padrão
SHELL ["/bin/zsh", "-c"]

# Entrypoint padrão com Zsh
CMD ["zsh"]