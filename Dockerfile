# Dockerfile — reproducible environment for the jet / EEC analysis (analysis + ML layer)
#
#   Build : docker build -t jet-eec .
#   Run   : docker run -it --rm -p 8888:8888 -v "$PWD":/work jet-eec
#   GPU   : docker run -it --rm --gpus all -p 8888:8888 -v "$PWD":/work jet-eec
#           (requires nvidia-container-toolkit on the host + a CUDA PyTorch build — see README)
#
# This container is the ANALYSIS + ML layer. For CMS AOD-tier data you additionally
# use the official CMSSW image (cmsopendata/cmssw_*) as a separate data-access layer.

FROM mambaorg/micromamba:1.5-jammy

# Solve the environment into base (the name: field in the yml is ignored on install)
COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/environment.yml
RUN micromamba install -y -n base -f /tmp/environment.yml && \
    micromamba clean --all --yes

# Activate base for subsequent RUN/CMD layers
ARG MAMBA_DOCKERFILE_ACTIVATE=1

WORKDIR /work
EXPOSE 8888

# Launch JupyterLab (token disabled for LOCAL use only — do not expose this port publicly)
CMD ["micromamba", "run", "-n", "base", "jupyter", "lab", \
     "--ip=0.0.0.0", "--port=8888", "--no-browser", "--allow-root", \
     "--ServerApp.token="]
