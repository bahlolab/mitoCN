FROM mambaorg/micromamba:1.5.8
COPY --chown=$MAMBA_USER:$MAMBA_USER env.yaml /tmp/env.yaml
RUN micromamba install -y -n base -f /tmp/env.yaml && micromamba clean --all --yes
WORKDIR /app
COPY mitoCN.sh mtDNA_CN.R .
COPY hg19 hg19
COPY hg38 hg38
ENTRYPOINT ["/usr/local/bin/_entrypoint.sh", "./mitoCN.sh", "-d", "/app"]
