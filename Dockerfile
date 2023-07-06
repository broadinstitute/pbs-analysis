ARG R_VERSION=4.0.5

### BigWig tools

FROM gcr.io/cloud-builders/wget AS bigwig

COPY scripts/install_bigwig.sh .

RUN ./install_bigwig.sh /opt


### IGVTools

FROM gcr.io/cloud-builders/curl AS igvtools

WORKDIR /tmp

ENV VERSION="2.8.2"

RUN apt-get update -qq && \
    apt-get install -qq --no-install-recommends jq unzip && \
    curl -so IGV.zip "https://data.broadinstitute.org/igv/projects/downloads/${VERSION%.*}/IGV_${VERSION}.zip" && \
    unzip -q IGV.zip && \
    cd IGV_${VERSION} && \
    mv igvtools igv.args lib /opt


### R libraries

FROM r-base:${R_VERSION} AS r

RUN apt-get update -qq && \
    apt-get install -qq --no-install-recommends \
      libcurl4-openssl-dev \
      libssl-dev \
      libxml2-dev \
      zlib1g-dev

COPY scripts/install.R .

RUN ./install.R


### Final image

FROM r-base:${R_VERSION}

WORKDIR /

RUN apt-get update -qq && \
    apt-get install -qq --no-install-recommends \
      libxml2 \
      openjdk-17-jre-headless \
      samtools \
    && rm -rf /var/lib/apt/lists/*

COPY --from=bigwig /opt /bin
COPY --from=igvtools /opt /bin
COPY --from=r /usr/local/lib/R /usr/local/lib/R
COPY . .

ENTRYPOINT []
