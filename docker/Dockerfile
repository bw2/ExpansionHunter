FROM python:3.9-slim-bullseye

RUN apt update && apt-get install --no-install-recommends -y \
    ca-certificates \
    pkg-config \
    dpkg-dev \
    build-essential \
    libssl-dev \
    libncurses5-dev \
    libcurl4-openssl-dev \
    libbz2-dev \
    liblzma-dev \
    zlib1g-dev \
    gcc \
    g++ \
    autoconf \
    make \
    unzip \
    bzip2 \
    procps \
    curl \
    wget \
    git \
    jq \
    time

ENV SAMTOOLS_VERSION="1.17"
RUN wget https://github.com/samtools/htslib/releases/download/${SAMTOOLS_VERSION}/htslib-${SAMTOOLS_VERSION}.tar.bz2 \
    && tar xjf htslib-${SAMTOOLS_VERSION}.tar.bz2 \
    && rm htslib-${SAMTOOLS_VERSION}.tar.bz2 \
    && cd htslib-${SAMTOOLS_VERSION} \
    && ./configure --enable-gcs \
    && make \
    && make install \
    && make clean

RUN wget -q https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2 \
    && tar xjf samtools-${SAMTOOLS_VERSION}.tar.bz2 \
    && rm samtools-${SAMTOOLS_VERSION}.tar.bz2 \
    && cd samtools-${SAMTOOLS_VERSION} \
    && ./configure \
    && make \
    && make install \
    && make clean

RUN python3 -m pip install --upgrade cmake pip --no-cache-dir

# install gcloud
ENV GCLOUD_SDK_VERISON="458.0.1"
RUN wget https://dl.google.com/dl/cloudsdk/channels/rapid/downloads/google-cloud-sdk-${GCLOUD_SDK_VERISON}-linux-x86_64.tar.gz \
    && tar xzf google-cloud-sdk-${GCLOUD_SDK_VERISON}-linux-x86_64.tar.gz \
    && /google-cloud-sdk/install.sh -q \
    && /google-cloud-sdk/bin/gcloud components update --quiet

ENV PATH=/google-cloud-sdk/bin:$PATH

RUN apt update && apt-get install --no-install-recommends -y libbz2-dev libcurl4-openssl-dev liblzma-dev libssl-dev zlib1g-dev libboost-all-dev  libgmock-dev
RUN apt update && apt-get install --no-install-recommends -y libgtest-dev libspdlog-dev libabsl-dev

# install REViewer
RUN wget -q https://github.com/bw2/REViewer/archive/refs/heads/master.zip \
    && unzip master.zip \
    && cd REViewer-master \
	&& mkdir build \
	&& cd build \
	&& cmake .. \
	&& make \
	&& cp install/bin/REViewer /usr/local/bin/ \
	&& cd / \
	&& rm -r REViewer-master master.zip

COPY ehunter .

RUN mkdir build \
    && cd build \
    && cmake .. \
    && make \
    && make install 

ARG CACHEBUST=1 
RUN git clone https://github.com/broadinstitute/str-analysis.git \
    && cd str-analysis \
    && python3 -m pip install . --no-cache-dir

ENTRYPOINT ["/bin/bash"]