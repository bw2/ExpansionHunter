# syntax=docker/dockerfile:1
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

ENV SAMTOOLS_VERSION="1.21"
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
ENV GCLOUD_SDK_VERISON="512.0.0"
RUN wget https://dl.google.com/dl/cloudsdk/channels/rapid/downloads/google-cloud-sdk-${GCLOUD_SDK_VERISON}-linux-x86_64.tar.gz \
    && tar xzf google-cloud-sdk-${GCLOUD_SDK_VERISON}-linux-x86_64.tar.gz \
    && /google-cloud-sdk/install.sh -q \
    && /google-cloud-sdk/bin/gcloud components update --quiet

ENV PATH=/google-cloud-sdk/bin:$PATH

RUN apt update && apt-get install --no-install-recommends -y libbz2-dev libcurl4-openssl-dev liblzma-dev libssl-dev zlib1g-dev libgmock-dev
RUN apt update && apt-get install --no-install-recommends -y libgtest-dev libspdlog-dev

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

# install mimalloc (faster allocator; ~30% faster genotyping with lower peak RSS, byte-identical output).
# Built from source so its CMake config is installed and ExpansionHunter's find_package(mimalloc) links it.
ENV MIMALLOC_VERSION="2.1.7"
RUN wget -q https://github.com/microsoft/mimalloc/archive/refs/tags/v${MIMALLOC_VERSION}.tar.gz \
    && tar xzf v${MIMALLOC_VERSION}.tar.gz \
    && rm v${MIMALLOC_VERSION}.tar.gz \
    && cd mimalloc-${MIMALLOC_VERSION} \
    && mkdir build \
    && cd build \
    && cmake .. -DCMAKE_BUILD_TYPE=Release -DMI_BUILD_TESTS=OFF \
    && make -j$(nproc) \
    && make install \
    && cd / \
    && rm -r mimalloc-${MIMALLOC_VERSION}

COPY ehunter ehunter
COPY CMakeLists.txt CMakeLists.txt

# Bind-mount the host .git (read-only, BuildKit-only, never baked into an image layer) so cmake's
# `git rev-parse`/`git log` calls (ehunter/CMakeLists.txt) can embed the real commit SHA/date into the
# "Version" field of the output JSON, instead of falling back to "unknown". Requires BuildKit (the
# default builder since Docker 23+). The bind-mounted .git keeps the host user's UID while this RUN
# executes as root, which trips git's "dubious ownership" check (CVE-2022-24765) unless it's told to
# trust it, so add safe.directory before running cmake.
RUN --mount=type=bind,source=.git,target=/.git,ro \
    git config --global --add safe.directory '*' \
    && mkdir build \
    && cd build \
    && cmake .. \
    && make \
    && make install \
    && cp /build/install/bin/ExpansionHunter /usr/bin/ \
    && make clean

ARG CACHEBUST=1 
RUN git clone https://github.com/broadinstitute/str-analysis.git \
    && cd str-analysis \
    && python3 -m pip install . --no-cache-dir

ENTRYPOINT ["/bin/bash"]