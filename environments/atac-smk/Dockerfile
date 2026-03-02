FROM ubuntu:22.04

# Use bash
SHELL ["/bin/bash", "-c"]

# Environment variables
ENV PYTHONUNBUFFERED=1
ENV DEBIAN_FRONTEND=noninteractive
ENV BCFTOOLS_PLUGINS=/usr/local/bcftools/plugins
ENV PATH=/root/.local/bin:${PATH}
RUN echo -e $PATH | tr ":" "\n"

RUN apt-get -y update && \
    apt-get -y install \
        gcc-12 \
        g++-12 \
        rsync \
        python3-dev \
        python3-full \
        build-essential \
        mbuffer \
        curl \
        git \
        wget \
        lz4 \
        bc \
        jq \
        unzip \
        cmake \
        curl \
        pkg-config \
        autoconf \
        zlib1g-dev \
        libbz2-dev \
        liblzma-dev \
        libcurl4-openssl-dev \
        libssl-dev \
        lbzip2 \
        libexpat1-dev \
        libncurses5-dev \
        libfreetype6-dev \
        libfontconfig1-dev
RUN update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-12 100 && \
    update-alternatives --install /usr/bin/g++ g++ /usr/bin/g++-12 100

# Clean up
RUN apt-get -y clean && \
    apt-get -y autoremove && \
    rm -rf /var/lib/apt-get/lists/*

RUN echo -e "GCC Version: $(gcc --version | head -n 1)"

# Install python deps with uv
ENV VIRTUAL_ENV=/python-venv
RUN curl -LsSf https://astral.sh/uv/install.sh | sh
RUN uv venv ${VIRTUAL_ENV}
ENV PATH=${VIRTUAL_ENV}/bin:${PATH}
# Install all optional dependencies for polars EXCEPT deltalake and apache iceberg
RUN uv pip install --no-cache-dir \
    macs3 \
    beautifulsoup4 \
    pandas \
    polars \
    oxbow \
    pysam \
    numpy \
    matplotlib \
    tqdm \
    multiqc \
    bioframe

# Install bowtie2 from source
RUN git clone --recursive https://github.com/BenLangmead/bowtie2.git /opt/bowtie2 && \
    cd /opt/bowtie2 && \
    git submodule init && git submodule update && \
    make libsais -j$(nproc) USE_SAIS_OPENMP=1 && \
    make bowtie2-build-s -j$(nproc) USE_SAIS_OPENMP=1 && \
    make -j$(nproc) && \
    cd /usr/local/bin && \
    ln -s /opt/bowtie2/bowtie2

# Environment variables (compiler)
# NOTE: These must come after python installs, as they can interfere with certain packages (i.e. polars)
ENV CFLAGS="-O3 -march=x86-64-v4"
ENV CXXFLAGS="-O3 -march=x86-64-v4"

# Install fastp
RUN wget http://opengene.org/fastp/fastp && \
    mv fastp /usr/local/bin/ && \
    chmod a+x /usr/local/bin/fastp

# Install pigz (parallel gzip implementation)
RUN curl -sSLf https://zlib.net/pigz/pigz.tar.gz | tar xzf - && \
    cd /pigz && \
    make && \
    cp pigz unpigz /usr/local/bin
RUN rm -rf /pigz

# Install samtools stuff
RUN curl -s https://api.github.com/repos/samtools/htslib/releases/latest | \
    jq -r '.assets[] | select(.name | test("htslib-.*\\.tar\\.bz2")) | .browser_download_url' | \
    xargs curl -fsSL | tar xjf - && \
    cd /htslib-* && \
    ./configure --prefix=/usr && \
    make -j$(nproc) && \
    make install
RUN curl -s https://api.github.com/repos/samtools/samtools/releases/latest | \
    jq -r '.assets[] | select(.name | test("samtools-.*\\.tar\\.bz2")) | .browser_download_url' | \
    xargs curl -fsSL | tar xjf - && \
    cd /samtools-* && \
    ./configure && \
    make -j$(nproc) && \
    make install
RUN rm -rf /samtools-* /htslib-*

# Install bedtools
RUN curl -s https://api.github.com/repos/arq5x/bedtools2/releases/latest | \
    jq -r '.assets[] | select(.name | test("bedtools-.*\\.tar\\.gz")) | .browser_download_url' | \
    xargs curl -fsSL | tar xzf - && \
    cd bedtools2 && \
    make static && \
    mv bin/bedtools.static /usr/local/bin/bedtools
RUN rm -rf bedtools2

# Download ucsc tools
RUN rsync -aP rsync://hgdownload.soe.ucsc.edu/genome/admin/exe/linux.x86_64/bedClip /usr/local/bin/bedClip && \
    chmod a+x /usr/local/bin/bedClip
RUN rsync -aP rsync://hgdownload.soe.ucsc.edu/genome/admin/exe/linux.x86_64/bedGraphToBigWig /usr/local/bin/bedGraphToBigWig && \
    chmod a+x /usr/local/bin/bedGraphToBigWig

# Download static sambamba binary
RUN curl -fSsL https://github.com/biod/sambamba/releases/download/v1.0.1/sambamba-1.0.1-linux-amd64-static.gz | unpigz > /usr/local/bin/sambamba && \
    chmod a+x /usr/local/bin/sambamba

# Test installations
RUN which samtools \
    which bedtools \
    which bowtie2 \
    which bedClip \
    which bedGraphToBigWig \
    which sambamba

CMD ["/bin/bash"]
