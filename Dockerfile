FROM nfcore/base:1.7
LABEL authors="mpozud00@gmail.com" \
      description="Docker image containing all requirements for the mpozud00/demultiplexing pipeline"

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/mpozuelo-demultiplexing/bin:$PATH
