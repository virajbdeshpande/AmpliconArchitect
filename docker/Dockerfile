# Use an official Python runtime as a parent image
FROM ubuntu:20.04

# Build in non-interactive mode for online continuous building
ARG DEBIAN_FRONTEND=noninteractive

# Set the working directory to /app
WORKDIR /home/

#Copy AA and mosek to image
RUN mkdir -p /home/programs

#Download libraries for AA
RUN apt-get update
RUN apt-get install software-properties-common -y
RUN add-apt-repository universe -y
RUN apt-get install -y python2
ADD https://bootstrap.pypa.io/pip/2.7/get-pip.py /home/programs
RUN python2 /home/programs/get-pip.py
RUN pip2 --version
RUN apt-get update && apt-get install -y
RUN apt-get install libcurl4-gnutls-dev libssl-dev python-dev libbz2-dev liblzma-dev gfortran zlib1g-dev samtools wget unzip curl -y
#RUN update-alternatives --install /usr/bin/python python /usr/bin/python2 10
#RUN update-alternatives --config python
RUN pip2 install Cython
RUN pip2 install future numpy scipy matplotlib pysam==0.15.2 Flask

RUN cd /home/programs && wget http://download.mosek.com/stable/8.0.0.60/mosektoolslinux64x86.tar.bz2
RUN cd /home/programs && tar xf mosektoolslinux64x86.tar.bz2
# ADD mosek.lic /home/programs/mosek/8/licenses/mosek.lic

RUN mkdir -p /home/output/
RUN mkdir -p /home/input/
RUN mkdir -p /home/programs/mosek/8/licenses/
ADD run_aa_script.sh /home/

#Set environmental variables

RUN echo export MOSEKPLATFORM=linux64x86 >> ~/.bashrc
RUN export MOSEKPLATFORM=linux64x86
RUN echo export PATH=$PATH:/home/programs/mosek/8/tools/platform/$MOSEKPLATFORM/bin >> ~/.bashrc
RUN echo export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/programs/mosek/8/tools/platform/$MOSEKPLATFORM/bin >> ~/.bashrc
RUN echo export MOSEKLM_LICENSE_FILE=/home/programs/mosek/8/licenses >> ~/.bashrc
RUN cd /home/programs/mosek/8/tools/platform/linux64x86/python/2/ && python2 setup.py install
RUN echo export AA_DATA_REPO=/home/data_repo >> ~/.bashrc
ADD https://github.com/jluebeck/AmpliconArchitect/archive/master.zip /home/programs
RUN cd /home/programs && unzip master.zip
