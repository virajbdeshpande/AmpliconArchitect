# Use an official Python runtime as a parent image
FROM ubuntu:20.04

# Build in non-interactive mode for online continuous building
ARG DEBIAN_FRONTEND=noninteractive

# Set the working directory to /app
WORKDIR /home/

# Make "requirements" directory and copy requirements over
RUN mkdir -p /home/requirements
COPY requirements/* /home/requirements/

#Copy AA and mosek to image
RUN mkdir -p /home/programs

#Download libraries for AA
RUN apt-get update
RUN apt-get install -f software-properties-common=0.99.9.8 -y
RUN add-apt-repository universe -y
RUN apt-get install -y python2=2.7.17-2ubuntu4
ADD https://bootstrap.pypa.io/pip/2.7/get-pip.py /home/programs
RUN python2 /home/programs/get-pip.py
RUN pip2 --version
RUN apt-get update && apt-get install -y
#RUN apt-get install libcurl4-gnutls-dev libssl-dev python-dev libbz2-dev liblzma-dev gfortran zlib1g-dev samtools wget unzip curl -y
RUN apt-get install -y --fix-missing \
bcftools=1.10.2-2 \
bwa=0.7.17-4 \
fontconfig=2.13.1-2ubuntu3 \
gfortran=4:9.3.0-1ubuntu2 \
libbz2-dev=1.0.8-2 \
liblzma-dev \
python-dev \
python3-dev=3.8.2-0ubuntu2 \
samtools=1.10-3 \
ttf-mscorefonts-installer=3.7ubuntu6 \
unzip=6.0-25ubuntu1 \
wget=1.20.3-1ubuntu2 \
zlib1g-dev
RUN fc-cache -f

RUN pip2 install -r /home/requirements/pip2_requirements.txt

# To install Mosek8 (Mosek 9 is installed in the pip2 requirements).
#RUN cd /home/programs && wget http://download.mosek.com/stable/8.0.0.60/mosektoolslinux64x86.tar.bz2
#RUN cd /home/programs && tar xf mosektoolslinux64x86.tar.bz2
#RUN echo export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/programs/mosek/8/tools/platform/$MOSEKPLATFORM/bin >> ~/.bashrc
# RUN echo export MOSEKPLATFORM=linux64x86 >> ~/.bashrc
# RUN echo export PATH=$PATH:/home/programs/mosek/8/tools/platform/$MOSEKPLATFORM/bin >> ~/.bashrc
# RUN cd /home/programs/mosek/8/tools/platform/linux64x86/python/2/ && python2 setup.py install


RUN mkdir -p /home/output/
RUN mkdir -p /home/input/
RUN mkdir -p /home/mosek/
ADD run_aa_script.sh /home/

#Set environmental variables
RUN echo export AA_DATA_REPO=/home/data_repo >> ~/.bashrc
ADD https://github.com/jluebeck/AmpliconArchitect/archive/master.zip /home/programs
RUN cd /home/programs && unzip master.zip
