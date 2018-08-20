 # Use an official Python runtime as a parent image
FROM ubuntu

# Set the working directory to /app
WORKDIR /home/

#Copy AA and mosek to image
RUN mkdir -p /home/programs
ADD AmpliconArchitect.v1.tar.gz /home/programs
ADD $MOSEK_DIR /home/programs/mosek

#Download libraries for AA
RUN apt-get update && apt-get install -y
RUN apt-get install python-dev gfortran python-numpy python-scipy python-matplotlib python-pip zlib1g-dev samtools -y

#Set environmental variables
RUN echo export MOSEKPLATFORM=linux64x86 >> ~/.bashrc
RUN export MOSEKPLATFORM=linux64x86
RUN echo export PATH=$PATH:/home/programs/mosek/8/tools/platform/$MOSEKPLATFORM/bin >> ~/.bashrc
RUN echo export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/programs/mosek/8/tools/platform/$MOSEKPLATFORM/bin >> ~/.bashrc
RUN echo export MOSEKLM_LICENSE_FILE=/home/programs/mosek/8/licenses >> ~/.bashrc
RUN cd /home/programs/mosek/8/tools/platform/linux64x86/python/2/ && python setup.py install
RUN echo export AA_DATA_REPO=/home/data_repo >> ~/.bashrc

RUN pip install pysam

RUN mkdir -p /home/output/
RUN mkdir -p /home/input/
ADD run_aa_script.sh /home/
