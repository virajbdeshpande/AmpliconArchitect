 # Use an official Python runtime as a parent image
FROM ubuntu

# Set the working directory to /app
WORKDIR /home/

#Copy AA and mosek to image
RUN mkdir -p /home/programs

#Download libraries for AA
RUN apt-get update && apt-get install -y
RUN apt-get install python-dev gfortran python-numpy python-scipy python-matplotlib python-pip zlib1g-dev samtools unzip -y
RUN pip install pysam

RUN mkdir -p /home/output/
RUN mkdir -p /home/input/
ADD run_aa_script.sh /home/

#Set environmental variables
ADD $MOSEK_DIR /home/programs/
RUN echo export MOSEKPLATFORM=linux64x86 >> ~/.bashrc
RUN export MOSEKPLATFORM=linux64x86
RUN echo export PATH=$PATH:/home/programs/mosek/8/tools/platform/$MOSEKPLATFORM/bin >> ~/.bashrc
RUN echo export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/programs/mosek/8/tools/platform/$MOSEKPLATFORM/bin >> ~/.bashrc
RUN echo export MOSEKLM_LICENSE_FILE=/home/programs/mosek/8/licenses >> ~/.bashrc
RUN cd /home/programs/mosek/8/tools/platform/linux64x86/python/2/ && python setup.py install
RUN echo export AA_DATA_REPO=/home/data_repo >> ~/.bashrc
ADD https://github.com/virajbdeshpande/AmpliconArchitect/archive/master.zip /home/programs
RUN cd /home/programs && unzip master.zip