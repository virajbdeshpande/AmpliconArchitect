# AmpliconArchitect Docker
We provide instructions below in order to build/run the AmpliconArchitect Docker.

# Running Docker:
1. Install Docker:
  * `https://docs.docker.com/install/`

2. Pull image:
  * `docker pull virajbdeshpande/ampliconarchitect`

3. Set environment variables AA_DATA_REPO and MOSEKLM_LICENSE_FILE as described in the README.md for AmpliconArchitect

4. Use run_aa_docker.sh to run AmpliconArchitect
  * `bash run_aa_docker.sh <OPTIONS>`
  * where `OPTIONS` are all options in the same format as provided to AmpliconArchitect.py

# Building Docker image from scratch:
`docker build . -t aa`
