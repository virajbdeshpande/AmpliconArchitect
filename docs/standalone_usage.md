## Installing AA as a standalone tool

This is not recommended, as the AA module will be seperated from the AmpliconSuite-pipeline tools that prepare and filter its inputs.

### Installing from GitHub source code:


>`git clone https://github.com/jluebeck/AmpliconArchitect.git`

* Set `AmpliconArchitect/src` as `$AA_SRC`:
```bash
    cd AmpliconArchitect
    echo export AA_SRC=$PWD/src >> ~/.bashrc
```

#### Dependencies

1. Python 2.6+ or 3.5+
2. Ubuntu libraries and tools:
```bash
sudo apt-get install software-properties-common -y
sudo add-apt-repository universe -y
sudo apt-get update && sudo apt-get install -y
sudo apt-get install build-essential python-dev gfortran zlib1g-dev samtools wget -y

# below steps only required if using python2
sudo apt-get install python2
wget https://bootstrap.pypa.io/pip/2.7/get-pip.py
sudo python2 get-pip.py
```
3. **Python3 packages**: Note that [pysam](https://github.com/pysam-developers/pysam) verion 0.9.0  or higher is required. Flask is optional.

`pip3 install pysam Cython numpy scipy matplotlib future mosek Flask`

**... or for python 2:**

`pip2 install 'pysam==0.15.2' Cython numpy scipy matplotlib future mosek Flask` 

Note that 0.15.2 is the last version of pysam which appears to support pip2 installation, however AA itself supports the more recent versions.

4. Configure the Mosek optimization tool:
```bash
mkdir -p $HOME/mosek/
# Then please obtain license from https://www.mosek.com/products/academic-licenses/ or https://www.mosek.com/try/ and place in $HOME/mosek/
```
If you happen to be using the commerical version of the Mosek license (this is uncommon as Mosek is free for academic use), you will need the version which supports both PTON and PTS functions. 

5. (Optional) Arial font for matplotlib:

To get Microsoft fonts on Ubuntu:
```
sudo apt-get install fontconfig ttf-mscorefonts-installer
sudo fc-cache -f
```

### Standalone usage of AA
Standalone usage of AA requires many more manual steps than using [AmpliconSuite-pipeline](https://github.com/AmpliconSuite/AmpliconSuite-pipeline), and does not include best practices for seed region identification.
### Input data
AA requires 2 input files:

1. Coordinate-sorted, indexed BAM file:
    * Align reads to a reference present in the `data_repo`.
    * Recommended depth of coverage for WGS data is 5X-10X (higher is also fine but may run a bit more slowly).
    * Bamfile may be downsampled using `$AA_SRC/downsample.py` or when running AA with the option `--downsample` (default is `--downsample 10`. 
    * If sample has multiple reads groups with very different read lengths and fragment size distributions, then we recommend downsampling the bam file be selecting the read groups which have similar read lengths and fragment size distribution.
2. BED file with seed intervals:
    * We recommend generating this using [AmpliconSuite-pipeline](https://github.com/AmpliconSuite/AmpliconSuite-pipeline) as it automates all the steps described below.
    * AA has been tested on seed intervals generated as follows:
        - CNVs from CNV caller CNVkit, ReadDepth (with parameter file `$AA_SRC/src/read_depth_params`), or Canvas.
        - Prefilter the CNV calls using `cnv_prefilter.py` from AmpliconSuite-pipeline to check CNV calls against arm-level copy number.
        - Select CNVs with copy number > 4.5x and size > 50kbp (default), apply sequence context filters to filter repetitive and low complexity elements. Finally it merges adjacent CNVs into a single interval using:

            `python2 $AA_SRC/amplified_intervals.py --bed {bed file of cnv calls} --out {outFileNamePrefix} --bam {BamFileName} --ref {ref}`
        - **Note that these preprocessing steps are critical to AA as it removes low-mappability and low-complexity regions. AmpliconSuite-pipeline provides additional filters for karyotypic abnormalities not provided by `amplified_intervals.py` alone**.
        - Optional argument `--ref` should match the name of the folder in `data_repo` which corresponds to the version of human reference genome used in the BAM file.


### Example command
>`python $AA_SRC/AmpliconArchitect.py --bam {input_bam} --bed {seeds bed file} --out {prefix_of_output_files} <optional arguments>`

