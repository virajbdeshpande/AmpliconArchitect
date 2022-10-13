# AmpliconArchitect (AA)
![GitHub release (latest by date)](https://img.shields.io/github/v/release/jluebeck/AmpliconArchitect)
![GitHub commits since tagged version](https://img.shields.io/github/commits-since/jluebeck/AmpliconArchitect/v1.3.r1/master)
![GitHub all releases](https://img.shields.io/github/downloads/jluebeck/AmpliconArchitect/total)


### Recent updates:

### October 2022 update:
Version `1.3.r2` sets a deterministic seed for read downsampling, which can be disabled by setting `--random_seed`.

### April 2022 update:
Version `1.3.r1` released which supports both `python3` and `python2`. 
This release also corrects a bug in the visualization of amplicon copy numbers for 
genome segments < 50kbp and includes small improvements to seed filtering and 
merging of nearby segments. Complete changelog available [here](https://docs.google.com/document/d/1tGwKcHgaU36OzMZ02S5ky5JR-67WvRtXI3-Hv0FYkm4/edit?usp=sharing). 

### January 2022 update:
We have released a testing version of the mm10 mouse genome data repo [here](https://datasets.genepattern.org/?prefix=data/module_support_files/AmpliconArchitect/). 
To use with an existing installation please extract and place the mm10 directory with the other reference build directories in your AA data repo.
Upstream and downstream tools (AmpliconSuite-pipeline, AmpliconClassifier, CycleViz) are also enabled to accept the `--ref mm10`
argument. 


**[Older update descriptions are available here.](https://docs.google.com/document/d/1jqnCs46hrpYGBGrZQFop31ezskyludxNJEQdZONwFdc/edit?usp=sharing)**

## Introduction
Focal oncogene amplification and rearrangements drive tumor growth and evolution in multiple cancer types. Proposed mechanisms for focal amplification include extrachromosomal DNA (ecDNA) formation, breakage-fusion-bridge (BFB) mechanism, tandem duplications, chromothripsis and others.
Focally amplified regions are often hotspots for genomic rearrangements. As a result, the focally amplified region may undergo rapid copy number changes and the structure of the focally amplified region may evolve over time contributing to tumor progression. 
Furthermore, ecDNA elements may reintegrate back into the genome to form HSRs. The inter-cell heterogeneity in copy number of ecDNA as well as the interchangeability between ecDNA and HSR may allow the tumor to adapt to changing environment, e.g. targetted drug application. As a result, understanding the architecture of the focal amplifications is important to gain insights into cancer biology. AmpliconArchitect (AA) is a tool which can reconstruct the structure of focally amplified regions in a cancer sample using whole genome sequence short paired-end data.

Please check out the **detailed guide** on running AA [available here](https://github.com/jluebeck/PrepareAA/blob/master/GUIDE.md) to learn about best practices and see some FAQs.

**AmpliconArchitect was originally developed by Viraj Deshpande**, and is maintained by Jens Luebeck, Viraj Deshpande, and others in Vineet Bafna's lab. A full description of the method can be found in the following manuscript. Please cite the following reference if using AmpliconArchitect in your work:

*Deshpande, V. et al. Exploring the landscape of focal amplifications in cancer using AmpliconArchitect. Nat. Commun. 10, 392 (2019).* [(Article)](https://www.nature.com/articles/s41467-018-08200-y)

## Table of contents:
1. [Quickstart](#quickstart)
2. [Installation](#installation)
3. [Usage](#ampliconarchitect-reconstruction)
4. [The AA Algorithm](#the-aa-algorithm)
5. [File formats](#file-formats)
6. [Checkpointing and modular integration with other tools](#checkpointing-and-modular-integration-with-other-tools)

## Quickstart:
#### Prerequisites:
1. Docker (only if intending to use docker image):
    * Install docker: `https://docs.docker.com/install/`
    * (Optional): Add user to the docker group and relogin:
        `sudo usermod -a -G docker $USER`
2. License for Mosek optimization tool:
    * Obtain license file `mosek.lic` (`https://www.mosek.com/products/academic-licenses/` or `https://www.mosek.com/try/`)
    * `export MOSEKLM_LICENSE_FILE=<Parent directory of mosek.lic> >> ~/.bashrc && source ~/.bashrc`
3. Download AA data repositories and set environment variable AA_DATA_REPO:
    * Download from ` https://datasets.genepattern.org/?prefix=data/module_support_files/AmpliconArchitect/`

#### Usage:

`$AA --bam {input_bam} --bed {bed file} --out {prefix_of_output_files} <optional arguments>`

The execution script `$AA` may point to `AA=AmpliconArchitect/docker/run_aa_docker.sh`.

See instructions below for manual installation without docker.

## Installation
AA can be installed in 2 ways:
1. Docker image: This will automatically pull the latest build including necessary dependencies.
2. Github source code.

### Option 1: Docker image:
#### Prerequisites:
1. Docker:
    * Install docker: `https://docs.docker.com/install/`
    * (Optional): Add user to the docker group and relogin:
        `sudo usermod -a -G docker $USER`
2. License for Mosek optimization tool:
    * Obtain license file `mosek.lic` (`https://www.mosek.com/products/academic-licenses/` or `https://www.mosek.com/try/`)
    * `export MOSEKLM_LICENSE_FILE=<Parent directory of mosek.lic> >> ~/.bashrc && source ~/.bashrc`
3. Download AA data repositories and set environment variable AA_DATA_REPO:
    * Download from ` https://datasets.genepattern.org/?prefix=data/module_support_files/AmpliconArchitect/`

    * Set enviroment variable AA_DATA_REPO to point to the data_repo directory:
        ```bash
        mkdir data_repo && cd data_repo
        # copy or download files into data_repo directory
        tar zxf [hg19/GRCh37/GRCh38/mm10].tar.gz
        echo export AA_DATA_REPO=$PWD >> ~/.bashrc
        touch coverage.stats && chmod a+r coverage.stats
        source ~/.bashrc
        ```
#### Obtain AmpliconArchitect image and execution script:
1. Pull docker image:
    * `docker pull jluebeck/ampliconarchitect`

2. Clone script `run_aa_docker.sh` from Github:
    * `git clone https://github.com/jluebeck/AmpliconArchitect.git`

### Option 2: Github source code:
`git clone https://github.com/jluebeck/AmpliconArchitect.git`

**Note: In the rest of this document, we will refer to the path of the parent directory `AmpliconArchitect/src` as `$AA_SRC`**

```
cd AmpliconArchitect
echo export AA_SRC=$PWD/src >> ~/.bashrc
```

#### Prerequisites:
1. Python 2.6+ or 3.5+
2. Ubuntu libraries and tools:
```bash
sudo apt-get install software-properties-common -y
sudo add-apt-repository universe -y
sudo apt-get update && sudo apt-get install -y
sudo apt-get install build-essential python-dev gfortran zlib1g-dev samtools python2 wget -y
# Last two steps only required if using python2
wget https://bootstrap.pypa.io/pip/2.7/get-pip.py
sudo python2 get-pip.py
```
3. [Pysam](https://github.com/pysam-developers/pysam) verion 0.9.0  or higher is required. Flask is optional.

`sudo pip3 install pysam Cython numpy scipy matplotlib Flask future`

 **... or for python 2:**

`sudo pip2 install pysam==0.15.2 Cython numpy scipy matplotlib Flask future` 


Note that 0.15.2 is the last version of pysam which appears to support pip2 installation, however AA itself supports the more recent versions.

4. Mosek optimization tool version 8.x (https://www.mosek.com/). **Due to breaking changes in the newer versions of Mosek, we require version 8 to be used**:
```bash
wget http://download.mosek.com/stable/8.0.0.60/mosektoolslinux64x86.tar.bz2
tar xf mosektoolslinux64x86.tar.bz2
echo Please obtain license from https://www.mosek.com/products/academic-licenses/ or https://www.mosek.com/try/ and place in $PWD/mosek/8/licenses
echo export MOSEKPLATFORM=linux64x86 >> ~/.bashrc
export MOSEKPLATFORM=linux64x86
echo export PATH=$PATH:$PWD/mosek/8/tools/platform/$MOSEKPLATFORM/bin >> ~/.bashrc
echo export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$PWD/mosek/8/tools/platform/$MOSEKPLATFORM/bin >> ~/.bashrc
echo export MOSEKLM_LICENSE_FILE=$PWD/mosek/8/licenses >> ~/.bashrc
cd $PWD/mosek/8/tools/platform/linux64x86/python/3/
# OR FOR PYTHON2 
# cd $PWD/mosek/8/tools/platform/linux64x86/python/2/ # if using python2
# sudo python2 setup.py install #(--user)

sudo python3 setup.py install #(--user) [can also build locally with "pip install -e ."]

cd -
source ~/.bashrc
```
- Arial font for matplotlib (optional installation):

To get the Microsoft fonts on Ubuntu
```
sudo apt-get install fontconfig ttf-mscorefonts-installer
sudo fc-cache -f
```

#### Data repositories:
Set annotations directory and environment variable AA_DATA_REPO:
```bash
mkdir -p data_repo
echo export AA_DATA_REPO=$PWD/data_repo >> ~/.bashrc
cd $AA_DATA_REPO && touch coverage.stats && chmod a+r coverage.stats
source ~/.bashrc
```
Download and uncompress AA data repo matching the version of the reference genome used to generate the input BAM file in the $AA_DATA_REPO directory. You may have multiple annotations in the same directory, where the name of the subdirectory matches the version of the reference indicated by `--ref` argument to AA.
```
cd $AA_DATA_REPO
tar zxf $ref.tar.gz
```
The annotations may be downloaded here:
` https://datasets.genepattern.org/?prefix=data/module_support_files/AmpliconArchitect/`
Available annotations: 
* hg19
* GRCh37
* GRCh38 (hg38)
* mm10

In the data repo files, "`indexed`" indicates the BWA index is packaged as well, which is only needed if also using for alignment.

### Option 3: Other platforms:
AA can also be run through Nextflow, using the [nf-core/circdna pipeline](https://nf-co.re/circdna) constructed by [Daniel Schreyer](https://github.com/DSchreyer).

## Running AmpliconArchitect

#### AmpliconSuite-pipeline:
We provide a wrapper for jumping off at any intermediate step including generating the prerequisite BAM alignments with BWA, CNV calls for seeding and CNV seed selection. PrepareAA is available at https://github.com/jluebeck/AmpliconSuite-pipeline. We recommmend this for users who are less experienced with AA as it greatly simplifies the process of selecting CNV seed regions to feed to AA.
AmpliconSuite-pipeline can directly invoked AA if installed. 


### 1) Input data:
AA requires 2 input files:

1. Coordinate-sorted, indexed BAM file:
    * Align reads to a reference present in the `data_repo`.
    * AA has been tested on `bwa mem` on `hg19` and `GRCh38` reference genomes.
    * Recommended depth of coverage for WGS data is 5X-10X.
    * Bamfile may be downsampled using `$AA_SRC/downsample.py` or when running AA with the option `--downsample`. 
    * If sample has multiple reads groups with very different read lengths and fragment size distributions, then we recommend downsampling the bam file be selecting the read groups which have similar read lengths and fragment size distribution.
2. BED file with seed intervals:
    * One or more intervals per amplicon in the sample
    * AA has been tested on seed intervals generated as follows:
        - CNVs from CNV caller ReadDepth (with parameter file `$AA_SRC/src/read_depth_params`), Canvas and CNVkit
        - Select CNVs with copy number > 5x and size > 100kbp (default) and merge adjacent CNVs into a single interval using:

            `python2 $AA_SRC/amplified_intervals.py --bed {read_depth_folder}/output/alts.dat --out {outFileNamePrefix} --bam {BamFileName} --ref {ref}`
        - ***Note that this preprocessing step is critical to AA as it removes low-mappability and low-complexity regions. 
        - Optional argument `--ref` should match the name of the folder in `data_repo` which corresponds to the version of human reference genome used in the BAM file.

### 2) Usage:
`$AA --bam {input_bam} --bed {bed file} --out {prefix_of_output_files} <optional arguments>`

The execution script `$AA` is provided within the Github source code and the exact path depends on the installation option used:
1. Docker image: `AA=AmpliconArchitect/docker/run_aa_docker.sh`
2. Github source: `AA=python2 AmpliconArchitect/src/AmpliconArchitect.py`


#### Outputs:
AA generates informative output at each step in the algorithm (details below):
1. Summary file: List of amplicons and corresponding intervals are listed in a summary file.
2. SV view: A PNG/PDF image for each amplicon displaying all rearrangement signatures. Underlying data is provided in text format as intermediate files.
3. Graph file: For each amplicon, a text file describing the graph and predicted copy count.
4. Cycles file: For each amplicon: a text file describing the list of simple cycles predicted.
5. Cycle view: A web interface with operations for visualizing and modifying the simple cycles.

The user may provide intermediate files as a way to either kickstart AA from an intermediate step or to use alternative intermediate data (e.g. from external tools) for reconstruction.

#### Required Arguments:

| Argument | Type | Description                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         |
| ---------- | ---- |---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `--bed`         |FILE| Bed file with putative list of amplified intervals                                                                                                                                                                                                                                                                                                                                                                                                                                                                  | 
| `--bam`         |FILE| Coordinate sorted BAM file with index mapped to provided reference genome                                                                                                                                                                                                                                                                                                                                                                                                                                           | 
| `--out`         |PATH| Prefix for output files                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             | 
| `--ref`         |STR | Values: [`hg19`, `GRCh37`, `GRCh38`, `mm10`, `<CUSTOM>`, `None`]. Pick reference annotations to use from the AA_DATA_REPO directory. BAM and BED files match these annotations. <br> - `hg19`/`GRCh38` : chr1,, chr2, .. chrM etc <br> - `GRCh37` : '1', '2', .. 'MT' etc<br> - `<CUSTOM>` : User provided annotations in AA_DATA_REPO directory. <br> - `None` : do not use any annotations. AA can tolerate additional chromosomes not stated but accuracy and annotations may be affected. <br> - Default: `hg19` |

**NOTE1:** Optional argument `--ref` should match the name of the folder in `data_repo` which corresponds to the version of human reference genome used in the BAM file.

**NOTE2:** The user should be aware that AA uses intermediate files with prefix provided to `--out`. If these files are already present then AA will reuse these files. If the user intends to run AA without using any prior data, then the user should ensure that these files are not already present in the provided path (See checkpoint and modularity section).

**NOTE3:** The current docker script cannot accept paths with special characters including spaces. For paths with special characters, please install AA from the github source.

#### Optional Arguments:

| Argument | Type | Description                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                |
| ---------- | ---- |------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `-h`, `--help`    |    | Show this help message and exit                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            |
| `-v`, `--version` |    | Print program version and exit.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            |
| `--runmode`     | STR| Values: [`FULL`/`BPGRAPH`/`CYCLES`/`SVVIE`W]. This option determines which stages of AA will be run. <br> - `FULL`: Run the full reconstruction including breakpoint graph, cycles as well as SV visualization. <br> - `BPGRAPH`: Only reconstruct the breakpoint graph and estimate copy counts, but do not reconstruct the amplicon cycles. <br> - `CYCLES`: Only reconstruct the breakpoint graph and cycles, but do not create the output for SV visualization. <br> - `SVVIEW`: Only create the SV visualization, but do not reconstruct the breakpoint graph or cycles. <br> Default: `FULL`         | 
| `--extendmode`  |STR | Values: [`EXPLORE`/`CLUSTERED`/`UNCLUSTERED`/`VIRAL`]. This determines how the input intervals in bed file are treated.<br> - `EXPLORE` : Search for all intervals in the genome that may be connected to input seed intervals.<br> - `CLUSTERED` : Input intervals are treated as part of a single connected amplicon and no new connected itervals are added. <br> - `UNCLUSTERED` : Each input interval is treated as a distinct single interval amplicon and no new intervals are added.<br> Default: `EXPLORE`                                                                                        | 
| `--sensitivems` | STR| Values: [`True`, `False`]. Set `True` only if copy counts are expected to vary by an order of magnitude, e.g. viral integration. Default: `False`                                                                                                                                                                                                                                                                                                                                                                                                                                                          | 
| `--plotstyle` | STR | Values: [`small`, `large`, `all_amplicons`]. `large`: large font, `all_amplicons`: display a large number of intervals in a single plot, recommeded for visualizing multiple amplicons in CLUSTERED mode. Default: `small`                                                                                                                                                                                                                                                                                                                                                                                 |
| `--downsample`  |FLOAT| Values: [`-1`, `0`, `C`(>0)]. Decide how to downsample the bamfile during reconstruction. Reads are automatically downsampled in real time for speedup. Alternatively pre-process bam file using $AA_SRC/downsample.py. For focal amplifications with copy number < 10, setting a downsample value of 30 or 40 may yield better results. <br> - `-1` : Do not downsample bam file, use full coverage. <br> - `0` : Downsample bamfile to 10X coverage if original coverage larger then 10. <br> - `C` (>0) : Downsample bam file to coverage `C` if original coverage larger than `C`. <br> - Default: `0` | 
| `--cbam`        |FILE| Use alternative bamfile to use for coverage calculation                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    | 
| `--cbed`        |FILE| Use provided bed file for coverage calculation. Bed file defines 1000 10kbp genomic windows.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               |
| `--insert_sdevs`|FLOAT| Values: `> 0`  Number of standard deviations around the insert size. May need to increase for sequencing runs with high variance after insert size selection step. (default 3.0)                                                                                                                                                                                                                                                                                                                                                                                                                           |
| `--pair_support_min`|INT| Values: `> 1`  Number of read pairs for minimum breakpoint support (default 2 but typically becomes higher due to coverage-scaled cutoffs). (default 2)                                                                                                                                                                                                                                                                                                                                                                                                                                                    |
| `--no_cstats`   |FLAG| Values: `--no_cstats`  Do not re-use coverage statistics from `$AA_DATA_REPO/coverage.stats` file. Set this if trying multiple different values of `--insert_sdevs` or `--pair_support_min`. (default not set)                                                                                                                                                                                                                                                                                                                                                                                             |


### 3) Output description
The software generates 4 types of output files. 1 summary file and 3 files per amplicon:

| File name | Description |
| --------- | ----------- |
| `{out}_summary.txt` | This file includes a summary for all amplicons detected by AA.|
| `{out}_amplicon{id}_graph.txt` | A text file for each amplicon listing the edges in the breakpoint graph, their categorization (sequence, discordant, concordant, source) and their copy counts.|
| `{out}_amplicon{id}_cycle.txt` | A text file for each amplicon listing the simple cycles and their copy counts.|
| `{out}_amplicon{id}.png/pdf` | A PNG/PDF image file displaying the SV view of AA.|

### 4) Interpreting the output
A common question after running AA is, **"How do I know if these reconstructions represent ecDNA?"**

To aid in answering that question we have separately developed amplicon classification methods which can be run on AA output to predict the type(s) of focal amplification present. Check out [AmpliconClassifier](https://github.com/jluebeck/AmpliconClassifier).

### 5) Visualizing reconstruction:
The outputs can be visualized in 3 different formats:
- SVVIEW: This corresponds to the output file `{out}_amplicon{id}.png/pdf` and roughly shows, the discordant edges, copy number segments and genomic coverage.
- CYCLEVIEW: This is an interactive format to visualize the structure of the amplicon described by the file `{out}_amplicon{id}_cycle.txt`. The file {out}_amplicon{id}_cycle.txt and optionally {out}_amplicon{id}.png may be uploaded to `genomequery.ucsd.edu:8800` to visualize and interactively modify the cycle. Alternatively, the user may run the visualization tool locally on port 8000 using the following commands:
```bash
export FLASK_APP=$AA_SRC/cycle_visualization/web_app.py
flask run --host=0.0.0.0 --port=8000
```
This format of visualization makes it easy to discern the segments in the structure in the context of their genomic position and copy number segments as well as identify multiple copies of the segments within the same structure.

- CYCLEVIZ: We developed a python program called [CycleViz](https://github.com/jluebeck/CycleViz) to visualize elements of AA's decompositions in Circos-style plots.  

### Instructions for web interface:
- Choose and upload a cycles files generated by AA.
- Optional but recommended: Upload corresponding png file generated AA.
- Operations:
    - Show-hide cyles by name, copy count threshold
    - Merge cycles with a common segment: Select cycle names and ments rank (NOT the ID displayed) in the order of segments played. For closed cycles, first segment rank is 0; for open walks, For n walks, first segment rank is 1 and so on.
    - Pivot cycle around inverted duplication: Pivot the portion of  cycle that connects two reversed occurences of a duplicated ment without changing any genomic connections.
    - Undo last edit: Go back to previous state.

## The AA Algorithm
A full description of the methods and detailed characterization of copy number amplifications and ecDNA can be found in the manuscript referenced in the introduction.
### Definitions:
1. Amplicon: A set of genomic intervals connected together and amplified in copy number
2. Amplicon structure(s): Ordered list(s) of segments from the amplicon intervals present in the sample.
### Inputs:
1. BAM file: WGS reads mapped to the human genome
2. BED file: Set of seed intervals to be used for searching and reconstructing the amplicons in the sample. User should provide at least 1 seed per amplicon.
### Algorithm:
AA implements various steps to predict the structure of the amplicons:
1. Interval set determination: Determine the list of intervals for each amplicon to be reconstructed.
2. SV detection: Detect copy number changes and structural variations using coverage and discordant read pairs within each amplicon.
3. Breakpoint graph construction: Construct a breakpoint graph consisting of sequence egdes (genomic segments), breakpoint edges (pairs of connected genomic positions) and optionally a source vertex and predict copy counts for all edges.
4. Cycle decomposition: Decompose the breakpoint graph into simple cycles which provide a simple representation of predicted amplicon structures.
5. Interactive cycle merging: Provide a web interface to interactively merge and modify the cycles to explore candidate structures.

## File Formats
### 1. Summary file (`{out}_summary.txt`)
The summary file includes a list of all amplicons and intervals within each amplicon including size, average copy number and oncogenes amplified. The first word of each line contains the amplicon ID, the second word describes the type of information in the rest of the file. (This format will be deprecated in the future to be replaced by JSON format.)

### 2. Graph file (`{out}_amplicon{id}_graph.txt`)
The graph file for each amplicon consists of 2 sections: (i) Sequence edges and (ii) Breakpoint edges. The first line in each section starts with a header for the section. The rest of the lines start with a keyword corresponding to the edge category followed by the edge coordinates, copy numbers, coverage/number of reads and other information.

Breakpoint vertices corresponding to genomic positions are listed in the format:
`{CHROM}:{POS}{STRAND}` where:
1. `{CHROM}`: chromosome name
2. `{POS}`: 0-based coordinate where the basepair of said coordinate is included in participating segment. If the actual position is unknown or outside the amplicon interval, then `{POS}` is set to `-1` corresponding to the source vertex.
3. `{STRAND}`: `+`/`'`. `+`/`-` indicates that the participating segment ends on the right or begins on the left of the corresponding basepair respectively.

* Sequence edges section:
The edges in this section correspond to a segment on the reference genome. Sequence edges are a list of non-overlapping segments covering all the intervals in the amplicon. Tab-separated fields:
    1. Sequence Edge: keyword `sequence`
    2. StartPosition: Breakpoint vertex corresponding to the coordinate of the first base in the segment. 
    3. EndPosition: Breakpoint vertex correspondong to the coordiante of the last base in the segment.
    4. PredictedCopyCount: Predicted copy number of segment.
    5. AverageCoverage: Average depth of read coverage.
    6. Size: Number of basepairs in the segment.
    7. NumberReadsMapped: Total number of reads aligning to the segment.
* Breakpoint edges section:
The edges in the section correspond to 2 positions (breakpoint vertices) in the reference genome that are connected together in the amplicon. Tab-separated fields:
    1. BreakpointEdge: keyword `discordant` - connection between two non-consecutive positions, `concordant`: connection between consecutive positions already connected in the reference genome, `source`: connection between a known genomic position and either an unknown position or a position outside amplicon interval set.
    2. StartPosition->EndPosition: `{POSITION1}->{POSITION2}` where `{POSTION1}` and `{POSITION2}` are the coordinates and strands of the connected breakpoint vertices.
    3. PredictedCopyCount: Predicted copy number of the breakpoint edge.
    4. NumberofReadPairs: Number of discordant read pairs mapping across the breakpoint edge.
    5. HomologySizeIfAvailable(<0ForInsertions): Size of homology at the breakpoint edge in terms of number of basepairs detected using the split read alignment with the largest homology. If an insertion is detected then this field is set to negative of the size of the insertion. If split reads are not found, this column is set to `None`.
    6. Homology/InsertionSequence: Sequence of the homologous sequence or insertion at the breakpoint. If split reads are not found, this column is set to `None`. If the size is `0`, then this column is empty.

### 3. Cycles files (`{out}_amplicon{id}_cycle.txt`)
This file describes the amplicon structure predicted by AA in the form of simple cycles. This file consists of 3 sections where the first word is a keyword represting the section:
1. Intervals section: List of intervals in amplicon. Tab-separated fields:
    * `Interval`: keyword
    * `{IntervalID}`: ID of the interval 1, 2, etc
    * `{CHROM}`: Chromosome name
    * `{START}`: Coordinate of the first basepair
    * `{END}`: Coordinate of the last basepair
2. Segments section: List of segments used in the structures. Tab-separated fields:
    * `Segment`: keyword
    * `{IntervalID}`: ID of the segment 1, 2, etc. 
    * `{CHROM}`: Chromosome name
    * `{START}`: Coordinate of the first basepair
    * `{END}`: Coordinate of the last basepair
3. Cycles section: List of predicted cycles. Semi-colon separated fields:
    * `Cycle={CycleID}`: ID of the cycle, largest copy numbers first
    * `Copy_count={CopyCount}`: Copy count of the cycle
    * `Segments={Segment1}{Strand1},{Segment2}{Strand2},...`: List of segments and their orientation in the structure. The last segment loops back to connect to the first segment.


**NOTE:** Segment ID `0` is reserved for connections vertex. A cycle may contain segment ID either 0 or 2 times. A cycle containing the segment 0 indicates a linear contig where the end-points are connected either to undetermined positions or to positions outside the amplicon interval set.

### 4. The SV view ({out}_amplicon{id}.png/pdf)
The SV view file is a PNG/PDF file displaying the underlying sequence signatures of the amplicon. This image consists of:
    * The set of amplicon intervals (x-axis)
    * Window-based depth of coverage across the intervals represented as histogram (grey vertical bars)
    * Segmentation of the intervals based on coverage and copy number estimate of these segments represented by (horizontal black lines spanning the segment where the y-position of the line represents the copy number)
    * Discordant read pair clusters represented as arcs where color represents orienation of the reads. Red: Length discordant in expected orientation (forward-reverse), Brown: Everted read pairs (reverse-forward), Teal: Both reads map to forward strand and Magenta: Both reads map to the reverse strand. Vertical blue lines indicate connections to source vertex.
    * Bottom panel may represent various annotations on the amplicon intervals where the default view displays oncogene annotations.

The SV view file may be uploaded to web interface for Cycle view to visualize the cycles in conjunction with the SV view.


### 5. [Intermediate] Copy number segmentation file (`{out}_{CHROM}_{START}_{END}_cnseg.txt`)
This files provides the segmentation of an interval based on coverage alone. Here `{CHROM}_{START}_{END}` represent the coordinates of the interval. First line represents the header. Tab-separated fields:
    * `{CHROM}`: Chromosome name
    * `{START}`: Coordinate of the first basepair in the segment
    * `{END}`: Coordinate of the last basepair in the segment
    * `{CN}`: Predicted copy number of segment
    * `{start_refined}`: Whether the start position of the segment could be determined at a fine resolution (300bp)
    * `{end_refined}`: Whether the end position of the segment could be determined at a fine resolution (300bp)

If the copy number segmentation is present, AA will use this file instead of calculating the interval segmentation.

### 6. [Intermediate] Unrefined edges file (`{out}_amplicon{ampliconid}_edges.txt`)
This file provides the list of discordant edges in the amplicon. If this file is present but the refined edges file (See 7) is absent, then AA will use these edges in the amplicon reconstruction and also predict additional edges at boundaries of copy number segments (from 5) which do not have a matching discordant edge. Tab-separated fields:
    * `{POSITION1}->{POSITION2}`: Breakpoint vertices corresponding to the discordant edge.
    * NumberofReadPairs: Number of discordant read pairs mapping across the breakpoint edge.
    * HomologySizeIfAvailable(<0ForInsertions): Size of homology at the breakpoint edge in terms of number of basepairs detected using the split read alignment with the largest homology. If an insertion is detected then this field is set to negative of the size of the insertion. If split reads are not found, this column is set to `None`.
    * Homology/InsertionSequence: Sequence of the homologous sequence or insertion at the breakpoint. If split reads are not found, this column is set to `None`. If the size is `0`, then this column is empty.

### 7. [Intermediate] Refined edges file (`{out}_amplicon{ampliconid}_edges_cnseg.txt`)
This file provides the list of discordant edges in the amplicon. If this file is present, then AA will use these edges as the final set of edges in the amplicon reconstruction. Tab-separated fields:
    * `{POSITION1}->{POSITION2}`: Breakpoint vertices corresponding to the discordant edge.
    * NumberofReadPairs: Number of discordant read pairs mapping across the breakpoint edge.
    * HomologySizeIfAvailable(<0ForInsertions): Size of homology at the breakpoint edge in terms of number of basepairs detected using the split read alignment with the largest homology. If an insertion is detected then this field is set to negative of the size of the insertion. If split reads are not found, this column is set to `None`.
    * Homology/InsertionSequence: Sequence of the homologous sequence or insertion at the breakpoint. If split reads are not found, this column is set to `None`. If the size is `0`, then this column is empty.


## Checkpointing and modular integration with other tools:
The user may force AA to use existing data in order to recalculate information from previous runs or use data from external tools. 
Here are instructions for using prior or external data at various stages:
1. Interval selection: The user may select intervals from each amplicon and provide them in a separate BED file. The user may then run AA separately with each amplicon (bed file) using the option `--extendmode CLUSTERED` to indicate that AA should use all intervals within the provided as a single amplicon).
2. Copy number segmentation: The user may place the copy number segmentation file in the output directory in format described in File format section 5.
3. Unrefined discordant edges: The user may place the unrefined edges file in the output directory in format described in the File format section 6.
4. Refined discordant edges: The user may place the refined edges file in the output directory in format described in the File format section 7.
