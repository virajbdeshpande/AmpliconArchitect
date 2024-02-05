# AmpliconArchitect (AA)
![GitHub release (latest by date)](https://img.shields.io/github/v/release/jluebeck/AmpliconArchitect)
![License](https://img.shields.io/badge/license-UCSD_software_license-blue)

### AmpliconArchitect is best used through [AmpliconSuite-pipeline](https://github.com/AmpliconSuite/AmpliconSuite-pipeline)
Installation instructions for AmpliconArchitect are provided here, but to prepare the inputs, invoke AA and classify the outputs, please do so by using [AmpliconSuite-pipeline](https://github.com/AmpliconSuite/AmpliconSuite-pipeline). 

### Recent updates:

### February 2024 update
- `1.3.r8` adds breakpoint microhomology detection from alignments marked as supplementary (not just secondary). Also tweaks sashimi plot visualization of CN. 

### December 2023 update
- `1.3.r7` refines the CN segmentation shown in the visualizations to prevent mismatches between displayed CN and coverage. Also allows SV VCF to use "." in the FILTER field instead of only "PASS". MOSEK convergence criteria relaxed slightly to prevent rare termination issues. 

### July 2023 update
- `1.3.r6` adds multiple new features:
  - `--sv_vcf` argument which allows users to augment AA's SV detection with their own SV calls provided in a VCF format.
  - Automated protection against improperly-formatted inputs
  - Reduces bugs created when AA is rerun into the same directory with existing files having the same sample name but different input files.
  - Bugfix for edge case where AA does not properly expand a newly discovered interval if a discovered SV lands exactly on the endpoint of the explored interval. 

### March 2023 updates:
- `1.3.r5` provides better compatibility with the AmpliconSuite-pipeline Singularity image and versions of Mosek installed via pip/conda.

- `1.3.r4` adds a bugfix to coverage plotting, some code reorganization to provide a modest speedup (approx 20% in the average case), automatic testing of the MOSEK license status, and better handling of the coverage stats lookup file.

### January 2023 update:
Version `1.3.r3` adds support for Mosek versions 9 and 10. Many thanks to the Mosek team for adding these changes (especially Michal Adamaszek). Our testing revealed that usage of different Mosek versions 
will slightly change AA copy number estimates between versions (typical difference < 0.02 copies).
`1.3.r3` makes text objects in the PDF amplicon plots editable - as a text object instead of text outline (thank you to Kaiyuan Zhu for proposing this improvement).
Now adjusting font type and size on AA output figures can be done with much more ease.
This update also adds improvements to cached coverage stats lookup and more control when using `downsample.py` manually.

**[Older update descriptions are available here.](https://docs.google.com/document/d/1jqnCs46hrpYGBGrZQFop31ezskyludxNJEQdZONwFdc/edit?usp=sharing)**

## Introduction
Focal oncogene amplification and rearrangements drive tumor growth and evolution in multiple cancer types. Proposed mechanisms for focal amplification include extrachromosomal DNA (ecDNA) formation, breakage-fusion-bridge (BFB) mechanism, tandem duplications, chromothripsis and others.
Focally amplified regions are often hotspots for genomic rearrangements. As a result, the focally amplified region may undergo rapid copy number changes and the structure of the focally amplified region may evolve over time contributing to tumor evolution. 
Furthermore, ecDNA elements may reintegrate back into the genome to form HSRs. The inter-cell heterogeneity in copy number of ecDNA as well as the interchangeability between ecDNA and HSR may allow the tumor to adapt to changing environment, e.g. targetted drug application. As a result, understanding the architecture of the focal amplifications is important to gain insights into cancer biology. AmpliconArchitect (AA) is a tool which can reconstruct the structure of focally amplified regions (>10kbp) in a cancer sample using whole genome sequence short paired-end data.

Please check out the **detailed guide** on running AA [available here](https://github.com/AmpliconSuite/AmpliconSuite-pipeline/blob/master/GUIDE.md) to learn about best practices and see some FAQs.

**AmpliconArchitect was originally developed by Viraj Deshpande**, and is maintained by Jens Luebeck, Viraj Deshpande, and others in Vineet Bafna's lab. A full description of the method can be found in the following publication:

*Deshpande, V. et al., Exploring the landscape of focal amplifications in cancer using AmpliconArchitect. Nat. Commun. 10, 392 (2019).* PMID: 30674876. [(Article)](https://www.nature.com/articles/s41467-018-08200-y)

## Table of contents:
1. [AmpliconSuite-pipeline](#recommended-way-to-run-aa-ampliconsuite-pipeline)
2. [Installation](#installation)
3. [Usage](#running-ampliconarchitect)
4. [AA outputs and arguments](#ampliconarchitect-output-files-and-command-line-arguments)
5. [The AA Algorithm](#the-aa-algorithm)
6. [Checkpointing and modular integration with other tools](#checkpointing-and-modular-integration-with-other-tools)

## Recommended way to run AA: [AmpliconSuite-pipeline](https://github.com/AmpliconSuite/AmpliconSuite-pipeline)

We provide an end-to-end wrapper, which supports entry from any intermediate step, so users may start with fastqs, or a bam file, and the wrapper enables generation of the 
 CNV calls and amplicon seed regions before running AA. After invoking AA, AmpliconSuite-pipeline calls AmpliconClassifier to enable predictions of ecDNA status, and other modes of focal amplification. AmpliconSuite-pipeline is available at
https://github.com/AmpliconSuite/AmpliconSuite-pipeline.

*Importantly, AmpliconSuite-pipeline uses all our recommended best practices*, and simplifies both upstream preparation and downstream interpretation of results. We *strongly* recommend AmpliconSuite-pipeline be used to invoke AmpliconArchitect.

**Singularity and Docker images containing AmpliconArchitect can be found on the [AmpliconSuite-pipeline GitHub page](https://github.com/AmpliconSuite/AmpliconSuite-pipeline)**

### Installation-free ways to use AA (via AmpliconSuite-pipeline):

### - GenePattern Web Interface
In collaboration with the [GenePattern](https://genepattern-notebook.org/) team, AmpliconSuite-pipeline 
can now be used from your web browser. No tool installation required. Visit https://genepattern.ucsd.edu/ to register. 
After registering and signing in, search for the "AmpliconSuite" module. 

### - Nextflow
AmpliconSuite can also be run through Nextflow, using the [nf-core/circdna pipeline](https://nf-co.re/circdna) constructed by [Daniel Schreyer](https://github.com/DSchreyer).

## Installation
AA can be installed in three ways:
1. Conda installation of AmpliconSuite-pipeline (which includes AA and all recommended modules). 
2. Obtain a containerized image of AmpliconSuite-pipeline (Docker or Singularity).
3. Manual installation from GitHub source code and manual management of dependencies.

### Option 1: Conda
**[Follow the instructions here.](https://github.com/AmpliconSuite/AmpliconSuite-pipeline#option-b-install-with-conda)**

### Option 2: Containerized images
**[Follow the instructions here.](https://github.com/AmpliconSuite/AmpliconSuite-pipeline#option-d-singularity--docker-images)**

### Option 3: Standalone installation

AmpliconSuite-pipeline (including AmpliconArchitect and AmpliconClassifier modules) can be installed manually following the instructions [here](https://github.com/AmpliconSuite/AmpliconSuite-pipeline#option-c-standalone-installation-using-the-installer-script).

AmpliconArchitect can be installed as a standalone module using the instructions [here](docs/standalone_usage.md). Note that this is not recommended, as it removes AA from the modules that prepare and filter the input bed file.
Failure to properly filter inputs can lead to extreme runtimes and false-positive calls.
    

### Setting up the AA data repo
**This is required regardless of the installation option selected above**

1. To set annotations directory and environment variable `AA_DATA_REPO`:
```bash
mkdir -p data_repo
echo export AA_DATA_REPO=$PWD/data_repo >> ~/.bashrc
cd $AA_DATA_REPO && touch coverage.stats && chmod a+r coverage.stats
source ~/.bashrc
```
2. Download and uncompress AA data repo files matching the reference genome(s) needed. Data repo files are available here: https://datasets.genepattern.org/?prefix=data/module_support_files/AmpliconArchitect.

If you have [AmpliconSuite-pipeine](https://github.com/AmpliconSuite/AmpliconSuite-pipeline) installed, you can simply do
>`AmpliconSuite-pipeline.py --download_repo [ref names]`

If not, you can download the data repo files by doing:
```bash
cd $AA_DATA_REPO
wget [url for data repo [hg19/GRCh37/GRCh38/mm10].tar.gz]
tar -xzf [hg19/GRCh37/GRCh38/mm10].tar.gz
```
Available data repo annotations: 
* hg19
* GRCh37
* GRCh38 (hg38)
* GRCh38_viral (includes oncoviral sequences)
* mm10 (GRCm38)

On the data repo download page, the suffix `indexed` indicates the BWA index is packaged as well, which is only needed if also using the packaged fasta for alignment.

## Running AmpliconArchitect
**[Please see the example commands here.](https://github.com/AmpliconSuite/AmpliconSuite-pipeline#running-ampliconsuite-pipeline)**

## AmpliconArchitect output files and command-line arguments

### Outputs
AA generates informative output at each step in the algorithm (details below):
1. Summary file: List of amplicons and corresponding intervals are listed in a summary file.
2. SV view: A PNG/PDF image for each amplicon displaying all rearrangement signatures. Underlying data is provided in text format as intermediate files.
3. Graph file: For each amplicon, a text file describing the graph and predicted copy count.
4. Cycles file: For each amplicon: a text file describing the list of simple cycles predicted.
5. Cycle view: A web interface with operations for visualizing and modifying the simple cycles.

The user may provide intermediate files as a way to either kickstart AA from an intermediate step or to use alternative intermediate data (e.g. from external tools) for reconstruction.

### Required Arguments

| Argument | Type | Description                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        |
| ---------- | ---- |----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `--bed`         |FILE| Bed file with putative list of amplified intervals                                                                                                                                                                                                                                                                                                                                                                                                                                                 | 
| `--bam`         |FILE| Coordinate sorted BAM file with index mapped to provided reference genome                                                                                                                                                                                                                                                                                                                                                                                                                          | 
| `--out`         |PATH| Prefix for output files                                                                                                                                                                                                                                                                                                                                                                                                                                                                            | 
| `--ref`         |STR | Values: [`hg19`, `GRCh37`, `GRCh38`, `mm10`, `<CUSTOM>`, `None`]. Pick reference annotations to use from the AA_DATA_REPO directory. BAM and BED files match these annotations. <br> - `hg19`/`GRCh38` : chr1,, chr2, .. chrM etc <br> - `GRCh37` : '1', '2', .. 'MT' etc<br> - `<CUSTOM>` : User provided annotations in AA_DATA_REPO directory. <br> - `None` : do not use any annotations. AA can tolerate additional chromosomes not stated but accuracy and annotations may be affected. <br> |


**NOTE1:** Optional argument `--ref` should match the name of the folder in `data_repo` which corresponds to the version of human reference genome used in the BAM file.

**NOTE2:** The user should be aware that AA uses intermediate files with prefix provided to `--out`. If these files are already present then AA will reuse these files. If the user intends to run AA without using any prior data, then the user should ensure that these files are not already present in the provided path (See checkpoint and modularity section).

**NOTE3:** The current docker script cannot accept paths with special characters including spaces. For paths with special characters, please install AA from the github source.

### Optional Arguments

| Argument | Type  | Description                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                |
| ---------- |-------|------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `-h`, `--help`  |       | Show this help message and exit                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            |
| `-v`, `--version` |       | Print program version and exit.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            |
| `--runmode`     | STR   | Values: [`FULL`/`BPGRAPH`/`CYCLES`/`SVVIE`W]. This option determines which stages of AA will be run. <br> - `FULL`: Run the full reconstruction including breakpoint graph, cycles as well as SV visualization. <br> - `BPGRAPH`: Only reconstruct the breakpoint graph and estimate copy counts, but do not reconstruct the amplicon cycles. <br> - `CYCLES`: Only reconstruct the breakpoint graph and cycles, but do not create the output for SV visualization. <br> - `SVVIEW`: Only create the SV visualization, but do not reconstruct the breakpoint graph or cycles. <br> Default: `FULL`         | 
| `--extendmode`  | STR   | Values: [`EXPLORE`/`CLUSTERED`/`UNCLUSTERED`/`VIRAL`]. This determines how the input intervals in bed file are treated.<br> - `EXPLORE` : Search for all intervals in the genome that may be connected to input seed intervals.<br> - `CLUSTERED` : Input intervals are treated as part of a single connected amplicon and no new connected itervals are added. <br> - `UNCLUSTERED` : Each input interval is treated as a distinct single interval amplicon and no new intervals are added.<br> Default: `EXPLORE`                                                                                        | 
| `--sensitivems` | STR   | Values: [`True`, `False`]. Set `True` only if copy counts are expected to vary by an order of magnitude, e.g. viral integration. Default: `False`                                                                                                                                                                                                                                                                                                                                                                                                                                                          | 
| `--plotstyle` | STR   | Values: [`small`, `large`, `all_amplicons`]. `large`: large font, `all_amplicons`: display a large number of intervals in a single plot, recommeded for visualizing multiple amplicons in CLUSTERED mode. Default: `small`                                                                                                                                                                                                                                                                                                                                                                                 |
| `--downsample`  | FLOAT | Values: [`-1`, `0`, `C`(>0)]. Decide how to downsample the bamfile during reconstruction. Reads are automatically downsampled in real time for speedup. Alternatively pre-process bam file using $AA_SRC/downsample.py. For focal amplifications with copy number < 10, setting a downsample value of 30 or 40 may yield better results. <br> - `-1` : Do not downsample bam file, use full coverage. <br> - `0` : Downsample bamfile to 10X coverage if original coverage larger then 10. <br> - `C` (>0) : Downsample bam file to coverage `C` if original coverage larger than `C`. <br> - Default: `0` | 
| `--cbam`        | FILE  | Use alternative bamfile to use for coverage calculation                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    | 
| `--cbed`        | FILE  | Use provided bed file for coverage calculation. Bed file defines 1000 10kbp genomic windows.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               |
| `--insert_sdevs`| FLOAT | Values: `> 0`  Number of standard deviations around the insert size. May need to increase for sequencing runs with high variance after insert size selection step. (default 3.0)                                                                                                                                                                                                                                                                                                                                                                                                                           |
| `--pair_support_min`| INT   | Values: `> 1`  Number of read pairs for minimum breakpoint support (default 2 but typically becomes higher due to coverage-scaled cutoffs). (default 2)                                                                                                                                                                                                                                                                                                                                                                                                                                                    |
| `--no_cstats`   | FLAG  | Values: `--no_cstats`  Do not re-use coverage statistics from `$AA_DATA_REPO/coverage.stats` file. Set this if trying multiple different values of `--insert_sdevs` or `--pair_support_min`. (default not set)                                                                                                                                                                                                                                                                                                                                                                                             |
| `--sv_vcf`      | FILE  | Provide a VCF file of externally-called SVs to augment SVs identified by AA internally. Supports SV calls from GRIDSS, DELLY, LUMPY, Manta, SvABA. Sniffles is not supported as its SV calls do not appear to follow the VCF breakend specification.                                                                                                                                                                                                                                                                                                                                                       |
| `--sv_vcf_no_filter` | FLAG | Use all external SV calls from the --sv_vcf arg, even those without "PASS" in the FILTER column.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           | 


### Output description
The software generates 4 types of output files. 1 summary file and 3 files per amplicon:

| File name | Description |
| --------- | ----------- |
| `{out}_summary.txt` | This file includes a summary for all amplicons detected by AA.|
| `{out}_amplicon{id}_graph.txt` | A text file for each amplicon listing the edges in the breakpoint graph, their categorization (sequence, discordant, concordant, source) and their copy counts.|
| `{out}_amplicon{id}_cycle.txt` | A text file for each amplicon listing the simple cycles and their copy counts.|
| `{out}_amplicon{id}.png/pdf` | A PNG/PDF image file displaying the SV view of AA.|


## File Formats
### 1. Summary file `{out}_summary.txt`
The summary file includes a list of all amplicons and intervals within each amplicon including size, average copy number and oncogenes amplified. The first word of each line contains the amplicon ID, the second word describes the type of information in the rest of the file. (This format will be deprecated in the future to be replaced by JSON format.)

### 2. Graph file `{out}_amplicon{id}_graph.txt`
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

### 3. Cycles files `{out}_amplicon{id}_cycle.txt`
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
3. Cycles section: List of predicted cycles. Semicolon separated fields:
    * `Cycle={CycleID}`: ID of the cycle, largest copy numbers first
    * `Copy_count={CopyCount}`: Copy count of the cycle
    * `Segments={Segment1}{Strand1},{Segment2}{Strand2},...`: List of segments and their orientation in the structure. The last segment loops back to connect to the first segment.


**NOTE:** Segment ID `0` is reserved for connections vertex. A cycle may contain segment ID either 0 or 2 times. A cycle containing the segment 0 indicates a linear contig where the end-points are connected either to undetermined positions or to positions outside the amplicon interval set.

#### Interpreting the AA cycles files:
A common question after running AA is, **"How do I know if these reconstructions represent ecDNA?"**

- [AmpliconClassifier](https://github.com/AmpliconSuite/AmpliconClassifier): To aid in answering that question we have separately developed amplicon classification methods which can be run on AA output to predict the type(s) of focal amplification present.

- [CycleViz](https://github.com/AmpliconSuite/CycleViz): We developed a python program called [CycleViz](https://github.com/AmpliconSuite/CycleViz) to visualize elements of AA's decompositions in Circos-style plots.  


Except in structurally simple cases, the decompositions reported by AA in the cycles file may represent computational substructures of a larger complete ecDNA of unknown structure. Because the signatures of ecDNA are still present in these substructures, for predictions on which entries are ecDNA-like, and for the predicted genome intervals captured on ecDNA please use [AmpliconClassifier](https://github.com/AmpliconSuite/AmpliconClassifier) (AC) on your AA outputs to predict ecDNA status. Annotated versions of these cycles files indicating which cycles appear ecDNA-like, bed files and additional summary tables are directly produced by AC, simplifying interpretation. Most importantly however, individual AA cycles should not be interpreted as complete ecDNA reconstructions without first ruling out that these are substructures of the amplified regions.


### 4. The SV view `{out}_amplicon{id}.png/pdf`
The SV view file is a PNG/PDF file displaying the underlying sequence signatures of the amplicon. This image consists of:
- The set of amplicon intervals (x-axis)
- Window-based depth of coverage across the intervals represented as histogram (grey vertical bars)
- Segmentation of the intervals based on coverage and copy number estimate of these segments represented by (horizontal black lines spanning the segment where the y-position of the line represents the copy number)
- Discordant read pair clusters represented as arcs where color represents orienation of the reads. Red: Length discordant in expected orientation (forward-reverse), Brown: Everted read pairs (reverse-forward), Teal: Both reads map to forward strand and Magenta: Both reads map to the reverse strand. Vertical colored lines going to the top of the plot indicate connections to source vertex.
Thickness of the arc qualitatively depicts the amount of paired-end read support.
- Bottom panel may represent various annotations on the amplicon intervals where the default view displays oncogene annotations.

The SV view file may be uploaded to web interface for Cycle view to visualize the cycles in conjunction with the SV view.


### 5. [Intermediate] Copy number segmentation file `{out}_{CHROM}_{START}_{END}_ws10000_cnseg.txt`
This files provides the segmentation of an interval based on coverage alone. The value after `ws` indicates the window size used for the segmentation. Here `{CHROM}_{START}_{END}` represent the coordinates of the interval. First line represents the header. Tab-separated fields:
- `{CHROM}`: Chromosome name
- `{START}`: Coordinate of the first basepair in the segment
- `{END}`: Coordinate of the last basepair in the segment
- `{CN}`: Predicted copy number of segment
- `{start_refined}`: Whether the start position of the segment could be determined at a fine resolution (300bp)
- `{end_refined}`: Whether the end position of the segment could be determined at a fine resolution (300bp)

If the copy number segmentation is present, AA will use this file instead of calculating the interval segmentation.

### 6. [Intermediate] Unrefined edges file `{out}_amplicon{ampliconid}_edges.txt`
This file provides the list of discordant edges in the amplicon. If this file is present but the refined edges file (See 7) is absent, then AA will use these edges in the amplicon reconstruction and also predict additional edges at boundaries of copy number segments (from 5) which do not have a matching discordant edge. Tab-separated fields:
- `{POSITION1}->{POSITION2}`: Breakpoint vertices corresponding to the discordant edge.
- NumberofReadPairs: Number of discordant read pairs mapping across the breakpoint edge.
- HomologySizeIfAvailable(<0ForInsertions): Size of homology at the breakpoint edge in terms of number of basepairs detected using the split read alignment with the largest homology. If an insertion is detected then this field is set to negative of the size of the insertion. If split reads are not found, this column is set to `None`.
- Homology/InsertionSequence: Sequence of the homologous sequence or insertion at the breakpoint. If split reads are not found, this column is set to `None`. If the size is `0`, then this column is empty.

### 7. [Intermediate] Refined edges file `{out}_amplicon{ampliconid}_edges_cnseg.txt`
This file provides the list of discordant edges in the amplicon. If this file is present, then AA will use these edges as the final set of edges in the amplicon reconstruction. Tab-separated fields:
- `{POSITION1}->{POSITION2}`: Breakpoint vertices corresponding to the discordant edge.
- NumberofReadPairs: Number of discordant read pairs mapping across the breakpoint edge.
- HomologySizeIfAvailable(<0ForInsertions): Size of homology at the breakpoint edge in terms of number of basepairs detected using the split read alignment with the largest homology. If an insertion is detected then this field is set to negative of the size of the insertion. If split reads are not found, this column is set to `None`.
- Homology/InsertionSequence: Sequence of the homologous sequence or insertion at the breakpoint. If split reads are not found, this column is set to `None`. If the size is `0`, then this column is empty.


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


## Checkpointing and modular integration with other tools
The user may force AA to use existing data in order to recalculate information from previous runs or use data from external tools. 
Here are instructions for using prior or external data at various stages:
1. Interval selection: The user may select intervals from each amplicon and provide them in a separate BED file. The user may then run AA separately with each amplicon (bed file) using the option `--extendmode CLUSTERED` to indicate that AA should use all intervals within the provided as a single amplicon).
2. Copy number segmentation: The user may place the copy number segmentation file in the output directory in format described in File format section 5.
3. Unrefined discordant edges: The user may place the unrefined edges file in the output directory in format described in the File format section 6.
4. Refined discordant edges: The user may place the refined edges file in the output directory in format described in the File format section 7.
