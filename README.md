# AmpliconArchitect (AA)
Focal oncogene amplification and rearrangements drive tumor growth and evolution in multiple cancer types. Proposed mechanisms for focal amplification include extrachromosomal DNA (ecDNA) formation, breakage-fusion-bridge (BFB) mechanism, tandem duplications, chromothripsis and others. Focally amplified regions are often hotspots for genomic rearrangements. As a result, the focally amplified region may undergo rapid copy number changes and the structure of the focally amplified region may evolve over time contributing to tumor progression. ecDNA originating from distinct genomic regions may recombine to form larger ecDNA elements bringing together multiple oncogenes for simultaneous amplification. Furthermore, ecDNA elements may reintegrate back into the genome to form HSRs. The inter-cell heterogeneity in copy number of ecDNA as well as the interchangeability between ecDNA and HSR may allow the tumor to adapt to changing environment, e.g. targetted drug application. As a result, understanding the architecture of the focal amplifications is important to gain insights into tumor progression as well as response to treatment. AmpliconArchitect (AA) is a tool which can reconstruct the structure of focally amplified regions in a cancer sample using whole genome sequence short paired-end data.



# Installation:

## AA download:
`git clone https://github.com/virajbdeshpande/AmpliconArchitect.git`

**Note: In the rest of this document, we will refer to the path of the parent directory `AmpliconArchitect/src` as `$AA_SRC`**

## Dependencies:
1. Python 2.7
2. Ubuntu libraries and tools:
`sudo apt-get install build-essential python-dev gfortran python-numpy python-scipy python-matplotlib python-pip zlib1g-dev samtools`
3. Pysam verion 0.9.0 or higher and Flask (optional)(https://github.com/pysam-developers/pysam):
`sudo pip install pysam Flask`
4. Mosek optimization tool (https://www.mosek.com/):
```bash
wget http://download.mosek.com/stable/8.0.0.60/mosektoolslinux64x86.tar.bz2
tar xf mosektoolslinux64x86.tar.bz2
echo Please obtain license from https://mosek.com/resources/academic-license or https://mosek.com/resources/trial-license and place in $PWD/mosek/8/licenses
echo export MOSEKPLATFORM=linux64x86 >> ~/.bashrc
export MOSEKPLATFORM=linux64x86
echo export PATH=$PATH:$PWD/mosek/8/tools/platform/$MOSEKPLATFORM/bin >> ~/.bashrc
echo export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$PWD/mosek/8/tools/platform/$MOSEKPLATFORM/bin >> ~/.bashrc
echo export MOSEKLM_LICENSE_FILE=$PWD/mosek/8/licenses >> ~/.bashrc
cd $PWD/mosek/8/tools/platform/linux64x86/python/2/
sudo python setup.py install #(--user)
cd -
source ~/.bashrc
```

## Data repositories:
Download data repositories for hg19/GRCh37 from:

https://drive.google.com/uc?export=download&confirm=V4Wy&id=0ByYcg0axX7udUDRxcTdZZkg0X1k

```bash
tar zxf data_repo.tar.gz
echo export AA_DATA_REPO=$PWD/data_repo >> ~/.bashrc
source ~/.bashrc
```


# AmpliconArchitect reconstruction:

## 1) Input data:

1. Coordinate-sorted, indexed BAM file:
    * Align reads to a reference present in the `data_repo`.
    * AA has been tested on `bwa mem` on `hg19` reference genome.
    * Recommended depth of coverage for WGS data is 5X-10X.
    * Bamfile may be downsampled using `$AA_SRC/downsample.py` or when running AA with the option `--downsample`. 
    * If sample has multiple reads groups with very different read lengths and fragment size distributions, then we recommend downsampling the bam file be selecting the read groups which have similar read lengths and fragment size distribution.
2. BED file with seed intervals:
    * One or more intervals per amplicon in the sample
    * AA has been tested on seed intervals generated as follows:
        - CNVs from CNV caller ReadDepth with parameter file `$AA_SRC/src/read_depth_params`
        - Select CNVs with copy number > 5x and size > 100kbp and merge adjacent CNVs into a single interval using:

            `python $AA_SRC/amplified_intervals.py --bed {read_depth_folder}/output/alts.dat > {outName}.bed`

## Amplicon reconstruction
`python $AA_SRC/AmpliconArchitect.py --bam {input_bam} --bed {bed file} --out {prefix_of_output_files} --ref hg19/GRCh37/{custom}`

Here `--ref` should match the name of the folder in `data_repo` which corresponds to the version of human reference genome used in the BAM file.

Optional Arguments:

| Argument | Type | Description |
| ---------- | ---- | ----------- |
| -h, --help    |    |   show this help message and exit| 
| --bed         |FILE|   Bed file with putative list of amplified intervals| 
| --bam         |FILE|   Coordinate sorted BAM file with index mapped to provided reference genome| 
| --out         |FILE|   Prefix for output files| 
| --runmode     | STR|   Values: [FULL/BPGRAPH/CYCLES/SVVIEW]. This option determines which stages of AA will be run. FULL: Run the full reconstruction including breakpoint graph, cycles as well as SV visualization. BPGRAPH: Only reconstruct the breakpoint graph and estimate copy counts, but do not reconstruct the amplicon cycles. CYCLES: Only reconstruct the breakpoint graph and cycles, but do not create the output for SV visualization. SVVIEW: Only create the SV visualization, but do not reconstruct the breakpoint graph or cycles. Default: FULL| 
| --extendmode  |STR |   Values: [EXPLORE/CLUSTERED/UNCLUSTERED/VIRAL]. This determines how the input intervals in bed file are treated. EXPLORE : Search for all connected intervals genome that may be connected to input intervals. CLUSTERED : Input intervals are treated as part of a single connected amplicon and no new connected itervals are added. UNCLUSTERED : Input intervals are treated part of a single connected amplicon and no new connected intervals are added. Default: EXPLORE| 
| --sensitivems | STR|   Values: [True, False]. Set "True" only if expected copy counts to vary by orders of magnitude, .e.g viral integration. Default: False| 
| --ref         |STR | Values: [hg19, GRCh37, None]. "hg19"(default) : chr1,, chr2, .. chrM etc / "GRCh37" : '1', '2', .. 'MT' etc/ "None" : do not use any annotations. AA can tolerate additional chromosomes not stated but accuracy and annotations may be affected. Default: hg19| 
| --downsample  |FLOAT|  Values: [-1, 0, C(>0)]. Decide how to downsample the bamfile during reconstruction. Reads are automatically downsampled in real time for speedup. Alternatively pre-process bam file using $AA_SRC/downsample.py. -1 : Do not downsample bam file, use full coverage. 0 : Downsample bamfile to 10X coverage if original coverage larger then 10. C (>0) : Downsample bam file to coverage C if original coverage larger than C. Default: 0| 
| --cbam        |FILE| Optional bamfile to use for coverage calculation| 
| --cbed        |FILE| Optional bedfile defining 1000 10kbp genomic windows for coverage calcualtion| 


## Output description
The software generates 4 types of output files. 1 summary file and 3 files per amplicon:

| File name | Description |
| --------- | ----------- |
| {out}_summary.txt | This file includes a summary for all amplicons detected |
| {out}_amplicon{id}_graph.txt | For each amplicon, breakpoint graph describing sequence edges and akpoint edges in the graph Edge each consists of 2 vertices in the format romosome>:{position}{orientation} e.g. A sequence edge chr1:10001-->chr1:20000+ represents the 00bp genomic sequence. A breakpoint edge chr1:10001-->chr1:20000+ lies that on chr1, the forward strand from position 20000 is tinues onto the forward strand at position 10001 in the 5' to 3' ection to create a 10000bp loop. Breakpoint vertices with position are reserved for source vertex at the end of linear contig. |
| {out}_amplicon{id}_cycle.txt | For each amplicon, this file represents the list of structures tained in the amplicon. Segments represent genomic intervals (with possible overlaps) and h segment is assigned a unique id. Segment id 0 is reserved for source vertices. Cycles are represented as ordered lists of connected segments along h orientation of the segments. The last segment connects back to the first segment of the cycle. Linear contigs have 0+ and 0- as their first and last segments. May be  uploaded to web interface for visualization and operations cycles (See below). |
| {out}_amplicon{id}.png | A single image file displaying the set of intervals, underlying erage, segmentation of coverage, discordant read pairs and oncogene annotations. May be uploaded to web interface for visualization and operations on cycles (See below). |

# Visualizing reconstruction:
The file {out}_amplicon{id}_cycle.txt and optionally {out}_amplicon{id}.png may be uploaded to genomequery.ucsd.edu:8800 to visualize and interactively modify the cycle.

Alternatively, the use may run the visualization tool locally on port 8000 using the following command:
```bash
export FLASK_APP=cycle_visualization/web_app.py
flask run --host=0.0.0.0 --port=8000
```
## Instructions for web interface:
- Choose and upload a cycles files generated by AA.
- Optional but recommended: Upload corresponding png file generated AA.
- Operations:
    - Show-hide cyles by name, copy count threshold
    - Merge cycles with a common segment: Select cycle names and ments rank (NOT the ID displayed) in the order of segments played. For closed cycles, first segment rank is 0; for open walks, For n walks, first segment rank is 1 and so on.
    - Pivot cycle around inverted duplication: Pivot the portion of  cycle that connects two reversed occurences of a duplicated ment without changing any genomic connections.
    - Undo last edit: Go back to previous state.

