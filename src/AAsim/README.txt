This folder includes scripts to simulate amplicon structures and rearrangements and generate short read sequence data for these structures.


Prerequisites:
Python2
ART_ILLUMINA read simulator


To run:

Open ViraGE_benchmark_2_1.py:
	Update paths in the following lines:
		VIRAL_SEQ_FILENAME="/home/jluebeck/viral_integration/data/ViralReferenceGenomes/hpv16.fasta"
		VIRAL_NAME="hpv16ref_1"
		GENOME_BUILD_FASTA="/home/jluebeck/ref_genomes/hg19.fa"
		CHROMS="/home/jluebeck/ref_genomes/hg19_chr_names.txt"

Open gen_epi_reads_2.py:
	Update paths in the following lines:
		ARTPath = "/pedigree2/projects/cancer_viral_integration/simulation/libs/art_bin_MountRainier/"


simulate_episomes.sh




