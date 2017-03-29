#!/bin/bash

for rep in 4 5 6 7; do
	mkdir -p bm$rep
	mkdir -p bm$rep/sim_metadata
	mkdir -p bm$rep/fasta
	mkdir -p bm$rep/fastq
	for covV in 1 4 16 32; do
	 	for copyNV in 4 16 32; do
	 		for bpointV in 0 4 8 16 32; do
	 			for episSize in 40000 160000 640000 2400000; do
	 				python $AA_SRC/episome_simulation/ViraGE_benchmark_2_1.py $rep $covV $copyNV $bpointV $episSize
	 			done
	 		done
	 	done
	done
	mv bm$rep/sim_metadata/*fasta bm$rep/fasta
	python $AA_SRC/gen_epi_reads_2.py bm$rep/fasta
	mv bm$rep/fasta/*.fastq  bm$rep/fastq/



done

for rep in 4l 5l 6l 7l; do
	mkdir -p bm$rep
	mkdir -p bm$rep/sim_metadata
	mkdir -p bm$rep/fasta
	mkdir -p bm$rep/fastq
	for covV in 32; do
	 	for copyNV in 32; do
	 		for bpointV in 0 8 32; do
	 			for episSize in 5000000 10000000; do
	 				python $AA_SRC/ViraGE_benchmark_2_1.py $rep $covV $copyNV $bpointV $episSize
	 			done
	 		done
	 	done
	done
	mv bm$rep/sim_metadata/*fasta bm$rep/fasta
	python $AA_SRC/gen_epi_reads_2.py bm$rep/fasta
	mv bm$rep/fasta/*.fastq  bm$rep/fastq/
done




