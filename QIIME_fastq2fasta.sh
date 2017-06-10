#!/bin/bash

#Since QIIME produces qual files after filtering which cannot be directly used by convert_fastaqual_fastq.py some editing is required.

#I am using a for loop to go over 77 samples, change the incorrect format and then use convert_fastaqual_fastq.py to convert Fasta and Qual Files back to FastQ files.


for i in {100..176}
do
	cd sample_SP$i ;

	cat seqs.qual | sed 's/\[/ /g' > seq1.qual ;

	cat seq1.qual | sed 's/\]/ /g' > seq2.qual ;

	rm seq1.qual ;

	cd .. ;

   	convert_fastaqual_fastq.py -f sample_SP$i/seqs.fna -q sample_SP$i/seq2.qual -o filtered_fastq_files ;

   	cd filtered_fastq_files ;

	mv seqs.fastq SP$i.fastq ;
	
	cd .. ;
   
done


