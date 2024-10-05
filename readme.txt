Author: Joshua Topper

The purpose of this script is to analyze a large amount of pooled sequences of 
patient samples and identify the mutations that may result in color change 
of a specific species of mold.

NOTE: The following files must be in your id_mutations directory.
dgorgon_reference.fa
harrington_clincal_data.txt
hawkins_pooled_sequences.fastq
pipeline.py

NOTE: You may notice that there are other items in the id_mutations directory.
These are not necessary to run the pipeline.py script. This includes:

reset.sh - this script will reset the id_mutations directory to it's original state before
the pipeline.py script was executed.

necessary_scripts - this folder contains template scripts and fake bam files
that were used to test logic in the pipeline. This folder can be ingored.

##########
To execute pipeline.py, type "python3 pipeline.py" on the command line while in the 
id_mutations directory and hit enter.
This will create the following:


Directories: 
fastqs - this will contain the fastq files for each patient and their identified reads
based on their barcode.

bams - this will contain the sorted and indexed bam files that were algined to the 
reference file.

Files:
report.txt - A complete report of mutations found in each patient fastqs that may be 
responsible for the color of their corresponding mold. 

NOTE: You will notice several different variations of the dgorgon_reference.fa file.
These files are generated throughout the script and can be ignored. 
For example: 
dgorgon_reference.fa.amb
dgorgon_reference.fa.ann