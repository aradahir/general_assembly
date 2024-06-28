# general_assembly
the pipeline for paired end influenza segment MP assembly

This snakemake pipeline is a tool for assembly paired-end sequences, preferably illumina sequences. The tool consists of the quality control, trimming adapters, aligning to the given reference genome, and producing the consensus sequences with depth coverage information for post analysis quality control. Feel free to change the references for adapting the pipeline with other segments or studies. Example files are available in the 'data' folder

# Installation
The pipeline can be installed using the anaconda environment with 'req.txt'.
```
  mamba create -n general_assembly -f req.txt
```
# Dependencies
1. snakemake
2. python
3. samtools
4. bwa-mem2
5. fastp
6. R::argparse
   
# Uses
1. Head to the working directory and wake up the snakemake environment
```
  conda activate general_assembly
```
2. Open the snakefile inside the directory and change the input_dir and output_dir to point out the working folders. Full-path is required.
3. Run the pipeline, c is for adjusting the cores.
```
  snakemake -c 40
```

# Results
The pipeline will create the folder 'result', which consists of 
1.  trimmed: for storing the sequences quality control and the pre-processed sequences (ie. filtering and adapter trimming).
2.  aligned: for storing .bam file after aligning with the references.
3.  Perbaseddepth: for storing the after processing files including depth coverage tables and figures.
4.  results: for storing the consensus sequences.
5.  log: for storing the processed informations.
6.  status: for creating the check-point of each step.

# Comment
References folder can chaged regard to the segment usages. New references need to pre-indexing with bwa-mem2 before using with the pipeline.
Other requirements may be updated later..
