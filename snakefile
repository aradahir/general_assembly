##MP assembly pipeline :MP_assembly:

#run this script from the /molbio/bin/MP_assembly using the MP_assembly environment
input_fq = 'data/'# directory of the data
workspace = 'result/'# directory of the output folders
min_depth = 20 #minimum cutoff for confident reads mapping

################no touch section##############
import subprocess, sys, os, glob, snakemake

pattern_file = '/{sample}_L001_R1_001.fastq.gz'
samples = glob_wildcards(input_fq + pattern_file).sample

rule all:
	input:
		expand(workspace + 'status/qc_{sample}_R1.txt', sample = samples),
		expand(workspace + 'status/aligned_{sample}.txt', sample = samples),
		expand(workspace + 'status/consensus_{sample}.txt', sample = samples),
		expand(workspace + 'status/depthtable_{sample}.txt', sample = samples),
		expand(workspace + 'status/depthplot_{sample}.txt', sample = samples)

rule qc:
	input:
		faR1 = input_fq + '{sample}_L001_R1_001.fastq.gz',
		faR2 = input_fq + '{sample}_L001_R2_001.fastq.gz'
	output:
		qcR1 = workspace + 'status/qc_{sample}_R1.txt', 
		outR1 = workspace + 'trimmed/trimmed_{sample}_R1.fastq.gz',
		outR2 = workspace + 'trimmed/trimmed_{sample}_R2.fastq.gz'
	params:
		html = workspace + 'trimmed/qc_{sample}.html', 
		json = workspace + 'trimmed/qc_{sample}.json' 
	log: workspace + 'log/qc_{sample}.log'
	message: 'Qc and trimming on {wildcards.sample} ..'
	shell: """
		fastp -i {input.faR1} -I {input.faR2}  -f 20 -T 20 --html {params.html} --json {params.json} -o {output.outR1} -O {output.outR2}  2> {log}
		touch {output.qcR1}
		"""

rule align_BWA:
	input:
		status = workspace + 'status/qc_{sample}_R1.txt'
	output:
		out = workspace + 'status/aligned_{sample}.txt',
		outbam = workspace + 'aligned/mp.sorted_{sample}.bam'
	params:
		faR1 = workspace+ 'trimmed/trimmed_{sample}_R1.fastq.gz',
		faR2 = workspace + 'trimmed/trimmed_{sample}_R2.fastq.gz'
	log: workspace + 'log/aligned_{sample}.log'
	message: 'Aligning BWA-mem on {wildcards.sample} ..'
	shell: """
		bwa-mem2 mem -t 10 reference/pHW2000-M_CA04.fasta {params.faR1} {params.faR2} | samtools sort | samtools view -F 4 -o {output.outbam} 2> {log} 
		touch {output.out}
		"""

rule concensus:
	input:
		status = workspace + 'status/aligned_{sample}.txt'
	output:
		con_out = workspace + 'results/mp_concensus_{sample}.fa',
		out = workspace + 'status/consensus_{sample}.txt'
	params:
		in_bam = workspace + 'aligned/mp.sorted_{sample}.bam',
		min_depth = min_depth
	log: workspace + 'log/concensus_{sample}.log'
	message: 'Making concensus sequences on {wildcards.sample} ..'
	shell: """
		samtools consensus -a -f fasta {params.in_bam} -d {params.min_depth} -o {output.con_out}
		touch {output.out}
	"""

rule depth_table:
	input:
		status = workspace + 'status/consensus_{sample}.txt'
	output:
		out = workspace + 'status/depthtable_{sample}.txt',
		out_depth = workspace + 'PerbasedDepth/depth_table_{sample}.txt'
	params:
		in_bam = workspace + 'aligned/mp.sorted_{sample}.bam'
	log: 'log/Post_analysis_{sample}.log'
	message: 'Create depth qc tables on {wildcards.sample} ..'
	shell: """
		samtools depth {params.in_bam} > {output.out_depth}
		touch {output.out}
	"""
# for plotting the depth coverage
rule depth_plot:
	input:
		status = workspace + 'status/depthtable_{sample}.txt'
	output:
		out = workspace + 'status/depthplot_{sample}.txt'
	params:
		in_txt = workspace + 'PerbasedDepth/depth_table_{sample}.txt',
		out_fig = workspace + 'PerbasedDepth/depth_plot_{sample}.pdf',
		name = '{sample}'
	log: 'log/Post_analysis_plot_{sample}.log'
	message: 'Create depth plot on {wildcards.sample}'
	shell:"""
		Rscript tools/depthPlot.R --input {params.in_txt} --output {params.out_fig} --name {params.name}
		touch {output.out}
	"""