##MP assembly pipeline :MP_assembly:

#run this script from the /molbio/bin/MP_assembly using the MP_assembly environment "conda activate MP_assembly" may need to deactivate wfi first
input_fq = 'data/ont/'# directory of the data
workspace = 'result/ont/'# directory of the output folders
min_depth = 20 #minimum cutoff for confident reads mapping
technology = 'ont' #illumina or ont

################no touch sections##############
import glob, snakemake, os, pandas, os.path, sys

df_plasmid_gene = pandas.read_csv('config.csv', usecols = ['plasmid', 'gene_ref'], na_filter = False)
plasmid_gene_combination = df_plasmid_gene["plasmid"] + df_plasmid_gene["gene_ref"]

df_gene = pandas.read_csv('config.csv', usecols = ['gene_naming', 'gene_ref'], na_filter = False)
dictData_gene = dict(zip(df_gene.gene_naming, df_gene.gene_ref))


if technology == 'illumina':
	pattern_file = '/{sample}_L001_R1_001.fastq.gz'
	pattern_gene = '/{plasmid}-{gene}_{number}_L001_R1_001.fastq.gz'
	# example of filename : pBris60-HA_S32_L001_R1_001.fastq.gz
if technology == 'ont':
	pattern_file = '/{sample}.fastq.gz'
	pattern_gene = '/{plasmid}-{gene}_{number}_{number2}.fastq.gz'
	#ex. example of filename : pBris60-HA_SDM144_1.fastq.gz


segments = glob_wildcards(input_fq + pattern_gene).gene
renamed_segment = [dictData_gene[x] for x in segments]
plasmids = glob_wildcards(input_fq + pattern_gene).plasmid
samples = glob_wildcards(input_fq + pattern_file).sample

#auto create a new reference; giving .fasta file in reference path
for pg in plasmid_gene_combination:
	file = 'reference/'+ str(pg) + '.fasta.pac'
	if os.path.isfile(file):
		print('reference' + pg + ' is already created!')
	else: 
		new_ref = 'reference/'+ str(pg)+ '.fasta'
		if os.path.isfile(new_ref):
			os.system('bwa-mem2 index %s' % new_ref)
			print("reference " + new_ref + " is created!")
		else:
			sys.exit('There is no reference.fasta provided!')

#############################pipeline beginning##############################
rule all:
	input:
		expand(workspace + 'status/{gene}_{plasmid}_qc_{sample}_R1.txt',zip, sample = samples, gene = renamed_segment, plasmid = plasmids),
		expand(workspace + 'status/{gene}_{plasmid}_aligned_{sample}.txt',zip, sample = samples, gene = renamed_segment, plasmid = plasmids),
		expand(workspace + 'status/{gene}_{plasmid}_consensus_{sample}.txt',zip, sample = samples, gene = renamed_segment, plasmid = plasmids),
		expand(workspace + 'status/{gene}_{plasmid}_depthtable_{sample}.txt', zip,sample = samples, gene = renamed_segment, plasmid = plasmids),
		expand(workspace + 'status/{gene}_{plasmid}_depthplot_{sample}.txt', zip,sample = samples, gene = renamed_segment, plasmid = plasmids)

if technology == 'illumina':
	rule qc:
		input:
			faR1 = input_fq + '{sample}_L001_R1_001.fastq.gz',
			faR2 = input_fq + '{sample}_L001_R2_001.fastq.gz'
		output:
			qcR1 = workspace + 'status/{gene}_{plasmid}_qc_{sample}_R1.txt', 
			outR1 = workspace + 'trimmed/{gene}_{plasmid}_trimmed_{sample}_R1.fastq.gz',
			outR2 = workspace + 'trimmed/{gene}_{plasmid}_trimmed_{sample}_R2.fastq.gz'
		params:
			html = workspace + 'trimmed/{gene}_{plasmid}_qc_{sample}.html', 
			json = workspace + 'trimmed/{gene}_{plasmid}_qc_{sample}.json' 
		log: workspace + 'log/{gene}_{plasmid}_qc_{sample}.log'
		message: 'Qc and trimming on {wildcards.sample} ..'
		shell: """
			fastp -i {input.faR1} -I {input.faR2}  -f 20 -T 20 --html {params.html} --json {params.json} -o {output.outR1} -O {output.outR2}  2> {log}
			touch {output.qcR1}
			"""
	rule align_BWA:
		input:
			status = workspace + 'status/{gene}_{plasmid}_qc_{sample}_R1.txt'
		output:
			out = workspace + 'status/{gene}_{plasmid}_aligned_{sample}.txt',
			outbam = workspace + 'aligned/{gene}_{plasmid}_sorted_{sample}.bam'
		params:
			faR1 = workspace+ 'trimmed/{gene}_{plasmid}_trimmed_{sample}_R1.fastq.gz',
			faR2 = workspace + 'trimmed/{gene}_{plasmid}_trimmed_{sample}_R2.fastq.gz',
			ref_fa = 'reference/{plasmid}{gene}.fasta' 
		log: workspace + 'log/{gene}_{plasmid}_aligned_{sample}.log'
		message: 'Aligning BWA-mem on {wildcards.sample} ..'
		shell: """
			bwa-mem2 mem -t 10 {params.ref_fa} {params.faR1} {params.faR2} | samtools view -h | samtools sort | samtools view -F 4 -o {output.outbam} 2> {log} 
			touch {output.out}
			"""

if technology == 'ont':
	print(input_fq + '{sample}.fastq.gz')
	rule qc:
		input:
			fa = input_fq + '{sample}.fastq.gz'
		output:
			qc = workspace + 'status/{gene}_{plasmid}_qc_{sample}_R1.txt'
		message: 'Qc and trimming on {wildcards.sample} ..'
		shell: """
			touch {output.qc}
			"""
	rule align_minimap2:
		input:
			status = workspace + 'status/{gene}_{plasmid}_qc_{sample}_R1.txt'
		output:
			out = workspace + 'status/{gene}_{plasmid}_aligned_{sample}.txt',
			outbam = workspace + 'aligned/{gene}_{plasmid}_sorted_{sample}.bam'
		params:
			fa = input_fq + '{sample}.fastq.gz',
			ref_fa = 'reference/{plasmid}{gene}.fasta' 
		log: workspace + 'log/{gene}_{plasmid}_aligned_{sample}.log'
		message: 'Aligning minimap2 on {wildcards.sample} ..'
		shell: """
			bwa-mem2 mem -t 10 {params.ref_fa} {params.fa} | samtools view -h | samtools sort | samtools view -F 4 -o {output.outbam} 2> {log}
			touch {output.out}
			"""

rule concensus:
	input:
		status = workspace + 'status/{gene}_{plasmid}_aligned_{sample}.txt'
	output:
		con_out = workspace + 'results/{gene}_{plasmid}_concensus_{sample}.fa',
		out = workspace + 'status/{gene}_{plasmid}_consensus_{sample}.txt'
	params:
		in_bam = workspace + 'aligned/{gene}_{plasmid}_sorted_{sample}.bam',
		min_depth = min_depth
	log: workspace + 'log/{gene}_{plasmid}_concensus_{sample}.log'
	message: 'Making concensus sequences on {wildcards.sample} ..'
	shell: """
		samtools consensus -a -f fasta {params.in_bam} -d {params.min_depth} -o {output.con_out}
		touch {output.out}
	"""

rule depth_table:
	input:
		status = workspace + 'status/{gene}_{plasmid}_consensus_{sample}.txt'
	output:
		out = workspace + 'status/{gene}_{plasmid}_depthtable_{sample}.txt',
		out_depth = workspace + 'PerbasedDepth/{gene}_{plasmid}_depth_table_{sample}.txt'
	params:
		in_bam = workspace + 'aligned/{gene}_{plasmid}_sorted_{sample}.bam'
	log: 'log/{gene}_{plasmid}_Post_analysis_{sample}.log'
	message: 'Create depth qc tables on {wildcards.sample} ..'
	shell: """
		samtools depth {params.in_bam} > {output.out_depth}
		touch {output.out}
	"""
# for plotting the depth coverage
rule depth_plot:
	input:
		status = workspace + 'status/{gene}_{plasmid}_depthtable_{sample}.txt'
	output:
		out = workspace + 'status/{gene}_{plasmid}_depthplot_{sample}.txt'
	params:
		in_txt = workspace + 'PerbasedDepth/{gene}_{plasmid}_depth_table_{sample}.txt',
		out_fig = workspace + 'PerbasedDepth/{gene}_{plasmid}_depth_plot_{sample}.pdf',
		name = '{sample}'
	log: 'log/{gene}_{plasmid}_Post_analysis_plot_{sample}.log'
	message: 'Create depth plot on {wildcards.sample}'
	shell:"""
		Rscript tools/depthPlot.R --input {params.in_txt} --output {params.out_fig} --name {params.name}
		touch {output.out}
	"""
