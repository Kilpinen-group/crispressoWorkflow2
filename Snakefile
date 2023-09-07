import glob, os, sys, json, shutil

WORKING_DIR = os.path.dirname(workflow.snakefile)
if config=={}:
	print("Default config file loaded, from " + WORKING_DIR + "/config.json")
	configfile: WORKING_DIR+"/config.json"

## creation of the logs subdirectory
if not os.path.exists(WORKING_DIR+"/log"):
	os.mkdir(WORKING_DIR+"/log")

#put all config variable as variable in the snakefile
for configVar in config:
	if isinstance(config[configVar], str): exec(configVar+"= '"+config[configVar]+"'")
	else: exec(configVar+"="+str(config[configVar]))


## test of the path provided in the config.json file
if not os.path.exists(FASTQ_PATH):
	print("The directory " + FASTQ_PATH + " doesn't exist. Check the field FASTQ_PATH into the config.json file.")
	sys.exit(0)
else:
	## If the path ends by /, the / is suppressed
	if ( FASTQ_PATH[-1:] == "/" ):
		FASTQ_PATH = FASTQ_PATH[:-1]

INPUT_FASTQS = glob.glob(FASTQ_PATH+'/*.fastq.gz')
if len(INPUT_FASTQS) == 0:
	print("No fastq files found in " + FASTQ_PATH + ". Check the field FASTQ_PATH into the config.json file.")
	sys.exit(0)

SAMPLES = [os.path.basename(f)[:-len(FASTQ_EXTENSION)] for f in INPUT_FASTQS]

if(OUTPUT_PATH[-1] == "/") : OUTPUT_PATH = OUTPUT_PATH[:-1]

## suppress the .R1. and .R2. elements for paired-end fastq files
SAMPLES = [itemR2 for itemR2 in SAMPLES if (PAIR_END_FILE_PATTERN+"2") not in itemR2]	
SAMPLES = [itemR1.replace((PAIR_END_FILE_PATTERN+"1"),'') for itemR1 in SAMPLES]
PAIR_SUFFIX = [PAIR_END_FILE_PATTERN+"1",PAIR_END_FILE_PATTERN+"2"]

METADATA_PATH = OUTPUT_PATH+"/metadata"
if not os.path.exists(METADATA_PATH):
	os.mkdir(METADATA_PATH)

PLATES=[]
with open(PLATE_METADATA_PATH) as f:
	# skip header
	line=f.readline()
	for line in f:
		PLATES.append(line.split("\t")[0])

SAMPLES_PER_PLATE = {}
for plate in PLATES:
	SAMPLES_PER_PLATE[plate] = [sample for sample in SAMPLES if sample[0] == plate]

#############
### RULES ###
#############

rule all:
	input:
		OUTPUT_PATH + "/results/propAmpliconPerPlate.pdf",
		OUTPUT_PATH+"/results/multiqc_report.html"

rule WRITE_SAMPLE_LIST:
	output: METADATA_PATH+"/sampleNames.txt"
	run:
		with open(output[0], 'w') as f:
			for sample in SAMPLES:
				f.write(sample+"\n")

rule BUILD_SPECIFIC_ARG:
	input:
		PLATE_METADATA_PATH,
		METADATA_PATH+"/sampleNames.txt"
	output:
		expand(OUTPUT_PATH+"/specificArgs/{sample}_{tool}.txt",sample=SAMPLES,tool=["cutadapt","crispresso"])
	log:
		out=OUTPUT_PATH+"/log/BUILD_SPECIFIC_ARG.out",
		err=OUTPUT_PATH+"/log/BUILD_SPECIFIC_ARG.err"
	shell: """
	Rscript {WORKING_DIR}/buildSpecificArgs.R {OUTPUT_PATH} {PLATE_METADATA_PATH} 1> {log.out} 2> {log.err} 
	"""

rule CUTADAPT:
	input:
		R1=FASTQ_PATH+"/{sample}"+PAIR_SUFFIX[0]+FASTQ_EXTENSION,
		R2=FASTQ_PATH+"/{sample}"+PAIR_SUFFIX[1]+FASTQ_EXTENSION,
		specificArgs=OUTPUT_PATH+"/specificArgs/{sample}_cutadapt.txt"
	output:
		R1=OUTPUT_PATH+"/FASTQ_CUTADAPT/{sample}"+PAIR_SUFFIX[0]+".fastq.gz",
		R2=OUTPUT_PATH+"/FASTQ_CUTADAPT/{sample}"+PAIR_SUFFIX[1]+".fastq.gz"
	params:
		minLen = CUTADAPT_MIN_READ_LEN,
		cpu = THREAD_PER_SAMPLE
	log:
		out=OUTPUT_PATH+"/log/CUTADAPT_{sample}.out",
		err=OUTPUT_PATH+"/log/CUTADAPT_{sample}.err"
	shell: """
	cutadaptArg=$(cat {input.specificArgs} | tr -d '\\r\\n')
	cutadapt $cutadaptArg -j {params.cpu} --minimum-length {params.minLen} -o {output.R1} -p {output.R2} {input.R1} {input.R2} 1> {log.out} 2> {log.err} 
	"""


rule FASTQC:
	input: OUTPUT_PATH+"/FASTQ_CUTADAPT/{sample}{pair}.fastq.gz"
	output: multiext(OUTPUT_PATH+"/fastQC/{sample}{pair}_fastqc",".zip",".html")
	params:
		outpath=OUTPUT_PATH+"/fastQC",
		cpu = 1
	log:
		out=OUTPUT_PATH+"/log/FASTQC_{sample}{pair}.out",
		err=OUTPUT_PATH+"/log/FASTQC_{sample}{pair}.err"
	shell: """
	fastqc -o {params.outpath} {input} 1> {log.out} 2> {log.err}
	"""

for plate in PLATES:
	rule:
		input:
			R1=OUTPUT_PATH+"/FASTQ_CUTADAPT/{sample}"+PAIR_SUFFIX[0]+".fastq.gz",
			R2=OUTPUT_PATH+"/FASTQ_CUTADAPT/{sample}"+PAIR_SUFFIX[1]+".fastq.gz",
			specificArgs=OUTPUT_PATH+"/specificArgs/{sample}_crispresso.txt"
		output: OUTPUT_PATH+"/CRISPRESSO/"+plate+"/CRISPResso_on_{sample}.html"
		params:
			cpu = 1,
			plate = plate
		log:
			out=OUTPUT_PATH+"/log/CRISPRESSO_{sample}.out",
			err=OUTPUT_PATH+"/log/CRISPRESSO_{sample}.err"
		shell: """
		crispressoArg=$(cat {input.specificArgs} | tr -d '\\r\\n')
		CRISPResso --fastq_r1 {input.R1} --fastq_r2 {input.R2} -q 25 $crispressoArg -o {OUTPUT_PATH}/CRISPRESSO/{params.plate} -n {wildcards.sample} 1> {log.out} 2> {log.err}
		"""

	rule:
		input: expand(OUTPUT_PATH+"/CRISPRESSO/"+plate+"/CRISPResso_on_{sample}.html",sample=SAMPLES_PER_PLATE[plate])
		output:
			OUTPUT_PATH+"/CRISPAGGR/CRISPRessoAggregate_on_"+plate+"/CRISPRessoAggregate_mapping_statistics.txt",
			OUTPUT_PATH+"/CRISPAGGR/CRISPRessoAggregate_on_"+plate+"/CRISPRessoAggregate_quantification_of_editing_frequency_by_amplicon.txt"
		params:
			cpu = 1,
			plate = plate
		log:
			out=OUTPUT_PATH+"/log/CRISPRESSO_AGGR.out",
			err=OUTPUT_PATH+"/log/CRISPRESSO_AGGR.err"
		shell: """
		cd {OUTPUT_PATH}/CRISPAGGR
		CRISPRessoAggregate --place_report_in_output_folder --prefix {OUTPUT_PATH}/CRISPRESSO/{params.plate} -n {params.plate}  1> {log.out} 2> {log.err}
		"""

rule MULTIQC:
	input:
		fastqc=expand(OUTPUT_PATH+"/fastQC/{sample}{pair}_fastqc{ext}", sample=SAMPLES,pair=PAIR_SUFFIX,ext=[".zip",".html"])
	output: OUTPUT_PATH+"/results/multiqc_report.html"
	params:
		outpath = OUTPUT_PATH + "/results",
		cpu = 1
	shell: """
	multiqc -f -e general_stats -e tophat -e bowtie2 {OUTPUT_PATH} -o {params.outpath}
	"""

rule ANALYZE_RESULTS:
	input: 
		expand(OUTPUT_PATH+"/CRISPAGGR/CRISPRessoAggregate_on_{plate}/CRISPRessoAggregate_mapping_statistics.txt",plate=PLATES)
	output:
		OUTPUT_PATH + "/results/propAmpliconPerPlate.pdf"
	log:
		out=OUTPUT_PATH+"/log/ANALYZE_RESULTS.out",
		err=OUTPUT_PATH+"/log/ANALYZE_RESULTS.err"
	shell: """
	Rscript {WORKING_DIR}/analyzeResults.R {OUTPUT_PATH} {PLATE_METADATA_PATH} 1> {log.out} 2> {log.err} 
	"""