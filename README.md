# CRISPRESSO2 WORKFLOW #

The goal is two provide a wrapper around *CRISPRESSO2* to analyze multiplate plates at the same time of samples crisperised with HDR, using *snakemake* to manage the execution of the workflow.
Each plate can have a different set of parameters, which means different adapters, guides and amplicons.

Works with Linux/WSL environment.

## Installation

1. Install [miniconda](https://docs.anaconda.com/free/miniconda/index.html)
2. Install mamba in the base environment of conda, for example with 
```conda install conda-forge::mamba``` in a terminal.
3. Then run the follwing where you want to install the workflow:

```bash
git clone https://github.com/Kilpinen-group/crispressoWorkflow2.git
cd crispressoWorkflow2
mamba env create -f condaEnv.yml
```

## Inputs
- The compressed FASTQS (wellName_suffix_.fastq.gz). The first character of the file name must be the the plate name designated as a single letter (A, B, C, ...)
- The plate annotation file in .tsv format (see in example folder):
    - each row is a plate
    - Adapters sequences will be passed to *cutadapt*.
    - Original column: passed to *-a* argument in crispresso, amplicon sequences. Several can be provided, *e.g.* for dealing with heterozygous samples
    - afterHDR: expected sequence after HDR (*-e* in crispresso)
    - gRNA: sgRNA sequence (passed as flexiguide to crispresso, *-fg*)
    For more details go to [crispresso2 repository](https://github.com/pinellolab/CRISPResso2).
    To avoid any error, you can simply copy the [example plate annotation file](https://github.com/Kilpinen-group/crispressoWorkflow2/blob/main/example/plateAnnotCrispresso.tsv)
    
- config file (*config.json*). You can edit the one which is in the workflow directory
    Parameters:
    - THREAD_PER_SAMPLE: number of thread per task
	- FASTQ_PATH: folder path to input fastqs (they must be all in the same folder).
	- OUTPUT_PATH: folder where the outputs will be written.
	- PLATE_METADATA_PATH: path to the plate annotation file
	- PAIR_END_FILE_PATTERN: Pattern to differentiate read 1 and read 2, for example *_R* if the for A1_R1.fastq.gz A1_R2.fastq.gz.
	- FASTQ_EXTENSION: file name after the paired end pattern. for example, for A1_R1_001.fastq.gz A1_R2_001.fastq.gz, it will be _001.fastq.gz.
	- CUTADAPT_MIN_READ_LEN:: If read length less than that after cutadapt, discard read

Note, crispresso is run with parameter *-q 25*. you can modify the snakefile for altering this behavior or ad other parameters that will be applied to all samples.

## Launch
Go in the workflow directory, then:
 
```bash
conda activate crispresso #activate environment
snakemake -rp -j 5 # j --> number of max parralel task done at the same time
```

You do not have to use the config file in the workflow directory:
```bash
snakemake -rp -j 5 --configfile pathToYourFile/config.json
```

## Output

- metadata: by default, contain a generated sample list, but it's also a good folder to store a config file or a plate annotation file!
- log: log files for each rule executed by snakemake
- FASTQ_CUTADAPT: fastqs without adapters
- fastQC: results of fastqc on each fastq
- specificArgs: cutadapt/crispresso argument used for each individual samples
- CRISPRESSO: crispresso results on each sample
- CRISPAGGR: crispresso aggregate results on each plate
- results:
    - multiqc_data: data used for the multiQC report
    - multiqc_report.html: Synthesis of fastqc results
    - resultsAggr.tsv: Synthesis of crispresso results in the format used by crispresso (one row per amplicon per sample)
    - resultsFormatted.tsv: Synthesis of crispresso results with one row per sample
    - statsAggr.tsv: Per sample quality control values from crispresso
    - propAmpliconPerPlate.pdf: Plot with proportion of sequence type per sample. Page 1: Y axis = Percentage,  Y axis = count. Ambiguous = aligned but not attributed to an amplicon due to too many mismatches. Modified amplicon: Attributed to an amplicon but with some error in the sequence.

