# CRISPRESSO WORKFLOW #

## Step 1: build argument files
- Modify plate description file (*plateAnnotCrispresso.tsv*)
- Modify the input file list (*sampleName.txt*), please note that the plate number should be encoded in the file names.
- run in an R console *buildSpecificArgs.R*. Feel free to modify what you need to change! the script will ask for package installation the first time it will be launched.

## Step 2: Run crispresso
- Copy the directory of arguments files (*specificArgs*) where you will execute the workflow
- Install: snakemake, cutadapt, crispresso2, fastqc, multiqc. In theory, this can be easily done trough anaconda.org
- Modify the *config.json* 
- you can find in my RNA-Seq alignment workflow[https://github.com/DimitriMeistermann/rnaseq_align] a lot of help for this step, including a description of the arguments names in the config.json.
- Run snakemake!
- Don't forget to open the multiQC reports after the execution to control the sequencing quality

## Step 3: reanalyze crispresso outputs
- Copy the output data from the workflow
- Run the R script *analyzeResults.R*.

NB: This workflow is experimental, for example, the three steps could be reunited in one snakemake workflow.