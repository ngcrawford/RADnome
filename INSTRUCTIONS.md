## Super basic instructions

Some preliminary instructions follow that should get you started.

#### Preliminaries

1. Make sure you have forward and reverse PE-RAD fastq files. The reads in each file need to have the same length. For example, in the test data all the R1 reads are 100 bp and all the R1 reads are 95 bp.

1. You'll need to install git, rainbow, bowtie2, and samtools so that they can be run on the command line.

1. Your python version should be 2.7. Run `python --version` in terminal to check this.

1. Clone repo and install RADnome:

        git clone https://github.com/ngcrawford/RADnome
        cd RADnome
        python setup.py develop

#### Test the code

1. Run command on test data.

        cd /tests/data
        RADnome runPipeline \
        --fqs tests/data/fq1.10k.fq \
              tests/data/fq2.10k.fq \
              -c 3 \
              -r test_pipeline

1. Available commands.

        usage: RADnome runPipeline [-h] [--fqs FQS FQS] [-r RUN_NAME] [-c CORES]

        optional arguments:
          -h, --help            show this help message and exit
          --fqs FQS FQS         Forward and reverse fastq files.
          -r RUN_NAME, --run-name RUN_NAME
                                Name of run. Becomes fasta name.
          -c CORES, --cores CORES
                                Number of processor cores for bowtie2 to use.