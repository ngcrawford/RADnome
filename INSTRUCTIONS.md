### Super basic instructions:

1. Make sure you have to PE RAD fastq files.

1. You'll need to first install git, rainbow, bowtie2, and samtools.

1. Your python version should be 2.7. Run `python --version` in terminal to check this.

1. Pull repo and install RADnome:

    git clone git:github.com/ngcrawford/RADnome
    cd RADnome
    python setup.py develop

1. Available commands.

        usage: RADnome runPipeline [-h] [--fqs FQS FQS] [-r RUN_NAME] [-c CORES]

        optional arguments:
          -h, --help            show this help message and exit
          --fqs FQS FQS         Forward and reverse fastq files.
          -r RUN_NAME, --run-name RUN_NAME
                                Name of run. Becomes fasta name.
          -c CORES, --cores CORES
                                Number of processor cores for bowtie2 to use.

1. Run command with included test data.

        /RADnome.py runPipeline \
        --fqs tests/data/fq1.10k.fq \
              tests/data/fq2.10k.fq \
              -c 3 \
              -r test_pipeline
