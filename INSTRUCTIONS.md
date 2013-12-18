## Install Dependencies:

RADnome expects a number of dependencies. While we've tried to keep the number to the bare minimum, writing our own read aligners and flat-file indexers is beyond the scope of this project.

### Required dependencies:

- [python, version 2.7][10]
    - [pysam][11] module
- [git][6]
- [tokyo-cabinet][3]
- [samtools][4]
- [bowtie2][5]
- [rainbow][9]

Most linux distros as well as recent versions of OS X should have python and git pre-installed.

### OS X:

The easiest way to install the necessary dependencies on OS X is to use the [Homebrew package manager][1].

1. Before you begin install [Xcode][10].

3. You should also install the command line tools which you can do directly from the terminal.

        xcode-select --install

1. Open [terminal.app][7] and from the command line run the following commands:

1. Install Homebrew:

        ruby -e "$(curl -fsSL https://raw.github.com/mxcl/homebrew/go)"

    **Optional:** If you already have Homebrew installed, now would probably be a good time to upgrade your installation.

        brew upgrade

2. Install the Homebrew [science tap][2]. 'Taps' provide additional 'formula' for non standard linux tools. In this case the science tap provides formula for installing [samtools][4], and [bowtie2][5].

        brew tap homebrew/science

3. Use Homebrew to install [git][6], [tokyo cabinet][3], [samtools][4], [bowtie2][5], and [rainbow][9].

        brew install git
        brew install tokyo-cabinet
        brew install samtools
        brew install bowtie2
        brew install rainbow

    **Optional:** While you're at it I recommend installing tabix and wget.

        brew install tabix
        brew install wget


### Linux:

Installing the dependencies on a standard linux disro is a bit trickier than on OS X. You should be able to use `apt-get` to install [git][6]. For [tokyo cabinet][3], [samtools][4], [bowtie2][5], and [rainbow][9] you'll likely have to install from source. The instructions for [pysam][11] are the same as for OS X.

Tokyo Cabinet can be a pain to install. I recommend using `deb` or `RPM` packages to install it if you can. For this reason we're working on removing it as a dependency.

## Install or Update RADnome:

[RADnome][8] is still in *beta* so we're suggesting that you clone the repository from github and install it using the *develop* flag. Rather than installing the software into the python libraries folder this approach creates symbolic links to the git repo. This way you can easily install new versions by *pulling* updates from the repo.

    cd /to/the/directory/where/you/want/to/put/the/RADnome/Source/Code
    git clone https://github.com/ngcrawford/RADnome
    cd RADnome
    python setup.py develop
    
The setup.py script should automatically install two required python modules: pysam and py-tcdb.

### Run UnitTests:

UnitTests are located in `RADnome/tests.py` and can be run with `nosetests` or the following command:

    python RADnome/tests.py

### Run the pipeline command on test data:

    cd tests/data/
    RADnome runPipeline \
    --fqs fq1.10k.fq fq2.10k.fq \
    --run-name test_run \
    --cores 3

### Prepare your data:

We did our best to ensure that RADnome can work with reads straight off the sequencer. Your R1 and R2 reads should be in seperate files. The read IDs are expected to be idential except for the last character which should either be 1 or 2 depending upon . In most cases this is the default. 


### Updating RADnome:

To update the RADnome code to the newest version you just need to pull the changes from the repo. Since you installed RADnome using the `develop` flag these changes will immediately populate through out your system. There is no need to run setup.py again.

    cd /to/the/directory/where/you/put/the/RADnome/Source/Code
    git pull origin master

Once RADnome becomes more stable we will also tag the repo with versioned releases and also put the code on pypi.

## Post Analysis

We're still working on the post analysis code, but the basic steps are the following:

1. Align reads
2. Infer genotypes
3. Generate population stats

We typically use [bowtie2][5] as our aligner, the Broad Institute's genome analysis toolkit ([GATK][13]) as our genotyper, and [vcftools][14] and [PGDSpider][15] to generate population genomic statistics and input files for analysis packages.

### Aligning Reads:

Currently we're recommending that users align their reads with bowtie2 since bowtie2 provides a nice trade off between speed and sensitivity. However, you could certainly used BWA, Stampy, or any other short read aligner that produces valid SAM output.

To generate SAM alignments with bowtie2  you'll first need to index the RADnome file:

    bowtie2-build test_run.RADnome.fa test_run.RADnome

It's very important that the read group (RG) information is specified when you align your reads as as GATK won't run without it. To specify the RG flags properly, you'll need to know the plate ID and line for each sample. This info can usually be obtained by looking at the fastq sequence read IDs.

For example the following read ID is from plate D0T4UACXX lane 3.

    @9L6V3M1:312:D0T4UACXX:3:1101:5166:1969/2

A bowtie2 command that includes the minimum RG info for GATK might look like the following:

    bowtie2 \
    --very-sensitive-local \
    --rg ID:plate.lane \
    --rg SM:sample_id \
    --rg LB:library_id \
    --rg PL:ILLUMINA \
    --rg-id plate.lane \
    -x test_run.RADnome \
    -1 fq1.10k.fq \
    -2 fq2.10k.fq \
    | samtools view -Sb - > test_run.RADnome.bam

You'd replace "plate", "lane", "sample_id", and "library_id" with the appropriate values. Typically sample_id and library_id are identical unless you created multiple libraries from the same sample.

[TO DO: INSTRUCTIONS FOR RUNNING MANY INDIVIDUAL ALIGNMENTS.]

### Preparing alignments for GATK:

Finally you'll need to sort and index your bam files. 

    samtools sort aligned.bam aligned.sorted
    samtools index aligned.sorted.bam

[TO DO: INSTRUCTIONS FOR RUNNING GATK]

[TO DO: INSTRUCTIONS FOR RUNNING VCFTOOLS]

[TO DO: INSTRUCTIONS FOR RUNNING PGDSPIDER]

[1]: http://mxcl.github.io/homebrew/
[2]: https://github.com/Homebrew/homebrew-science
[3]: http://fallabs.com/tokyocabinet/
[4]: http://samtools.sourceforge.net/
[5]: http://bowtie-bio.sourceforge.net/bowtie2/index.shtml
[6]: http://git-scm.com/
[7]: http://en.wikipedia.org/wiki/Terminal_(OS_X)
[8]: radnome.org
[9]: http://sourceforge.net/projects/bio-rainbow/
[10]: https://developer.apple.com/xcode/‎
[11]: https://code.google.com/p/pysam/
[12]: http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#command-line
[13]: http://www.broadinstitute.org/gatk/
[14]: http://vcftools.sourceforge.net
[15]: http://www.cmpg.unibe.ch/software/Pgdspider/‎

