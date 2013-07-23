## Install Dependencies:

RADnome expects a number of dependancies. While we've tried to keep the number to the bare minimum, writing our own read aligners and flatfile indexers is beyond the scope of this project.

### Required dependancies: 

- [python, version 2.7][10]
    - pysam module
- [git][6]
- [tokyo-cabinet][3]
- [samtools][4]
- [bowtie2][5]
- [rainbow][9]

Most linux distros as well as OS X should have python and git pre-installed. 

### OS X:

The easiest way to install the necessary dependencies on OS X is to use the [Homebrew package manager][1].

1. Before you begin install [Xcode][10] 

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

4. Install [pysam][11]. Unfortunately installing pysam can be a bit of a PITA as the newer versions tend to contain compilation bugs. So, my recommendation (as of 7/19/13) is to install pysam version 0.6 which seems to install without errors.

        wget https://pysam.googlecode.com/files/pysam-0.6.tar.gz
        tar -xzf pysam-0.6.tar.gz
        cd pysam-0.6
        sudo python setup.py install

### Linux:

Installing the dependencies on a standard linux disro is a bit trickier than on OS X. You should be able to use `apt-get` to install [git][6]. For [tokyo cabinet][3], [samtools][4], [bowtie2][5], and [rainbow][9] you'll likely have to install from source. The instruction for [pysam][11] are the same as for OS X.

Tokyo Cabinet can be a pain to install. I recommend using `deb` or `RPM` packages to install it if you can. For this reason we're working on removing it as a dependancy. 

## Install or Update RADnome:

[RADnome][8] is still in *beta* so we're suggesting that you clone the repository from github and install it using the *develop* flag. Rather than installing the software into the python libraries folder this approach creates symbolic links to the git repo. This way you can easily install new versions by *pulling* updates from the repo.

    cd /to/the/directory/where/you/want/to/put/the/RADnome/Source/Code
    git clone https://github.com/ngcrawford/RADnome
    cd RADnome
    python setup.py develop

### Run UnitTests:

UnitTests are located in `RADnome/tests.py` and can be run with `nosetests` or the following command:

    python RADnome/tests.py

### Run command on test data:

    cd tests/data/
    RADnome runPipeline \
    --fqs fq1.10k.fq fq2.10k.fq \
    --run-name test_run \
    --cores 3

Because RADnome was installed using the develop flag these changes immediately populate throughout the system.

### Updating RADnome:

To update the RADnome code to the newest version you just need to pull the changes from the repo. Since you installed RADnome using the `develop` flag these changes will populate through out your system.

    cd /to/the/directory/where/you/put/the/RADnome/Source/Code
    git pull origin master

Once RADnome becomes more stable we will also tag the repo with releases as well as putting the code on pypi.

## Post Analysis

We're still working on the post analysis code, but the basic steps are the following

1. Align reads
2. Infer genotypes
3. Generate population stats

### Aligning Reads:

Currently we're recommending that users align their reads using bowtie2 since bowtie2 provides a nice trade off between speed and senstitivity.

**Example commands:** More details can be found [here][12].

    bowtie2-build RADnome.fa RADnome
    bowtie2 \
    --very-sensitive-local \
    --rg ID:Plate.lane \
    --rg SM:sample_id \
    --rg LB:library_id \
    --rg PL:ILLUMINA \
    --rg-id Plate.lane \
    -x RADnome \
    -1 fastq1.1.fa \
    -2 fastq1.2.fa \
    > aligned.sam

A couple of important things to note:

It's very important that the read group (RG) info be specified as GATK won't run without it. To specify the RG flags propperly, you'll need to know the plate ID and line for each sample. This info can usually be obtained by looking at the fastq sequence read IDs.

For example the following read ID is from plate D0T4UACXX lane 3.

> @9L6V3M1:312:**D0T4UACXX**:**3**:1101:5166:1969/2

You'll also need to run each sample


MORE HERE.

Prep sams for GATK:

    samtools view -Sb aligned.sam > aligned.bam
    samtools sort aligned.bam aligned.sorted
    samtools index aligned.sorted.bam






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


