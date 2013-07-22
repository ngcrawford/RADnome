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

## Install RADnome:

[RADnome][8] is still in *beta* so we're suggesting that you clone the repository from github and install it using the *develop* flag. Rather than installing the software into the python libraries folder this approach creates symbolic links to the git repo. This way you can easily install new versions by *pulling* updates from the repo.

        cd /to/the/directory/where/you/want/to/put/the/RADnome/Source
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

[1]: http://mxcl.github.io/homebrew/
[2]: https://github.com/Homebrew/homebrew-science
[3]: http://fallabs.com/tokyocabinet/
[4]: http://samtools.sourceforge.net/
[5]: http://bowtie-bio.sourceforge.net/bowtie2/index.shtml
[6]: http://git-scm.com/
[7]: http://en.wikipedia.org/wiki/Terminal_(OS_X)
[8]: radnome.org
[9]: http://sourceforge.net/projects/bio-rainbow/
[10]: http://www.python.org/download/releases/2.7.5/
[11]: https://code.google.com/p/pysam/


