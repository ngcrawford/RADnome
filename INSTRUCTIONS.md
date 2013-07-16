## Install Dependencies:

### OS X:

The easiest way to install the necessary dependencies on OS X is to use the [Homebrew package manager][1].

1. Open [terminal.app][7] and from the command line run the following commands:

1. Install [Xcode][10] 

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

### Linux:

Installing the dependencies on a standard linux disro is a bit trickier than on OS X. You should be able to use `apt-get` to install [git][6]. For [tokyo cabinet][3], [samtools][4], [bowtie2][5], and [rainbow][9] you'll likely have to install from source.

Tokyo Cabinet can be a pain to install. I recommend using `deb` or `RPM` packages to install it if you can. For this reason we're working on removing it as a dependancy. 

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
[10]: https://developer.apple.com/xcode/‎

