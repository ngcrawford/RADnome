## Install Dependencies:

### OS X:

The easiest way to install the necessary dependencies on OS X is to use the [Homebrew package manager][1].

1. Open [terminal.app][7] and from the command line run the following commands:

1. Install Homebrew:

        ruby -e "$(curl -fsSL https://raw.github.com/mxcl/homebrew/go)"

    **Optional:** If you already have Homebrew installed, now would probably be a good time to upgrade your installation.

        brew upgrade

2. Install the Homebrew [science tap][2]. 'Taps' provide additional 'formula' for non standard linux tools. In this case the science tap provides formula for installing [samtools][4], and [bowtie2][5].

        brew tap homebrew/science

3. Use Homebrew to install [git][6], [tokyo cabinet][3], [samtools][4], and [bowtie2][5].

        brew install git
        brew install tokyo-cabinet
        brew install samtools
        brew install bowtie2

4. Rainbow is not currently available on 

### Linux:

Installing the dependencies on a standard linux disro is a bit trickier than on OS X. You should be able to use `apt-get` to install [git][6]. For [tokyo cabinet][3], [samtools][4], and [bowtie2][5] you'll most likely have to install from source.

## Install RADnome:

[RADnome][8] is still in *beta* so we're suggesting that you clone the repository from github and install it using the *develop* flag. Rather than installing the software into the python libraries folder this approach creates symbolic links to the git repo. This way you can easily install new versions by *pulling* updates from the repo.

        cd /to/the/directory/where/you/want/to/put/the/RADnome/Source
        git clone https://github.com/ngcrawford/RADnome
        cd RADnome
        python setup.py develop

### Run UnitTests:

To be added.


[1]: http://mxcl.github.io/homebrew/
[2]: https://github.com/Homebrew/homebrew-science
[3]: http://fallabs.com/tokyocabinet/
[4]: http://samtools.sourceforge.net/
[5]: http://bowtie-bio.sourceforge.net/bowtie2/index.shtml
[6]: http://git-scm.com/
[7]: http://en.wikipedia.org/wiki/Terminal_(OS_X)
[8]: radnome.org

