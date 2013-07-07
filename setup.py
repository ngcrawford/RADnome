import distribute_setup
distribute_setup.use_setuptools()
from setuptools import setup

setup(
    name='RADnome',
    version='0.0.1',
    author='Nicholas G. Crawford',
    author_email='ngcrawford@gmail.com',
    # url='http://pypi.python.org/pypi/pypgen/',
    license='LICENSE.txt',
    description='Create pseudo-genomes from paired end RAD data.',
    long_description=open('README.md').read(),
    test_suite='nose.collector',
    tests_require=['nose'],
    install_requires=[
        "pysam == 0.6.0"   # Newer versions don't seem to compile correctly.
    ],

    packages=[
              'RADnome'
             ],

    package_data={
            '': ['*.txt',
                 '*.rst'
                 '*.md'],     # READMEs, etc
        },

    include_package_data=True,

    scripts=[
             'scripts/RADnome',
            ],
)