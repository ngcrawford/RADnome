
### Prepare reads:

1. Concatenate all reads into merged r1 and r2 files.

        cat sample1.r1.fq sample2.r1.fq > merged.r1.fq
        cat sample1.r2.fq sample2.r2.fq > merged.r2.fq

### Run rainbow:

2. Cluster reads.
        
        rainbow cluster -1 merged.r1.fq > r1.cluster.out 2> r1.cluster.log
        rainbow cluster -1 merged.r2.fq > r2.cluster.out 2> r2.cluster.log

3. Divide potential groups into haplotypes.

        rainbow div -i r1.cluster.out -o r1.div.out
        rainbow div -i r2.cluster.out -o r2.div.out

4. Locally assemble merged reads into contigs.

        rainbow merge -o r1.asm.out -a -i r1.div.out
        rainbow merge -o r2.asm.out -a -i r2.div.out

### Prepare RADnome:

1. Create pseudo-genome with 50 base-pairs of N's as padding between 'contigs'.

    **Example:** [contig1] + [N * 50] + [contig2] + [N * 50] + [contig3] ... etc.

        ./make_RAD_pseudo_genome.py -n 50 -r r1.asm.RADnome r1.asm.out r1.asm.fa
        ./make_RAD_pseudo_genome.py -n 50 -r r2.asm.RADnome r2.asm.out r2.asm.fa

2. Create bowtie2 indexes.

        bowtie2-build -f /Users/ngcrawford/Dropbox/JMP_CAS/Nick/RADnome/r1.asm.fa r1.asm
        bowtie2-build -f /Users/ngcrawford/Dropbox/JMP_CAS/Nick/RADnome/r2.asm.fa r2.asm



