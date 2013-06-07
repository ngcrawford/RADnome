
### Prepare Reads:

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