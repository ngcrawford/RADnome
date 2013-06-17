### Prepare work enviroment:

1. Create directory structure:

        mkdir Project_Folder
        cd Project_Folder
        mkdir {RADnome,bams,bowtie2,fastqs,rainbow}

    If everything worked correctly, running `tree` should produce the following directory diagram.

        Project_Folder
        ├── RADnome
        ├── bams
        ├── bowtie2
        ├── fastqs
        └── rainbow

2. Move your reads (e.g., fastq files) into the fastq directory.

### Prepare reads:

Here we have named our fastq file with a sample name and whether they contain the first or second read (e.g., r1 or r2). We use the the file ending fq to indicate that these are fastq files. 

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

5. At this point we recommend you gzip your fastq files to save space. 

        gzip *.fq

### Prepare RADnome:

1. Create pseudo-genome with 50 base-pairs of N's as padding between 'contigs'.

    *Shorthand Example:* [contig1] + [N * 50] + [contig2] + [N * 50] + [contig3] ... etc.

        ./make_RAD_pseudo_genome.py -n 50 -r r1.asm.RADnome r1.asm.out r1.asm.fa
        ./make_RAD_pseudo_genome.py -n 50 -r r2.asm.RADnome r2.asm.out r2.asm.fa

2. Create bowtie2 indexes.

        bowtie2-build -f RADnome/r1.asm.fa r1.asm
        bowtie2-build -f RADnome/r2.asm.fa r2.asm

3. Run bowtie2.
 
    *Note: You may need to update bowtie2 as older versions have a bug that prevents the decompression of unpaired gzipped files.*

        bowtie2 --very-sensitive --end-to-end -x bowtie2/r1.asm -U fastqs/r1.fq.gz > bams/r1.sam
        bowtie2 --very-sensitive --end-to-end -x bowtie2/r2.asm -U fastqs/r2.fq.gz > bams/r2.sam

4. Convert sam to bam, and then sort and index bam

    Make bams.

        samtools view -bS bams/r1.sam > bams/r1.bam
        samtools view -bS bams/r2.sam > bams/r2.bam

    Sort bams.

        samtools sort bams/r1.bam bams/r1.sorted
        samtools sort bams/r2.bam bams/r2.sorted

    Index bams. (This step is not strictly necessary, but it's not a bad idea to have a index around.)

        samtools index bams/r1.sorted.bam 
        samtools index bams/r2.sorted.bam 

    Make sorted sams.

        samtools view -h bams/r1.sorted.bam > bams/r1.sorted.sam
        samtools view -h bams/r2.sorted.bam > bams/r2.sorted.sam

5. Associate contigs using read pairs.

    Index RADnome files.

        samtools faidx RADnome/r1.asm.fa.gz
        samtools faidx RADnome/r2.asm.fa.gz






