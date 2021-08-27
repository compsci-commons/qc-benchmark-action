samtools view -u \
    ftp://ftp-trace.ncbi.nih.gov/ReferenceSamples/giab/data/NA12878/Nebraska_NA12878_HG001_TruSeq_Exome/NIST-hg001-7001-ready.bam \
    22 | \
samtools fastq -1 benchmark-data/1.fastq.gz -2 benchmark-data/2.fastq.gz -0 /dev/null -