set -xeuo pipefail
IFS=$'\n\t'

region="21:10270925-19823611"
repl_chr="s/chr//"

# download test data (only keep proper pairs to avoid reads pointing to other chromosomes)
samtools view -f 3 -u \
    ftp://ftp-trace.ncbi.nih.gov/ReferenceSamples/giab/data/NA12878/Nebraska_NA12878_HG001_TruSeq_Exome/NIST-hg001-7001-ready.bam \
    $region | \
samtools sort -n -u | \
samtools fastq -1 $read1 -2 $read2 -0 /dev/null -

# download truth
bcftools view \
    https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/NA12878_HG001/NISTv4.2.1/GRCh38/HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz \
    chr$region | sed $repl_chr > $truth

# download confidence bed
curl --insecure -L \
    https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/NA12878_HG001/NISTv4.2.1/GRCh38/HG001_GRCh38_1_22_v4.2.1_benchmark.bed \
    > $confidence_regions_all

# download liftover
curl --insecure -L \
    http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz \
    > $liftover

# download target regions
liftOver \
    <(curl --insecure -L ftp://ftp-trace.ncbi.nih.gov/ReferenceSamples/giab/data/NA12878/Nebraska_NA12878_HG001_TruSeq_Exome/TruSeq_exome_targeted_regions.hg19.bed) \
    $liftover \
    $target_regions \
    /dev/null

# intersect target and confidence regions
bedtools intersect -a $confidence_regions_all -b $target_regions | sed $repl_chr | bedtools intersect -a stdin -b <(echo $region | sed "s/[:-]/\t/g") > $confidence_regions

# download reference genome
curl --insecure -L \
    https://ftp.ensembl.org/pub/release-104/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.21.fa.gz | \
    gzip --decompress --stdout \
    > $reference
