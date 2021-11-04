repl_chr = "s/chr//"
reads = expand("reads.{read}.fq", read=[1, 2])


def get_cov_label(wildcards):
    if wildcards.cov == "low":
        return "1:5"
    if wildcards.cov == "callable":
        return "5:inf"


rule get_reads:
    output:
        reads,
    log:
        "logs/download-reads.log",
    conda:
        "../tools.yaml"
    shell:
        "samtools view -f3 -u "
        "ftp://ftp-trace.ncbi.nih.gov/ReferenceSamples/giab/data/NA12878/Nebraska_NA12878_HG001_TruSeq_Exome/NIST-hg001-7001-ready.bam "
        "21 | samtools sort -n -u | samtools fastq -1 {output[0]} -2 {output[1]} -0 /dev/null - 2> {log}"


rule get_truth:
    output:
        "truth.vcf",
    log:
        "logs/get-truth.log",
    conda:
        "../tools.yaml"
    shell:
        "bcftools view "
        "https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/NA12878_HG001/NISTv4.2.1/GRCh38/HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz "
        "chr21 | sed {repl_chr} > {output}"


rule get_confidence_bed:
    output:
        "confidence-regions.bed",
    log:
        "logs/get-confidence-regions.log",
    conda:
        "../tools.yaml"
    shell:
        "curl --insecure -L "
        "https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/NA12878_HG001/NISTv4.2.1/GRCh38/HG001_GRCh38_1_22_v4.2.1_benchmark.bed "
        "> {output}"


rule get_chromosome:
    output:
        "reference.fasta",
    params:
        species="homo_sapiens",
        datatype="dna",
        build="GRCh38",
        release="104",
        chromosome="21",
    log:
        "logs/get-genome.log",
    wrapper:
        "0.79.0/bio/reference/ensembl-sequence"


rule bwa_index:
    input:
        "reference.fasta",
    output:
        multiext("reference", ".amb", ".ann", ".bwt", ".pac", ".sa"),
    log:
        "logs/bwa-index.log",
    params:
        prefix="reference",
    wrapper:
        "0.79.0/bio/bwa/index"


rule bwa_mem:
    input:
        reads=reads,
        index=multiext("reference", ".amb", ".ann", ".bwt", ".pac", ".sa"),
    output:
        "mapped.bam",
    log:
        "logs/bwa-mem.log",
    params:
        index="reference",
        sorting="samtools",  # Can be 'none', 'samtools' or 'picard'.
        sort_order="coordinate",  # Can be 'queryname' or 'coordinate'.
    threads: 8
    wrapper:
        "0.79.0/bio/bwa/mem"


rule samtools_index:
    input:
        "mapped.bam",
    output:
        "mapped.bam.bai",
    log:
        "logs/samtools-index.log",
    wrapper:
        "0.79.0/bio/samtools/index"


rule mosdepth:
    input:
        bam="mapped.bam",
        bai="mapped.bam.bai",
    output:
        "coverage.mosdepth.global.dist.txt",
        "coverage.quantized.bed.gz",
        summary="coverage.mosdepth.summary.txt",  # this named output is required for prefix parsing
    log:
        "logs/mosdepth.log",
    params:
        extra="--no-per-base",
        quantize="1:5:",
    wrapper:
        "0.77.0/bio/mosdepth"


rule stratify_regions:
    input:
        confidence="confidence-regions.bed",
        coverage="coverage.quantized.bed.gz",
    output:
        "test-regions.cov-{cov}.bed",
    log:
        "logs/stratify-regions/{cov}.log",
    params:
        cov_label=get_cov_label,
    conda:
        "../tools.yaml"
    shell:
        "bedtools intersect "
        "-a {input.confidence} "
        "-b <(zcat {input.coverage} | grep '{params.cov_label}') "
        "> {output} 2> {log}"
