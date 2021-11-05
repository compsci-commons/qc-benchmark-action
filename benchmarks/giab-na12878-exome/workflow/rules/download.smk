rule get_reads:
    output:
        r1=reads[0],
        r2=reads[1],
    log:
        "logs/download-reads.log",
    conda:
        "../tools.yaml"
    shell:
        "samtools view -f3 -u "
        "ftp://ftp-trace.ncbi.nih.gov/ReferenceSamples/giab/data/NA12878/Nebraska_NA12878_HG001_TruSeq_Exome/NIST-hg001-7001-ready.bam "
        "21 | samtools sort -n -u | samtools fastq -1 {output.r1} -2 {output.r2} -0 /dev/null - 2> {log}"


rule get_truth:
    output:
        "truth.vcf",
    log:
        "logs/get-truth.log",
    params:
        repl_chr=repl_chr,
    conda:
        "../tools.yaml"
    shell:
        "bcftools view "
        "https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/NA12878_HG001/NISTv4.2.1/GRCh38/HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz "
        "chr21 | sed {params.repl_chr} > {output} 2> {log}"


rule get_confidence_bed:
    output:
        "confidence-regions.bed",
    log:
        "logs/get-confidence-regions.log",
    params:
        repl_chr=repl_chr,
    conda:
        "../tools.yaml"
    shell:
        "curl --insecure -L "
        "https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/NA12878_HG001/NISTv4.2.1/GRCh38/HG001_GRCh38_1_22_v4.2.1_benchmark.bed | "
        "sed {params.repl_chr} > {output} 2> {log}"


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


rule samtools_faidx:
    input:
        "reference.fasta",
    output:
        "reference.fasta.fai",
    wrapper:
        "0.79.0/bio/samtools/faidx"


rule bwa_index:
    input:
        "reference.fasta",
    output:
        multiext("reference", ".amb", ".ann", ".bwt", ".pac", ".sa"),
    log:
        "logs/bwa-index.log",
    params:
        prefix=get_io_prefix(lambda input, output: output[0]),
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
        index=get_io_prefix(lambda input, output: input.index[0]),
        sorting="samtools",  # Can be 'none', 'samtools' or 'picard'.
        sort_order="coordinate",  # Can be 'queryname' or 'coordinate'.
    threads: 8
    wrapper:
        "0.79.0/bio/bwa/mem"


rule mark_duplicates:
    input:
        "mapped.bam"
    output:
        bam="mapped.dedup.bam",
        metrics="dedup.metrics.txt"
    log:
        "logs/picard-dedup.log"
    params:
        extra="REMOVE_DUPLICATES=true"
    resources:
        mem_mb=1024
    wrapper:
        "0.79.0/bio/picard/markduplicates"


rule samtools_index:
    input:
        "mapped.dedup.bam",
    output:
        "mapped.dedup.bam.bai",
    log:
        "logs/samtools-index.log",
    wrapper:
        "0.79.0/bio/samtools/index"


rule mosdepth:
    input:
        bam="mapped.dedup.bam",
        bai="mapped.dedup.bam.bai",
    output:
        "coverage.mosdepth.global.dist.txt",
        "coverage.quantized.bed.gz",
        summary="coverage.mosdepth.summary.txt",  # this named output is required for prefix parsing
    log:
        "logs/mosdepth.log",
    params:
        extra="--no-per-base",
        quantize="1:10:30:",
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
