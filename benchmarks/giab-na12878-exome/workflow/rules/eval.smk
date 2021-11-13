
rule stratifications:
    input:
        expand("benchmark/test-regions.cov-{cov}.bed", cov=coverages),
    output:
        "benchmark/stratifications.tsv",
    log:
        "logs/stratification-init.log",
    run:
        with open(output[0], "w") as out:
            for cov, f in zip(coverages, input):
                print(cov, f, sep="\t", file=out)


rule merge_test_regions:
    input:
        expand("benchmark/test-regions.cov-{cov}.bed", cov=coverages),
    output:
        "benchmark/test-regions.all.bed",
    log:
        "logs/merge-test-regions.log",
    conda:
        "../tools.yaml"
    shell:
        "bedtools sort -i <(cat {input}) > {output} 2> {log}"


rule benchmark_variants:
    input:
        truth="benchmark/truth.vcf",
        query=config["results"],
        truth_regions="benchmark/test-regions.all.bed",
        strats="benchmark/stratifications.tsv",
        genome="reference/reference.fasta",
        genome_index="reference/reference.fasta.fai",
    output:
        happy_report,
    params:
        prefix=get_io_prefix(lambda input, output: output[0]),
        engine="vcfeval",
    log:
        "logs/happy.log",
    wrapper:
        "0.79.0/bio/hap.py/hap.py"


# rule extract_false_positives:
#     input:
#         "report.vcf.gz"
#     output:
