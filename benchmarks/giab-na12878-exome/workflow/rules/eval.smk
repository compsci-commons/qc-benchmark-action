
rule stratifications:
    input:
        expand("test-regions.cov-{cov}.bed", cov=coverages),
    output:
        "stratifications.tsv",
    log:
        "logs/stratification-init.log",
    run:
        with open(output[0], "w") as out:
            for cov, f in zip(coverages, input):
                print(cov, f, sep="\t", file=out)


rule merge_test_regions:
    input:
        expand("test-regions.cov-{cov}.bed", cov=coverages),
    output:
        "test-regions.all.bed",
    conda:
        "../tools.yaml"
    shell:
        "bedtools sort -i <(cat {input}) > {output} 2> {log}"


rule benchmark_variants:
    input:
        truth="truth.vcf",
        query=config["results"],
        truth_regions="test-regions.all.bed",
        strats="stratifications.tsv",
        genome="reference.fasta",
        genome_index="reference.fasta.fai",
    output:
        happy_report,
    params:
        prefix=get_output_prefix,
    log:
        "logs/happy.log",
    wrapper:
        "0.79.0/bio/hap.py/hap.py"
