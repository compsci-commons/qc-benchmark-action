chromosome = "15"
repl_chr = "s/chr//"
reads = expand("reads.{read}.fq", read=[1, 2])
bwa_index = multiext("reference", ".amb", ".ann", ".bwt", ".pac", ".sa")

coverages = ["low", "medium", "high"]


def get_io_prefix(getter):
    def inner(wildcards, input, output):
        return getter(input, output).split(".")[0]

    return inner


def get_cov_label(wildcards):
    if wildcards.cov == "low":
        return "1:10"
    if wildcards.cov == "medium":
        return "10:30"
    if wildcards.cov == "high":
        return "30:inf"
