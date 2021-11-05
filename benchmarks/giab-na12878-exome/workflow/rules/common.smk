repl_chr = "s/chr//"
reads = expand("reads.{read}.fq", read=[1, 2])
bwa_index = multiext("reference", ".amb", ".ann", ".bwt", ".pac", ".sa")

coverages = ["low", "callable"]

happy_report = multiext(
    "report",
    ".runinfo.json",
    ".vcf.gz",
    ".summary.csv",
    ".extended.csv",
    ".metrics.json.gz",
    ".roc.all.csv.gz",
    ".roc.Locations.INDEL.csv.gz",
    ".roc.Locations.INDEL.PASS.csv.gz",
    ".roc.Locations.SNP.csv.gz",
)


def get_io_prefix(getter):
    def inner(wildcards, input, output):
        return getter(input, output).split(".")[0]
    return inner


def get_cov_label(wildcards):
    if wildcards.cov == "low":
        return "1:5"
    if wildcards.cov == "callable":
        return "5:inf"
