set -xeuo pipefail
IFS=$'\n\t'

source=`dirname "$0"`

export SNAKEMAKE_OUTPUT_CACHE=$source/snakemake-cache
mkdir -p $SNAKEMAKE_OUTPUT_CACHE

snakemake download \
    --snakefile $source/workflow/Snakefile \
    --configfile $source/config/config.yaml \
    --directory $prefix --cores 2 --use-conda \
    --cache stratify_regions \
    --show-failed-logs