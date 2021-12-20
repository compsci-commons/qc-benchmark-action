set -xeuo pipefail
IFS=$'\n\t'

source=`dirname "$0"`

export SNAKEMAKE_OUTPUT_CACHE=$source/cache
mkdir -p $SNAKEMAKE_OUTPUT_CACHE

snakemake eval \
    --snakefile $source/workflow/Snakefile \
    --directory $prefix \
    --configfile $source/config/config.yaml \
    --config results=`realpath $results` \
    --cache stratify_regions \
    --cores 2 --use-conda