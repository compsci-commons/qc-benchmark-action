set -xeuo pipefail
IFS=$'\n\t'

snakemake --directory $prefix --cores 1 --use-conda