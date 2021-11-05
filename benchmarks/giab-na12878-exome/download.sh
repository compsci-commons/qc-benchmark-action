set -xeuo pipefail
IFS=$'\n\t'

source=`dirname "$0"`

snakemake download --snakefile $source/workflow/Snakefile --directory $prefix --cores 1 --use-conda --show-failed-logs