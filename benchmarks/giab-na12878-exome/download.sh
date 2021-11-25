set -xeuo pipefail
IFS=$'\n\t'

source=`dirname "$0"`

snakemake download --cores 2 --snakefile $source/workflow/Snakefile --directory $prefix --use-conda --show-failed-logs --all-temp