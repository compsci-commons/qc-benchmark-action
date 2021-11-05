set -xeuo pipefail
IFS=$'\n\t'

source=`dirname "$0"`

snakemake eval --snakefile $source/workflow/Snakefile --directory $prefix --config results=`realpath $results` --cores 1 --use-conda