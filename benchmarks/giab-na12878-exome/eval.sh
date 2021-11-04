set -xeuo pipefail
IFS=$'\n\t'

snakemake --directory $prefix --config results=`realpath $results` --cores 1 --use-conda