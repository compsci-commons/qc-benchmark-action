set -xeuo pipefail
IFS=$'\n\t'

samtools faidx $reference

hap.py $truth $results -f $confidence_regions -o $report -r $reference --write-vcf