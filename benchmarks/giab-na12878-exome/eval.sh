set -xeuo pipefail
IFS=$'\n\t'

hap.py $truth $results -f $confidence_regions -o $report -r $reference --write-vcf