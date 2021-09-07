set -xeuo pipefail
IFS=$'\n\t'

samtools faidx $reference

hap.py $truth $results -f $confidence -o $report -r $reference