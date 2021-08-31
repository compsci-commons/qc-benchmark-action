set -xeuo pipefail
IFS=$'\n\t'

hap.py $truth $results -f $confident -o $report -r $reference