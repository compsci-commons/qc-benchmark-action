set -xeuo pipefail
IFS=$'\n\t'

hap.py $truth $results -f $confidence -o $report -r $reference