# qc-benchmark-action
A generic github action for conducting benchmarks


## Benchmarks

### giab-na12878-exome

Benchmark with chr21 of GIAB NA12878 exome sequencing data.

#### Example workflow

```yaml
name: Benchmark

on:
  schedule:
    - cron: 0 13 * * 1 # every monday at 1PM UTC

jobs:
  testing:
    name: Testing
    runs-on: ubuntu-latest
    steps:
      - name: Checkout Code
        uses: actions/checkout@v2
      - name: Download benchmark
        id: benchmark_download
        uses: compsci-common/qc-benchmark-action@master
        with:
          task: download
          benchmark-name: giab-na12878-exome
      
      - name: Run your pipeline
        env:
          DATA: ${{ steps.benchmark_eval.outputs.data }}
        run: |
          # run your pipeline with the test data downloaded above
          # Read location: $DATA/reads.1.fq and $DATA/reads.2.fq
          # Reference genome: $DATA/reference.fa
          snakemake --cores 1 ... 

      - name: Eval benchmark
        uses: ./
        id: benchmark_eval
        with:
          task: eval
          benchmark-name: giab-na12878-exome
          results-path: calls.vcf

      - name: Show results
        env:
          REPORT: ${{ steps.benchmark_eval.outputs.report }}
        run: |
          cat $REPORT.summary.csv
```
