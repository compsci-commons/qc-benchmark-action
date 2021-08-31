# QC benchmark action
A generic github action for conducting scientific QC benchmarks.


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
        uses: compsci-commons/qc-benchmark-action@main
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
        uses: compsci-commons/qc-benchmark-action@main
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

## Contributing

To extend this action with additional benchmarks,

* add a new subfolder to the `benchmarks` directory,
* and copy the contents of `benchmarks/dummy` to the new subfolder.

Then, modify the contents according to your needs.
* The `download.sh` script shall contain the code to download the necessary data for conducting the benchmark. It is invoked within the Conda environment defined in `download-env.yaml`.
* The `eval.sh` script shall contain the code to evaluate results generated on the benchmark data. It is invoked within the Conda environment defined in `eval-env.yaml`.

Finally, extend the `README.md` file in the repository root with a description of and an example for using your benchmark.

For inspiration, have a look at `benchmarks/giab-na12878-exome`.