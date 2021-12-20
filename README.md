# QC benchmark action
A generic GitHub action for conducting scientific QC benchmarks.


## Benchmarks

### giab-na12878-exome

Benchmark with chromosome 21 of GIAB NA12878 exome sequencing data.

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
          benchmark_name: giab-na12878-exome
      
      - name: Run your pipeline
        env:
          DATA: ${{ steps.benchmark_download.outputs.data }}
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
          benchmark_name: giab-na12878-exome
          results_path: calls.vcf.gz

      - name: Upload report as artifact
        uses: actions/upload-artifact@v2
        with:
          name: benchmark-report
          path: ${{ steps.benchmark_eval.outputs.report }}
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

### Joining the team

Of course, we are keen on adding more benchmarks. If you like to become a maintainer of a particular benchmark, please add it via an initial PR from your fork at first. Then, we will ask you to join our team here, such that you receive write access to the repository for easier maintenance.
