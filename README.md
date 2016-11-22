# dge-analysis

### Terminal command

`Rscript <path_to_dge.R> <method> <counts_directory> <annotations_file>
<samples_config_file> -j <job1> <job2> ...`

### Available methods

- `--et` Fisher's exact test (suitable for experiments with a single factor)
- `--lrt` likelihood ratio test
- `--qlf` quasi-likelihood F-test

### Jobs

Format is `-j <baseline1>,<treatment1> <baseline2>,<treatment2> ...`
