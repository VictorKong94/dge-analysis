# dge-analysis

### Terminal command

`Rscript <path_to_dge.R> <experiment_name> <method> <counts_directory>
<gene_annotations_file> <samples_config_file> -j <job_1> <job_2> .. <job_n>`

### Experiment name

Should be a single string unbroken by whitespace.

### Available methods

- `--et` Fisher's exact test (suitable for experiments with a single factor)
- `--lrt` likelihood ratio test
- `--qlf` quasi-likelihood F-test

### Jobs

`<job_i>` = `<baseline_sample_i>,<treatment_sample_i>`
