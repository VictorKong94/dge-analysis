# dge-analysis

### Terminal command

`Rscript <path_to_dge.R> <counts_directory> <experiment> <method>
<annotations_file> <samples_config_file> <job1>,<job2>,...`

### Experiment name

Should be a single string unbroken by whitespace.

### Available methods

- `--et` Fisher's exact test (suitable for experiments with a single factor)
- `--lrt` likelihood ratio test
- `--qlf` quasi-likelihood F-test

### Jobs

Format is `<baseline1>:<treatment1>,<baseline2>:<treatment2>,...`
