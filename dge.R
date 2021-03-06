started = Sys.time()
library(limma)
library(edgeR)
library(locfit)
library(ggplot2)

#########
# SETUP #
#########

# Run from command line using:
#   `Rscript path/to/dge.R <counts_directory> <experiment> <method>
#    <annotations_file> <samples_config_file> <job1>,<job2>,...`
args = commandArgs(trailingOnly = F)
args = args[grep("--file=", args):length(args)]
args[1] = substring(args[1], 8)
args = args[-grep("--args", args)]
command = paste("Rscript", paste(args, collapse = " "))

# Pull parameters passed from command line
if (length(args) != 7) {
  stop(paste("Use `Rscript path/to/dge.R <counts_directory> <experiment>",
             "<method> <annotations_file> <samples_config_file>",
             "<job1>,<job2>,...`"))
}
# - directory containing *.counts files
dir = args[2]
if (substring(dir, nchar(dir)) != "/") dir = paste0(dir, "/")
# - experiment name
experiment = args[3]
# - method (either "--et", "--lrt", or "--qlf")
method = args[4]
if (!(method %in% c("--et", "--lrt", "--qlf"))) stop("Invalid method argument")
# - annotations file
anno = args[5]
# - samples configuration file
config = read.delim(args[6], row.names = "SAMPLE")
lib_sizes = config$LIBSIZE
config$LIBSIZE = NULL
groups = config$GROUP
groups = factor(groups, levels = unique(groups))
# - job specifications
jobs = unlist(strsplit(args[7], ","))

# Specify directory to which to save analysis files
outdir = sub("counted/", paste0("analysis/", experiment, "/"), dir)
logdir = sub("counted/", "logs/analysis/", dir)
outfiles = character()

# Find and import files with counts data
files = list.files(dir, pattern = "\\.counts$")
if (length(files) == 0) {
  stop("No `*.counts` files found")
} else if (length(files) %% length(groups) != 0) {
  stop("Number of `*.counts` files does not match number of group names given")
}
dge = readDGE(files, path = dir, group = groups)
rownames(dge$counts) = sub("\\.\\d+$", "", rownames(dge$counts))

# Remove genes with levels of expressions:
# - less than 10 counts per million
# - in at least 2 samples
dge$samples$lib.size = lib_sizes
keep = rowSums(cpm(dge$counts) > 10) >= 2
dge = dge[keep, ]

# Calculate normalization factors to account for varying library sizes
dge = calcNormFactors(dge)

# Estimate common dispersion and tagwise dispersions
dge = estimateCommonDisp(dge)
dge = estimateTagwiseDisp(dge)

###################
# Create MDS Plot #
###################

pdf(NULL)
MDS = as.data.frame(plotMDS(dge)$cmdscale.out)
colnames(MDS) = c("x", "y")
MDS$group = dge$samples$group

qplot(data = MDS, x = x, y = y,
      color = group, shape = group,
      xlim = c(-3, 3), ylim = c(-3, 3),
      xlab = "Dimension 1", ylab = "Dimension 2",
      size = I(4)) +
  ggtitle(expression("MDS Plot (Leading " * Log[2] * " fold change)"))
ggsave("MDS.png", path = outdir, width = 7.5, height = 6, units = "in")
outfiles = append(outfiles, paste0(outdir, "MDS.png"))

###################################
# Compile Counts Per Million Data #
###################################

cpm = cpm(dge)
annotations = read.delim(anno, quote = "", header = T,
                         row.names = "Ensembl.Transcript.ID")
annotations = annotations[rownames(cpm),]
cpm = data.frame("Transcript.ID" = rownames(cpm),
                 "Ensembl.Gene.ID" = annotations[, "Ensembl.Gene.ID"],
                 cpm,
                 annotations[, setdiff(colnames(annotations),
                                       "Ensembl.Gene.ID")])
write.table(cpm, paste0(outdir, "cpm.txt"),
            sep = "\t", quote = F, row.names = F)
outfiles = append(outfiles, paste0(outdir, "cpm.txt"))

#####################
# DGE Analysis Test #
#####################

# Create design matrix from vector of group names
# pass

# Complete the analysis for each job specified
for (job in jobs) {
  job = strsplit(job, ":")[[1]]

  # Perform analysis using method selected by user
  if (method == "--et") {
    results = exactTest(dge, pair = c(job[1], job[2]))
  } else if (method == "--lrt") {
    # pass
  } else if (method == "--qlf") {
    # pass
  }

  fc = data.frame("logFC" = results$table$logFC,
                  row.names = row.names(results$table))
  fc$FDR = p.adjust(results$table$PValue, method = "BH")
  
  # Determine names for output files
  prefix = paste0(outdir, "DGE.", job[2], ".", job[1])
  prefix = gsub("[ -]", "_", prefix)
  
  # Produce separate files for genes that are up- vs. down-regulated
  for (direction in c("up", "down")) {
    if (direction == "up") {
      de = subset(fc, FDR < 0.05 & logFC > 1)
    } else {
      de = subset(fc, FDR < 0.05 & logFC < -1)
    }
    de[, "Transcript.ID"] = rownames(de)
    de = merge(de, cpm, by = "Transcript.ID")
    de = data.frame(Transcript = de[, "Transcript.ID"],
                    Gene = de[, "Ensembl.Gene.ID"],
                    de[, setdiff(colnames(de),
                                 c("Transcript.ID", "Ensembl.Gene.ID"))])
    de = de[order(de$FDR), ]
    write.table(de, paste(prefix, direction, "txt", sep="."),
                quote = F, sep = "\t", row.names = F)
    outfiles = append(outfiles, paste(prefix, direction, "txt", sep="."))
  }
}
finished = Sys.time()
env = paste(Sys.info()[c("sysname", "release")], collapse = " ")
version = sub("R version ", "", R.Version()$version.string)
package = sapply(names(sessionInfo()$otherPkgs),
                 function(x) as.character(packageVersion(x)))

log = paste0(
  "---\n",
  "- COMMAND: ", command, "\n",
  "- STARTED: ", started, "\n",
  "- WORKING DIR: ", getwd(), "\n",
  "- R:\n",
  "\t- VERSION: ", version, "\n",
  "\t- PACKAGES:\n",
  paste0("\t\t- ", names(package), ": ", package, collapse = "\n"), "\n",
  "- ENV: ", env, "\n",
  "- OUTFILES:\n",
  paste0("\t- ", outfiles, collapse = "\n"), "\n",
  "- FINISHED: ", finished, "\n",
  "- ELAPSED: ", finished - started, " secs\n")
write(log, paste0(logdir,
                  "dge_analysis_log_",
                  gsub(" ", "_", date()),
                  "_",
                  experiment,
                  ".yml"))
