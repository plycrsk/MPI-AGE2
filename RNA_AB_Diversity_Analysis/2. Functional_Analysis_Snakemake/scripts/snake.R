## Loading Data ##

# Configure input paths

data <- read.csv(snakemake@input[[1]])

print(data)

# Configure output paths
output_path <- snakemake@output[[1]]

write.csv(data, output_path)