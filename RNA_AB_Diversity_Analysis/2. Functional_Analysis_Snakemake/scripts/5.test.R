ego_summary <- read.csv(snakemake@input[[1]])
ego_summary2 <- read.csv(snakemake@input[[2]])

print("hello")
#ego_summary
#ego_summary2

dev.copy2pdf(file = snakemake@output[[1]])