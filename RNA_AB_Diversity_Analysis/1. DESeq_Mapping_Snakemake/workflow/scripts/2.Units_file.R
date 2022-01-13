input = read.csv(snakemake@input[[1]], sep='\t')

cols = c("sample_name", "unit_name", "fq1", "fq2",
        "sra", "adapters", "strandedness")

srr_accession = input['sample_name']
write.table(srr_accession, snakemake@output[[1]], row.names = FALSE, quote = FALSE)

ncols = length(cols)
nrows = (nrow(input))

ncols

df <- data.frame(matrix(ncol=ncols, nrow=nrows))
colnames(df) <- cols
df$sample_name <- input$sample_name
df$sra <- input$sample_name
df$unit_name <- "lane1"
df$fq1 <- paste(df$sample_name, '_pass_1.fastq.gz', sep="")
df$fq2 <- paste(df$sample_name, '_pass_2.fastq.gz', sep="")

write.csv(df, snakemake@output[[2]], row.names = FALSE, quote = FALSE)
