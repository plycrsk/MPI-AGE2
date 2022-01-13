input = read.csv(snakemake@input[[1]], sep='\t')

cols = c("sample_name", "unit_name", "fq1", "fq2",
        "sra", "adapters", "strandedness")

ncols = length(cols)
nrows = (nrow(input))

ncols

row_list <- list()   
for(i in 1:nrow) {             
  row_list[[i]] <- input[i,1]
}

matrix(7,14)

df <- data.frame(matrix(ncol=ncols, nrow=nrows))
colnames(df) <- cols
df$sample_name <- input$sample_name
df$sra <- input$sample_name
df$unit_name <- "lane1"
df$fq1 <- paste(df$sample_name, '_pass_1.fastq.gz', sep="")
df$fq2 <- paste(df$sample_name, '_pass_2.fastq.gz', sep="")

write.csv(df, snakemake@output[1])
