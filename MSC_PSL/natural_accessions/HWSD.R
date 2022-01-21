library(Hmisc)

mdb_file <- HWSD.mdb

data <- mdb.get("HWSD.mdb", tables=NULL, lowernames=FALSE, allow=NULL,
        dateformat='%m/%d/%y', mdbexportArgs='-b strip')

as.data.frame(data)

dataframe <- data$HWSD_DATA
data$D_ROOTS

meta_data <- mdb.get("HWSD_META.mdb", tables=NULL, lowernames=FALSE, allow=NULL,
                    dateformat='%m/%d/%y', mdbexportArgs='-b strip')

write.csv(dataframe,"HWSD_dataframe.csv")
