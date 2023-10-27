#!/usr/bin/Rscript
#The file has been created to incorporate the DDG values for the 20 proteins in using teh database of stability_sort. 
#The following steps have been done on the virtual machine 
suppressPackageStartupMessages({
  library(data.table)
  library(duckdb)
  library(logger)
})
args <- commandArgs(trailingOnly = TRUE)
log_threshold(DEBUG) #change DEBUG to INFO if you don't want the debug lines 
hotspot.table<-fread("hotspot_table.csv")
con<-dbConnect(duckdb::duckdb(), dbdir=":memory:")
res0a<-dbExecute(con,paste0("pragma temp_directory='/mnt/tmp';"))
res0b<-dbExecute(con,"pragma memory_limit='30GB';")
res0c<-dbExecute(con,"pragma threads=8;")
duckdb::duckdb_register(con,"hotspot_table",hotspot.table)
res2<-dbExecute(con,"create or replace view asg as select uniprot, aa_pos, aa_alt, ddg FROM parquet_scan('/mnt/StabilitySort/AM_SS_G312.parquet') where fold='F1';")
log_debug('asg: {dbGetQuery(con,"select count(*) FROM asg;")}')
log_debug('master: {dbGetQuery(con,"select count(*) FROM hotspot_table;")}')
res3<-dbExecute(con,"create or replace view hotspot_joined as select distinct * FROM hotspot_table left join asg using (uniprot,aa_pos,aa_alt);")
log_debug('hotspot_joined: {dbGetQuery(con,"select count(*) FROM hotspot_joined;")}')
res4<-dbExecute(con,"copy hotspot_joined to 'hotspot_joined.csv'(format csv,header);")
dbDisconnect(con,shutdown=TRUE)