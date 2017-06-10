library(dada2); packageVersion("dada2")
library(ShortRead);

filtpath <- "~/data/LAB/New_Batch/filtered/" # CHANGE ME to the directory containing your filtered fastq files
fns <- list.files(filtpath, full.names = TRUE)

filts <- fns[grepl("fastq.gz$", fns)] # CHANGE if different file extensions
sample.names <- sapply(strsplit(basename(filts), ".f"), `[`, 1) # Assumes filename = samplename_XXX.fastq.gz
names(filts) <- sample.names

#------------------------------------------------
set.seed(100)
filts.learn <- sample(filts, 25) # Pick 25 samples (>100k reads/sample) to learn from
drp.learn <- derepFastq(filts.learn)
dd.learn <- dada(drp.learn, err=NULL, selfConsist=TRUE, multithread=TRUE)
err <- dd.learn[[1]]$err_out
rm(drp.learn);rm(dd.learn)
# Sample inference
dds <- vector("list", length(sample.names))
names(dds) <- sample.names
for(sam in sample.names) {
  cat("Processing:", sam, "\n")
  derep <- derepFastq(filts[[sam]])
  dds[[sam]] <- dada(derep, err=err, multithread=TRUE)
}
rm(derep)

# Construct sequence table and remove chimeras
seqtab <- makeSequenceTable(dds)
seqtab <- removeBimeraDenovo(seqtab, multithread=TRUE)
saveRDS(seqtab, "~/data/LAB/New_Batch/filtered/DADA_246.rds")
uniquesToFasta(seqtab, "~/data/LAB/New_Batch/filtered/New_DADA2_246.fasta")
seqtab <- t(seqtab)#transpose the table
seqtab <- cbind('#OTUID' = rownames(seqtab), seqtab)
write.table(seqtab, "~/data/LAB/New_Batch/filtered/dada2_OTU_table_246.txt", sep='\t', row.names=FALSE, quote=FALSE)
