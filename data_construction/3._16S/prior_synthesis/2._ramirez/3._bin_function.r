# binning functional groups in Ramirez dataset
library(dplyr)
library(purrr)
source("paths.r")
source("paths_fall2019.r")
source("NEFI_functions/crib_fun.r")

# set output path
output.path <- ramirez_tax_fun_abun.path

# read in taxonomic - functional group key
tax_fun <- readRDS(paste0(pecan_gen_16S_dir, "reference_data/bacteria_tax_to_function.rds"))

#### Read in relative abundance data ####
ramirez <- as.data.frame(readRDS(ramirez_raw_mapping_and_abundance.path))
# remove weird dataset
ramirez <- ramirez[ramirez$dataset != "X1",] # these data don't even sort of add up to 1.

tax <- ramirez[,grep("__bacteria", colnames(ramirez))]
rownames(tax) <- make.names(ramirez$sampleID)

# na_rows <- tax[rowSums(is.na(tax)) > 0,]
# na_cols <- tax[,colSums(is.na(tax)) > 0]
tax <- tax[,which(colSums(tax, na.rm=T) > 0)]
tax <- tax[which(rowSums(tax, na.rm=T) > 0),]
# add in dataset/sample identifiers
#rownames(tax) <- ramirez[which(rowSums(tax, na.rm=T) > 0),]$sampleID
tax$sampleID <- rownames(tax) 
tax$source <- "Ramirez"
# remove empty columns

tax <- tax[,c("sampleID", colnames(tax))]

# assign function to taxonomy
pathway_names <- colnames(tax_fun)[3:15]

# taxon-function assignments
tax.fun.out   <- as.data.table(tax[, c(1, grep("g__bacteria", colnames(tax)))])
genus   <- tax[, c(1, grep("g__bacteria", colnames(tax)))]
#tax.fun.out[, pathway_names] <- NA

for (i in 1:length(pathway_names)) {
  p <- pathway_names[i]
  #colname <- paste0("fg__", p)
  all.levels <- strsplit(colnames(genus), "__")
  # Classifications from literature search (multiple taxon levels)
  has_pathway <- tax_fun[tax_fun[,p] == 1,]
  levels <- c("Phylum", "Class", "Order", "Family", "Genus", "Species")
  has_pathway_taxa <- tolower(has_pathway$Taxon)
  listed <- sapply(all.levels, function(x) x %in% has_pathway_taxa)
  in_has_pathway <- sapply(listed, function(x) any(x == T))
  match <- colnames(genus)[in_has_pathway]
  # create new functional columns with the per-sample sum 
  tax.fun.out <- tax.fun.out %>% 
    dplyr::select(match) %>% 
    reduce(`+`) %>%
    mutate(tax.fun.out, newcol = .) %>% 
    dplyr::rename(!!p := newcol) 
}

rownames(tax.fun.out) <- rownames(tax)
# separate abundances and create "other" column, by functional group
fg <- tax.fun.out[, 2260:2272]
fg.list <- list()
for (i in 1:13) {
  f <- fg[,i, drop=F]
  f$other <- 1-f[,1]
  colnames(f) <- tolower(colnames(f))
  f$sampleID <- tax.fun.out$sampleID
  fg.list[[i]] <- f
}
names(fg.list) <- colnames(fg)


# subset by taxonomic level
phylum  <- tax[, c(1, grep("p__bacteria", colnames(tax)))]
class  <- tax[, c(1, grep("c__bacteria", colnames(tax)))]
order  <- tax[, c(1, grep("o__bacteria", colnames(tax)))]
family  <- tax[, c(1, grep("f__bacteria", colnames(tax)))]
genus   <- tax[, c(1, grep("g__bacteria", colnames(tax)))]

# assign taxonomic abundances by rank
tax_in <- list(phylum, class, order, family, genus)
tax_out <- list()
for (i in 1:5){
    df <- tax_in[[i]]
    df$dataset <- NULL
    df$sampleID <- NULL
    newnames <- strsplit(colnames(df), "__")
  newnames <- rapply(newnames, function(x) tail(x, 1))
  colnames(df) <- newnames
  colnames(df) <- gsub("_", " ", colnames(df))
  df <- as.data.frame(t(rowsum(t(df), group = rownames(t(df))))) # collapse all of the "unassigned" columns
  rowSums(df) > 1
  df$unassigned <- NULL
  df$other <- 0
  df[rowSums(df) < 1,]$other <- 1-rowSums(df[rowSums(df) < 1,])
  df$sampleID <- tax$sampleID
  tax_out[[i]] <- df
}
names(tax_out) <- c("phylum", "class", "order", "family", "genus")

out <- c(tax_out, fg.list)
saveRDS(out, output.path)
