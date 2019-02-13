#' common_group_quantification.r
#' Supply a sv/otu table, a taxonomy file, and a phylogenetic level.
#' Return sample abundances, relative abundances and sequence depth for all taxa at that level 
#' found in greater than X% of samples (default 50%).
#' Also returns the number of unique SVs/OTUs that went into each group.
#'
#' @param sv           #sv or otu table where columns are unique taxonomic assignments, rows are samples.
#' @param tax          #taxonomy for sv/otu table. 7 columns for k/p/c/o/f/g/s.
#' @param groups       #unique groups you are looking for. ex: c('Ascomycota','Basidiomycota')
#' @param phyla_level  #level of phylogeny in quotes. i.e. 'phylum'.
#'
#' @return
#' @export
#'
#' @examples
common_group_quantification <- function(sv, tax, groups, tax_level, samp_freq = 0.5, ref_filter = F){
  #some tests.
  if(ncol(sv) != nrow(tax)){
    stop('Number of columns of sv table does not match number of rows in taxonomy table.')
  }

  #merge taxonomy and sv file.
  k <- cbind(tax, t(sv))
  k <- data.table(k)
  to_change <- colnames(k)[colnames(k) %in% tax_level]
  setnames(k,to_change,'tax_level')
  #count sequence abundance in each group by sample, and number of samples in which it occurs.
  abundance <- list()
  frequency <- list()
  unique.sv <- list()
  diversity <- list()
   evenness <- list()
  for(i in 1:length(groups)){
    z <- k[tax_level == groups[i],]
    start <- ncol(tax) + 1
    seq.count <- colSums(z[,start:ncol(z)])
    seq.freq  <- length(seq.count[seq.count > 0]) / nrow(sv)
    #get diversity, richness and evenness.
    seq.div <- rowSums(z[,start:ncol(z)])
    seq.div <- seq.div[seq.div > 0]
    n.SVs   <- length(seq.div)
    div     <- vegan::diversity(seq.div / sum(seq.div)) #same whether proportional or not.
    n.SVs <- length(seq.div[seq.div > 0])
    even <- div/log(n.SVs)
    
    abundance[[i]] <- seq.count
    frequency[[i]] <- seq.freq
    unique.sv[[i]] <- n.SVs
    diversity[[i]] <- div
     evenness[[i]] <- even
  }
  abundance <- do.call(rbind,abundance)
  frequency <- do.call(rbind,frequency)
  unique.sv <- data.frame(unlist(unique.sv))
  diversity <- data.frame(unlist(diversity))
   evenness <- data.frame(unlist( evenness))
  
  
  #name some stuff.
     unique.sv$groups <- groups
     diversity$groups <- groups
      evenness$groups <- groups
  rownames(abundance) <- groups
  colnames(abundance) <- rownames(sv)
  colnames(unique.sv) <- c('N.SVs','groups')
  colnames(diversity) <- c('diversity','groups')
  colnames( evenness) <- c( 'evenness','groups')
  
  #put together frequency table and order.
  frequency <- data.frame(cbind(groups,frequency))
  colnames(frequency)[2] <- 'sample_frequency'
  frequency$sample_frequency <- as.character(frequency$sample_frequency)
  frequency$sample_frequency <- as.numeric(frequency$sample_frequency)
  #subset to those that are found in > sam_freq of samples (default 50%).
  if(ref_filter == F){
    ref.frequency <- frequency[frequency$sample_frequency > samp_freq,]
  }
  #if ref_filter is on, don't filter by frequency, filter by supplied unique groups.
  if(ref_filter == T){
    ref.frequency <- frequency[frequency$groups %in% groups,]
  }
  #kill unknown and anything unknown
  ref.frequency <- ref.frequency[!(ref.frequency$groups %in% c('unknown','Unknown')),]
  #merge in number of OTUs, diversity and evenness in each group to frequency table.
  frequency <- merge(frequency,unique.sv, all.x = T)
  frequency <- merge(frequency,diversity, all.x = T)
  frequency <- merge(frequency, evenness, all.x = T)
  frequency <- frequency[order(-frequency$sample_frequency),]
  
  #get sequence depth per sample.
  seq_total <- rowSums(sv)
  
  #only keep groups that are found in >samp_freq% of samples.
  abundance <- abundance[rownames(abundance) %in% ref.frequency$groups,]
  abundance <- t(abundance)
  
  #get abundances and relative abundances.
          other <- seq_total - rowSums(abundance)
      abundance <- cbind(other, abundance)
  rel.abundance <- abundance / rowSums(abundance)
  
  #return output: abundances, relative abundances, sequence depth, group sample frequencies.
         to_return  <- list(abundance,rel.abundance,seq_total,frequency)
   names(to_return) <- c('abundances','rel.abundances','seq_total','group_frequencies')
  return(to_return)
}