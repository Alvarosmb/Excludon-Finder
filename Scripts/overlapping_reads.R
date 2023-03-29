############################################################################
### CONVERT BAM TABLE TO SINGLE_READ COUNT OBJECT INCLUDING METADATA INFORMATION
###########################################################################
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# LOAD LIBRARIES AND PLOTTING FUNCTION
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
defaultW <- getOption("warn") 
options(warn = -1) 
#setwd("../../media/alvaro/Elements/Operon_finder/Automated_pipe/")
source("Scripts/load_libraries.R")

input_fasta <- commandArgs(trailingOnly = TRUE)[1]
input_gff <- commandArgs(trailingOnly = TRUE)[2]
input_bam <- commandArgs(trailingOnly = TRUE)[3]
sample <-commandArgs(trailingOnly = TRUE)[4]
mode <- commandArgs(trailingOnly = TRUE)[5]
print(mode)
split_ratio <- commandArgs(trailingOnly = TRUE)[6]
split_ratio <- as.integer(split_ratio)

#Get list of RNA types
gff <- rtracklayer::import(input_gff)
gff_df <- as.data.frame(gff)
list_types <- unique(gff_df$type)
list_types <- list_types[list_types != "gene"]
interesting_list <- list_types[list_types != "rRNA"]

# LOAD FUNCTIONS
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#...................................unlist bam to datatable
.unlist <- function (x){
  x1 <- x[[1L]]
  if (is.factor(x1)){
    structure(unlist(x), class = "factor", levels = levels(x1))
  } else {
    do.call(c, x)
  }
}

#...................................wrapper function for bam input (count with featurecounts and allow (!) multi mapping reads)
counts_wrapper <- function (input_bam_file, input_fasta_file, input_gff_file){
  
  # > read in fasta file
  fasta <- readDNAStringSet(input_fasta_file)
  
  # > use featurecounts to calculate counts of mapped reads to features (CDS, tRNA, rRNA) 
  datalist  <-  list()
  datalist2 <- list()
  datafile <- stringr::str_split(input_bam_file, "/")
  
  
  for (i in seq_along(interesting_list)){
    name <- interesting_list[i]
    dir.create(paste("featurecounts_data_",name, sep = ""),showWarnings = TRUE)
    datalist[[i]] <- featureCounts(allowMultiOverlap = T, files = input_bam_file, annot.ext = input_gff_file,
                                   isGTFAnnotationFile = T, GTF.featureType = name, GTF.attrType = "gene_id",
                                   isLongRead = T,nthreads = 8, reportReads = "CORE")
    datalist2[[i]] <-fread(paste(datafile[[1]][length(datafile[[1]])], ".featureCounts", sep = "")) %>%
      dplyr::rename(id = V1, gene = V4) %>%
      dplyr::select(id, gene) %>%
      dplyr::filter(!is.na(gene)) %>%
      mutate(mapped_type = name)
  }
  return(list(datalist, datalist2))
}

#...................................combine bam, summary, gff and featurecounts file
wrapper_bam_to_table <- function (input_bam_file, input_gff_file, input_fasta_file, datalist_input, output = c("read_ids", "gene_ids")){
  
  # > read in fasta file
  fasta <- readDNAStringSet(input_fasta_file)
  
  # > read in gff file and grep for feature numbers
  
  gff_temp <- rtracklayer::import(input_gff)
  keep <- gff_temp$type %in% interesting_list
  gff_table <- data.frame(id_name = gff_temp$gene_id[keep],
                          locus_name = gff_temp$gene_id[keep],
                          start_gene = start(gff_temp)[keep],
                          end_gene = end(gff_temp)[keep],
                          strand_gene = droplevels(as.factor(strand(gff_temp)[keep])))
  
  if(output == "gene_ids"){
    feature_list <- c(rep(interesting_list[1],length(datalist_input[[1]]$annotation$Chr)))
    
    big_data_all <- cbind(datalist_input[[1]]$annotation, datalist_input[[1]]$counts) %>%
      as_tibble() %>%
      mutate(type = feature_list) %>%
      dplyr::rename(counts = 7) %>%
      rowwise() %>%
      dplyr::filter(Strand == "+" | Strand == "-") 
    
    listofdfs <- list()
    # > enable calculation for different chromosomes
    for(i in 1:length(names(fasta))){
      names(fasta) <-  stringr::str_split_fixed(names(fasta), " ", 2)[,1]
      used_chr <- stringr::str_split_fixed(names(fasta[i]), " ", 2)[,1]
      
      
      df <- big_data_all %>%
        dplyr::filter(Chr == used_chr) %>%
        rowwise() %>%
        mutate(seq = ifelse(Strand == "+" & End < length(fasta[names(fasta) == used_chr][[1]]), as.character(fasta[names(fasta) == used_chr][[1]][Start:End]),
                            ifelse(Strand == "-" & End < length(fasta[names(fasta) == used_chr][[1]]),as.character(reverseComplement(fasta[names(fasta) == used_chr][[1]][Start:End])), NA)))
      listofdfs[[i]] <- df
    }
    
    full_table_big_data <- data.frame(Reduce(rbind, listofdfs))
    big_data_all_seq_names <- left_join(full_table_big_data, gff_table, by = c("GeneID" = "id_name"))
    return(big_data_all_seq_names)
  }
  
  # > single read table output
  if(output == "read_ids"){
    
    assigned_features <- do.call(rbind,datalist_input)
    
    # >calculate correct start end end positions of mapped reads | read in BAM file with NM tag
    allReads <- readGAlignments(input_bam_file, use.names = T, param = ScanBamParam(tag=c("NM"), what="mapq"))
    
    allReads_table <- GenomicAlignments::as.data.frame(allReads) %>%
      mutate(minion_read_name = names(allReads))# %>%
    #left_join(summary_table, by = c("minion_read_name" = "read_id")) 
    
    # > calculate number of aligned reads based on CIGAR operations (M,I)
    allReads_table$aligned_reads <- NA
    allReads_table$aligned_reads <- unlist(lapply(explodeCigarOpLengths(allReads_table$cigar, ops = c("M", "I")), function(x) sum(x)))
    
    # > join featurecounts table | calculate mapping identity | factorize strands
    allReads_table_filtered <- allReads_table %>%
      left_join(assigned_features, by = c("minion_read_name" = "id")) %>%
      mutate(identity = (1 - NM/aligned_reads)*100,
             length_read = qwidth)  %>%
      tidyr::separate_rows(gene, sep = ",") %>%
      left_join(gff_table, by = c("gene" = "id_name")) %>%
      mutate(shortest_distance_to_gene = ifelse(abs(start-start_gene) <= max(qwidth) & abs(end-end_gene)<=max(qwidth), T, F)) %>%
      dplyr::filter(shortest_distance_to_gene == T) %>%
      mutate(strand = factor(strand, levels = c("+","-"))) %>%
      dplyr::filter(strand == strand_gene) %>%
      dplyr::select(seqnames, strand, qwidth, start, end, width, mapq, NM, minion_read_name, 
                    #sequence_length_template, mean_qscore_template,
                    aligned_reads, identity, gene, mapped_type, length_read,
                    start_gene, end_gene, locus_name, strand_gene) 
    
    return(allReads_table_filtered)
  }
}

# > for BAM plotting (output from featurecounts)
full_counts_table <- counts_wrapper(input_bam_file = input_bam,
                                    input_fasta_file = input_fasta,
                                    input_gff_file = input_gff)

# > for BAM plotting (identity calculation, ...)
full_id_table <- wrapper_bam_to_table(input_bam_file = input_bam,
                                      input_gff_file = input_gff,
                                      input_fasta_file = input_fasta,
                                      datalist_input = full_counts_table[[2]],
                                      output = "read_ids")





# Coverage
depth_plus <- read.table(paste0("results/coverage_data/", sample, "_plus_depth.txt"),
                         sep = "\t", header = FALSE)
# depth_plus <- read.table("/home/alvaro/Navarrabiomed/Operon_finder/Baculovirus/results/coverage_data/ERR2300649_plus_depth.txt",
#                            sep = "\t", header = FALSE)
#                          


depth_minus <- read.table(paste0("results/coverage_data/", sample, "_minus_depth.txt"),
                          sep = "\t", header = FALSE)
# depth_minus <- read.table("/home/alvaro/Navarrabiomed/Operon_finder/Baculovirus/results/coverage_data/ERR2300649_minus_depth.txt",
#                          sep = "\t", header = FALSE)

gr_depth_plus <- GRanges(depth_plus$V1,
                         IRanges(start = depth_plus$V2, width = 1),
                         strand = "+",
                         depth = depth_plus$V3)

gr_depth_minus <- GRanges(depth_minus$V1,
                          IRanges(start = depth_minus$V2, width = 1),
                          strand = "-",
                          depth = depth_minus$V3)

gr_depth <- c(gr_depth_plus, gr_depth_minus)

rm(gr_depth_plus, gr_depth_minus)


gr_reads_all <- GRanges(full_id_table$seqnames,
                        IRanges(full_id_table$start, full_id_table$end),
                        strand = full_id_table$strand)
mcols(gr_reads_all) <- full_id_table[,c("minion_read_name", "identity")]

# Filter below 80% identity ####
gr_reads <- gr_reads_all[gr_reads_all$identity >= 80]


gff_keep <- gff[gff$type %in% interesting_list]



#First narrow down to any overlap with CDS or tRNA (middle of gene, not just edge)
fo <- findOverlaps(gr_reads, gff_keep,
                   type = "any", ignore.strand = FALSE)

gr_reads <- gr_reads[unique(queryHits(fo))]






# TSS for each gene ####
# get the annotated TSS and first 20 nucleotides; reads must overlap this
first20 <- flank(gff_keep, -20, start = TRUE)
fo <- findOverlaps(gr_reads_all, first20, type = "any", ignore.strand = FALSE)
# Get the median start site for all overlapping reads
gene_tss <-
  sapply(seq_along(gff_keep),
         function(x){
           thesereads <- gr_reads_all[queryHits(fo)[subjectHits(fo) == x]]
           if(as.character(strand(gff_keep[x])) == "-"){
             return(median(end(thesereads)))
           } else {
             return(median(start(thesereads)))
           }
         })

# compare to annotated TSS
anno_tss <- start(gff_keep)
negstrand <- as.character(strand(gff_keep)) == "-"
anno_tss[negstrand] <- end(gff_keep)[negstrand]

# TTS for each gene ####
# get the annotated TSS and last 20 nucleotides; reads must overlap this
last20 <- flank(gff_keep, -20, start = FALSE)
fo <- findOverlaps(gr_reads_all, last20, type = "any", ignore.strand = FALSE)
# Get the median start site for all overlapping reads
gene_tts <-
  sapply(seq_along(gff_keep),
         function(x){
           thesereads <- gr_reads_all[queryHits(fo)[subjectHits(fo) == x]]
           if(as.character(strand(gff_keep[x])) == "-"){
             return(median(start(thesereads)))
           } else {
             return(median(end(thesereads)))
           }
         })

anno_tts <- end(gff_keep)
anno_tts[negstrand] <- start(gff_keep)[negstrand]

# Export gene info ####
gene_out <- data.frame(refseq = seqnames(gff_keep),
                       gene_id = gff_keep$gene_id,
                       type = gff_keep$type,
                       strand = as.character(strand(gff_keep)),
                       annotated_tss = anno_tss,
                       annotated_tts = anno_tts,
                       median_read_start = round(gene_tss),
                       median_read_end = round(gene_tts))

gene_out$median_coverage <-
  sapply(seq_along(gff_keep),
         function(x){
           median(subsetByOverlaps(gr_depth, gff_keep[x])$depth)
         })





if (mode == "both" || mode == "Excludon"){
  ########################################
  #Identifying Overlapping UTRS
  ########################################
  print("################ IDENTIFYING OVERLAPING UTRS ################")
  
  final <- nrow(gene_out) - 3
  ## 5' UTRS
  index_5utr <- list()
  for (i in 1:final){
    if (gene_out[i,]$strand == "-"  & gene_out[i+1,]$strand == "+"){
      value <- gene_out[i,]$median_read_start - gene_out[i+1,]$median_read_start 
      if (value > 20 & is.na(value) == FALSE){
        index_5utr <<- append(index_5utr, row.names(gene_out[i,]))
        index_5utr <<- append(index_5utr, row.names(gene_out[i+1,]))
      }
    }
  }
  UTR5_df <- DataFrame()
  UTR5_df <- gene_out[unlist(index_5utr),]
  UTR5_df$Overlapping_type <- "5'UTR"
  UTR5_df$Pair <- rep(1:(nrow(UTR5_df)/2), each = 2)
  
  ## 3' UTRS
  index_3utr <- list()
  for (i in 1:final){
    if (gene_out[i,]$strand == "+"  & gene_out[i+1,]$strand == "-"){
      value <- gene_out[i,]$median_read_end - gene_out[i+1,]$median_read_end 
      if (value > 20 & is.na(value) == FALSE){
        index_3utr <<- append(index_3utr, row.names(gene_out[i,]))
        index_3utr <<- append(index_3utr, row.names(gene_out[i+1,]))
      }
    }
  }
  
  UTR3_df <- DataFrame()
  UTR3_df <- gene_out[unlist(index_3utr),]
  UTR3_df$Overlapping_type <- "3'UTR"
  UTR3_df$Pair <- rep(1:(nrow(UTR3_df)/2), each = 2)
  
  print(" #### DONE ####")
  write.csv(UTR3_df, file = paste0("results/operons/", sample, "Overlapping_3UTRs.csv"),
            sep = "\t", row.names = FALSE)
  
  write.csv(UTR5_df, file = paste0("results/operons/", sample, "Overlapping_5UTRs.csv"),
            sep = "\t", row.names = FALSE)
  
  
  
  
  
  
  
}


if (mode == "both" || mode == "Operon"){
  upstream_length <- 20
  print(paste("SPLIT RATIO: ", split_ratio))
  window_size <- 3
  # Merge overlapping reads ####
  gr_collapse <-  IRanges::reduce(gr_reads, ignore.strand = FALSE, min.gapwidth = 0L)
  
  nOperon <- length(gr_collapse)
  
  # Get gene IDs for each operon ####
  fo <- findOverlaps(gff_keep, gr_collapse,
                     ignore.strand = FALSE, type = "within")
  
  
  operon_genes <-
    GRangesList(lapply(seq_len(nOperon),
                       function(x){
                         ht <- queryHits(fo)[subjectHits(fo) == x]
                         gff_keep[ht]
                       }))
  # Determine whether to split operons based on depth ####
  # Loop through operons to see if they should be split
  
  
  #### set up parallel processing  #####
  n.cores <- parallel::detectCores() - 1
  #create the cluster
  my.cluster <- parallel::makeCluster(
    n.cores, 
    type = "PSOCK"
    
  )
  #register it to be used by %dopar%
  doParallel::registerDoParallel(cl = my.cluster)
  registerDoSNOW(my.cluster)
  pb <- txtProgressBar(max = nOperon, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  gr_operons_split <- GRanges()
  was_split <- logical(nOperon)
  split_last <- FALSE
  
  
  print("######### IDENTIFYING OPERONS #########")
  x <- foreach(
    i = 1:nOperon, 
    .combine = 'c', 
    .options.snow = opts,
    .packages=c("GenomicRanges", "GenomicAlignments")
    
    
  ) %dopar% {
    thesegenes <- operon_genes[[i]]
    ngenes <- length(thesegenes)
    thischr <- seqnames(gr_collapse[i])
    thisstrand <- as.character(strand(gr_collapse[i]))
    # order genes by TTS and initialize operon 5' end.
    # For genes on the + strand, the 5' end is the start.
    # On the - strand, the 5' end is the end.
    if(thisstrand == "-"){
      thesegenes <- thesegenes[order(start(thesegenes), decreasing = TRUE)]
      curr5 <- end(gr_collapse[i])
    } else {
      thesegenes <- thesegenes[order(end(thesegenes), decreasing = FALSE)]
      curr5 <- start(gr_collapse[i])
    }
    # look 20 nt (upstream_length) upstream, and compare depth to downstream through TTS of next gene
    tts_upstream <- flank(thesegenes, -upstream_length, start = FALSE)
    for(j in seq_len(ngenes)){
      if(thisstrand == "-"){
        next_gene_ends <- start(gff_keep)[start(gff_keep) < start(thesegenes[j]) & as.logical(strand(gff_keep) == "-")]
        if(length(next_gene_ends) == 0){
          this_end <- start(thesegenes[j]) - 100L
        } else {
          this_end <- max(next_gene_ends)
        }
        intergenic_space <- GRanges(seqnames(thesegenes[j]),
                                    IRanges(end = start(thesegenes[j]) - 1L,
                                            start = this_end),
                                    strand = "-")
      } else {
        next_gene_ends <- end(gff_keep)[end(gff_keep) > end(thesegenes[j]) & as.logical(strand(gff_keep) == "+")]
        if(length(next_gene_ends) == 0){
          this_end <- end(thesegenes[j]) + 100L
        } else {
          this_end <- min(next_gene_ends)
        }
        intergenic_space <- GRanges(seqnames(thesegenes[j]),
                                    IRanges(start = end(thesegenes[j]) + 1L,
                                            end = this_end),
                                    strand = "+")
      }
      
      depth_downstream <- subsetByOverlaps(gr_depth, intergenic_space)$depth
      depth_upstream <- mean(subsetByOverlaps(gr_depth, tts_upstream[j])$depth)
      # Determine whether to split after this gene
      base_fail <- depth_downstream < depth_upstream / split_ratio
      base_rle <- rle(base_fail)
      if(any(base_fail)){
        max_run <- max(base_rle$lengths[base_rle$values])
      } else {
        max_run <- 0
      }
      # Perform split
      if(max_run >= window_size || depth_upstream == 0){
        # Set 3' end of new operon to TTS of this gene
        curr3 <- ifelse(thisstrand == "-", start(thesegenes[j]), end(thesegenes[j]))
        # Set 5' end of new operon to TSS of this gene, if it isn't further upstream.
        curr5 <- ifelse(thisstrand == "-", max(curr5, end(thesegenes[j])),
                        min(curr5, start(thesegenes[j])))
        # Add operon to list (if it fully contains CDS, which it may not with overlapping CDS)
        newgr <- GRanges(thischr,
                         IRanges(start = min(curr5, curr3),
                                 end = max(curr5, curr3)),
                         thisstrand)
        if(length(findOverlaps(gff_keep, newgr, type = "within")) > 0){
          gr_operons_split <- c(gr_operons_split, newgr)
        }
        # Set the 5' end of the next operon immediately after 3' end of this one.
        curr5 <- ifelse(thisstrand == "-", curr3 - 1, curr3 + 1)
        # Record that this operon was split
        was_split[i] <- TRUE
        # If we "split" after the last gene, we don't want to add another operon
        if(j == ngenes) split_last <- TRUE
      }
    }
    # Add the remaining transcript (or whole transcript if not split)
    # assuming any genes are left on it.
    curr3 <- ifelse(thisstrand == "-", start(gr_collapse[i]), end(gr_collapse[i]))
    newgr <- GRanges(thischr,
                     IRanges(start = min(c(curr5, curr3)),
                             end = max(c(curr5, curr3))),
                     thisstrand)
    if(length(findOverlaps(newgr, thesegenes)) > 0 && !split_last){
      if(length(findOverlaps(gff_keep, newgr, type = "within")) > 0){
        gr_operons_split <- c(gr_operons_split, newgr)
      }
    }
    
    split_last <- FALSE
    gr_operons_split
    
  }
  
  x <- unique(x)
  gr_operons_split <- x
  # stop cluster
  stopCluster(my.cluster)

  
  # Detect discontiguous operons after split ####
  fo <- findOverlaps(gff_keep, gr_operons_split,
                     ignore.strand = FALSE, type = "within")
  
  operon_split_genes <-
    GRangesList(lapply(seq_along(gr_operons_split),
                       function(x){
                         ht <- queryHits(fo)[subjectHits(fo) == x]
                         gff_keep[ht]
                       }))
  
  operon_split_gene_ranges <- unlist(range(operon_split_genes))
  
  fo <- findOverlaps(operon_split_gene_ranges, ignore.strand = TRUE, type = "within")
  
  fo6 <- fo[queryHits(fo) != subjectHits(fo),]
  fo <- findOverlaps(gff_keep, gr_operons_split,
                     type = "within", ignore.strand = FALSE)
  
  stopifnot(anyDuplicated(queryHits(fo)) == 0)
  stopifnot(length(setdiff(seq_along(gr_operons_split), subjectHits(fo))) == 0)
  
  gene_out$operon <- subjectHits(fo)[match(seq_along(gff_keep), queryHits(fo))]
  tempstart <- start(gr_operons_split)[gene_out$operon]
  tempend <- end(gr_operons_split)[gene_out$operon]
  gene_out$operon_start[gene_out$strand == "+"] <- tempstart[gene_out$strand == "+"]
  gene_out$operon_end[gene_out$strand == "+"] <- tempend[gene_out$strand == "+"]
  gene_out$operon_start[gene_out$strand == "-"] <- tempend[gene_out$strand == "-"]
  gene_out$operon_end[gene_out$strand == "-"] <- tempstart[gene_out$strand == "-"]
  #gene_out$discontiguous_operon <- gene_out$operon %in% subjectHits(fo6)
  
  
  #######################################
  #Check if one read cover flanking genes
  #######################################
  
  list_index <- list()
  operon_num <- unique(gene_out$operon)
  operon_num <- na.omit(operon_num)
  for ( i in operon_num){
    genes <- which(gene_out$operon == i)
    gene_1 <- genes[1]
    gene_2 <- genes[length(genes)]
    keep_index <- gene_1:gene_2
    list_strand <- gene_out[gene_1:gene_2,]$strand
    if (gene_out[gene_1,]$strand == "+" && gene_out[gene_2,]$strand == "+" && (sum(list_strand == "-") >= 1)){
      list_index <<- append(list_index, keep_index)
    }else if (gene_out[gene_1,]$strand == "-" && gene_out[gene_2,]$strand == "-" && (sum(list_strand == "+") >= 1)){
      list_index <<- append(list_index, keep_index)
    }
  }
  
  df_filtered <- gene_out[unlist(list_index),]
  operones <- df_filtered$operon[duplicated(df_filtered$operon)]
  operones <- unique(operones)
  operones <- na.omit(operones)
  for (i in operones){
    genes <- which(df_filtered$operon == i)
    if (identical(genes, integer(0)) == F){
      gene_1 <- genes[1]
      gene_2 <- genes[length(genes)]
      list_genes <- df_filtered[gene_1:gene_2,]$gene_id
      list_genes <- na.omit(list_genes)
      if (any(duplicated(list_genes))) {
        duplicated_genes <- list_genes[duplicated(list_genes)][1]
        position <- which(df_filtered$gene_id == duplicated_genes)[2]
        df_filtered$operon[gene_1:(position-1)] <- i
        pos <- which(operones == i)
        df_filtered$operon[position:gene_2] <- operones[pos+1]
        df_filtered$operon_start[df_filtered$operon == i] <- df_filtered$annotated_tss[gene_1]
        df_filtered$operon_end[df_filtered$operon == i] <- df_filtered$annotated_tts[gene_2]
      }else if ( length(genes) > 1){
        df_filtered$operon[gene_1:gene_2] <- i
        df_filtered$operon_start[df_filtered$operon == i] <- df_filtered$annotated_tss[gene_1]
        df_filtered$operon_end[df_filtered$operon == i] <- df_filtered$annotated_tts[gene_2]
      }
      
    }
  }
  
  
  
  
  
  
  
  operones <- df_filtered$operon[duplicated(df_filtered$operon)]
  operones <- unique(operones)
  operones <- na.omit(operones)
  ## Just one reverse gene with minus strand
  for (i in operones){
    ##Get first and last gene of operon
    first_gene <- which(df_filtered$operon == i)[1]
    last_gene <- which(df_filtered$operon == i)[length(which(df_filtered$operon == i))]
    ### Get genes flanking the reverse gene
    operon <- df_filtered[first_gene:last_gene,]
    if (sum(operon$strand == "+") >= 2 && sum(operon$strand == "-") == 1 ){
      reverse_gene = which(operon$strand == "-")
      flanking_gene1 <- operon[reverse_gene -1,]$gene_id
      flanking_gene2 <- operon[reverse_gene +1,]$gene_id
      ###Obtain star and ends for these flanking genes
      df1 <- as.data.frame(operon_genes)
      f_gene1_end <- df1[which(df1$gene_id == flanking_gene1),]$end
      f_gene2_start <- df1[which(df1$gene_id == flanking_gene2),]$start
      reads <- full_id_table[full_id_table$start < f_gene1_end -20, ]
      readsfinal <- nrow(reads[reads$end > f_gene2_start +20,])
      position <- which(df_filtered$operon == i)[1]
      df_filtered[position, "Nº_Overlapping_reads*"] <- readsfinal
      
    }
  }
  
  ## Just one reverse gene with plus strand
  for (i in operones){
    ##Get first and last gene of operon
    first_gene <- which(df_filtered$operon == i)[1]
    last_gene <- which(df_filtered$operon == i)[length(which(df_filtered$operon == i))]
    ### Get genes flanking the reverse gene
    operon <- df_filtered[first_gene:last_gene,]
    if (sum(operon$strand == "-") >= 2 && sum(operon$strand == "+") == 1 ){
      reverse_gene = which(operon$strand == "+")
      flanking_gene1 <- operon[reverse_gene +1,]$gene_id
      flanking_gene2 <- operon[reverse_gene -1,]$gene_id
      ###Obtain star and ends for these flanking genes
      df1 <- as.data.frame(operon_genes)
      f_gene1_end <- df1[which(df1$gene_id == flanking_gene1),]$end
      f_gene2_start <- df1[which(df1$gene_id == flanking_gene2),]$start
      reads <- full_id_table[full_id_table$start < f_gene1_end +20, ]
      readsfinal <- nrow(reads[reads$end > f_gene2_start -20,])
      position <- which(df_filtered$operon == i)[1]
      df_filtered[position, "Nº_Overlapping_reads*"] <- readsfinal
    }
  }
  
  ## Several reverses genes with plus strand
  for (i in operones){
    #Get first and last gene of operon
    first_gene <- which(df_filtered$operon == i)[1]
    last_gene <- which(df_filtered$operon == i)[length(which(df_filtered$operon == i))]
    ### Get genes flanking the reverse gene
    operon <- df_filtered[first_gene:last_gene,]
    if (sum(operon$strand == "+") >= 2 && sum(operon$strand == "-") >= 2 ){
      plus_strand <- sum(operon$strand == "+")
      minus_strand <- sum(operon$strand == "-")
      if (plus_strand < minus_strand){
        x <- which(operon$strand == "+")
        reverse_gene <- split(x, cumsum(c(1, diff(x) != 1)))
        for (j in reverse_gene){
          if (length(j) == 1){
            flanking_gene1 <- operon[1 +1,]$gene_id
            flanking_gene2 <- operon[1 -1,]$gene_id
            if (identical(flanking_gene1, character(0) ) == FALSE && identical(flanking_gene2, character(0) ) == FALSE){
              ###Obtain star and ends for these flanking genes
              df1 <- as.data.frame(operon_genes)
              f_gene1_end <- df1[which(df1$gene_id == flanking_gene1),]$end
              f_gene2_start <- df1[which(df1$gene_id == flanking_gene2),]$start
              reads <- full_id_table[full_id_table$start < f_gene1_end +20, ]
              readsfinal <- nrow(reads[reads$end > f_gene2_start -20,])
              position <- which(df_filtered$operon == i)[1]
              df_filtered[position, "Nº_Overlapping_reads*"] <- readsfinal
            }
          }else if (length(j) > 1){
            first_reverse_gene <- operon[j[1],]$gene_id
            last_reverse_gene <- operon[j[length(j)],]$gene_id
            flanking_position1 <- which(df_filtered$gene_id == last_reverse_gene) +1
            flanking_position2 <- which(df_filtered$gene_id == first_reverse_gene) -1
            f_gene1 <- df_filtered[flanking_position1,]$gene_id
            f_gene2 <- df_filtered[flanking_position2,]$gene_id
            if (identical(f_gene1, character(0) ) == FALSE && identical(f_gene2, character(0) ) == FALSE){
              f_gene1_end <- df1[which(df1$gene_id == f_gene2),]$start
              f_gene2_start <- df1[which(df1$gene_id == f_gene1),]$start
              reads <- full_id_table[full_id_table$start < f_gene1_end +20, ]
              readsfinal <- nrow(reads[reads$end > f_gene2_start -20,])
              position <- which(df_filtered$operon == i)[1]
              df_filtered[position, "Nº_Overlapping_reads*"] <- readsfinal
            }
          }
        }
      }
    }
  }
  
  
  ## Several reverse genes with minus strand
  for (i in operones){
    #Get first and last gene of operon
    first_gene <- which(df_filtered$operon == i)[1]
    last_gene <- which(df_filtered$operon == i)[length(which(df_filtered$operon == i))]
    ### Get genes flanking the reverse gene
    operon <- df_filtered[first_gene:last_gene,]
    if (sum(operon$strand == "+") >= 2 && sum(operon$strand == "-") >= 2 ){
      plus_strand <- sum(operon$strand == "+")
      minus_strand <- sum(operon$strand == "-")
      if (plus_strand > minus_strand){
        x <- which(operon$strand == "-")
        reverse_gene <- split(x, cumsum(c(1, diff(x) != 1)))
        for (j in reverse_gene){
          if (length(j) == 1){
            flanking_gene1 <- operon[1 -1,]$gene_id
            flanking_gene2 <- operon[1 +1,]$gene_id
            if (identical(flanking_gene1, character(0) ) == FALSE && identical(flanking_gene2, character(0) ) == FALSE){
              ###Obtain star and ends for these flanking genes
              df1 <- as.data.frame(operon_genes)
              f_gene1_end <- df1[which(df1$gene_id == flanking_gene1),]$end
              f_gene2_start <- df1[which(df1$gene_id == flanking_gene2),]$start
              reads <- full_id_table[full_id_table$start < f_gene1_end -20, ]
              readsfinal <- nrow(reads[reads$end > f_gene2_start +20,])
              position <- which(df_filtered$operon == i)[1]
              df_filtered[position, "Nº_Overlapping_reads*"] <- readsfinal
            }
          }else if (length(j) > 1){
            first_reverse_gene <- operon[j[1],]$gene_id
            last_reverse_gene <- operon[j[length(j)],]$gene_id
            flanking_position1 <- which(df_filtered$gene_id == first_reverse_gene) -1
            flanking_position2 <- which(df_filtered$gene_id == last_reverse_gene) +1
            f_gene1 <- df_filtered[flanking_position1,]$gene_id
            f_gene2 <- df_filtered[flanking_position2,]$gene_id
            if (identical(f_gene1, character(0) ) == FALSE && identical(f_gene2, character(0) ) == FALSE){
              f_gene1_end <- df1[which(df1$gene_id == f_gene2),]$start
              f_gene2_start <- df1[which(df1$gene_id == f_gene1),]$start
              reads <- full_id_table[full_id_table$start < f_gene1_end -20, ]
              readsfinal <- nrow(reads[reads$end > f_gene2_start +20,])
              position <- which(df_filtered$operon == i)[1]
              df_filtered[position, "Nº_Overlapping_reads*"] <- readsfinal
            }
          }
        }
      }
    }
  }
  
  # Export operons to csv file ####
  if (nrow(df_filtered) == 0){
    print("###### ANY NON-CONTIGUOUS OPERON WAS IDENTIFIED ######")
  } else{
    write.csv(df_filtered, file = paste0("results/operons/", sample, "NonContiguous_operons.csv"),
              sep = "\t", row.names = FALSE)
    
    #   # Export operons to BED file ####
    #   list_start <- list()
    #   list_end <- list()
    #   list_resfseq <- list()
    #   list_strand <- list()
    #   for ( i in unique(df_filtered$operon)){
    #     i <-50
    #     first_gene <- which(df_filtered$operon == i)[1]
    #     last_gene <- which(df_filtered$operon == i)[length(which(df_filtered$operon == i))]
    #     start <- df_filtered$operon_start[first_gene]
    #     end <- df_filtered$operon_end[last_gene]
    #     seqname <- df_filtered$refseq[first_gene]
    #     strand <- df_filtered$strand[first_gene]
    #     
    #     list_start <<- append(list_start, start)
    #     list_end <<- append(list_end, end)
    #     list_resfseq <<- append(list_resfseq, seqname)
    #     list_strand <<- append(list_strand, strand)
    #   }
    #   
    #   df_NCoperons <- data.frame(chr=unlist(list_resfseq), start=unlist(list_start), end=unlist(list_end),
    #                              strand=unlist(list_strand))
    #   df_NCoperons_plus <- df_NCoperons[df_NCoperons$strand == "+", ]
    #   # subset the rows where strand == "-"
    #   swap_rows <- df_NCoperons[df_NCoperons$strand == "-", ]
    #   # swap the values in columns x and y in the selected rows
    #   swap_rows[, c("start", "end")] <- swap_rows[, c("end", "start")]
    #   
    #   df_NCoperons <- rbind(df_NCoperons_plus, swap_rows)
    #   gr <- makeGRangesFromDataFrame(df_NCoperons, keep.extra.columns=TRUE)
    #   rtracklayer::export.bed(gr,
    #                           con = paste0("results/operons/", sample, "_operons_PRUEBA",".bed"))
    #   
  }
}

if (mode == "both"){
  Noperon <- length(unique(df_filtered$operon))
  cat("###### ANALYSIS FINISHED ######\n\nNumber of reads identified\n\nNon-contiguous Operons: ",
      Noperon, "\nNumber of Overlapping 3' UTRs : ", nrow(UTR3_df)/2,
      "\nNumber of Overlapping 5' UTRs : ", nrow(UTR5_df)/2)
}

if (mode == "Operon"){
  Noperon <- length(unique(df_filtered$operon))
  cat("###### ANALYSIS FINISHED ######\n\nNumber of reads identified\n\nNon-contiguous Operons: ",
      Noperon)
}

if (mode == "Excludon"){
  cat("###### ANALYSIS FINISHED ######\n\nNumber of Overlapping 3' UTRs : ", nrow(UTR3_df)/2,
      "\nNumber of Overlapping 5' UTRs : ", nrow(UTR5_df)/2)
}

options(warn = defaultW)
