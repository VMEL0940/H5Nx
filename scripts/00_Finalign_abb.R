### Libraries and functions
library(Biostrings)
library(DECIPHER)
library(stringr)
library(dplyr)

cons_matrix <- function(filtered_seqs) {
  if (length(filtered_seqs) == 0) {
    return(DNAStringSet(character(0)))  # Return an empty DNAStringSet
  }
  
  cons_matrix <- consensusMatrix(filtered_seqs)
  consensus <- apply(cons_matrix, 2, function(x) names(which.max(x)))
  consensus_string <- paste(consensus, collapse = "")
  string_consensus <- consensus_string %>% DNAStringSet()
  return(string_consensus)
}

process_seq <- function(seq, consensus_string) {
  
  log <- character()
  numbers <- character()
  
  seq <- gsub("-", "N", as.character(seq)) %>% DNAStringSet()
  result <- data.frame(sequence_number = integer(), position = integer(), change = character())
  
  for (i in 1:length(seq)) {
    cat("\r", paste0("Repairing single errors ", i, "/", length(seq)), sep = "")
    
    gap_positions <- which(strsplit(as.character(seq[[i]]), "")[[1]] == "N")
    
    extracted_gap <- c()
    gap_diff <- diff(gap_positions)
    
    
    
    if(length(gap_positions) == 1){
      extracted_gap <- append(extracted_gap, gap_positions[1])
    }else if(length(gap_positions)>1){
      if(length(gap_diff) == 1){
        if(gap_diff[1] > 1){
          extracted_gap <- append(extracted_gap, gap_positions[1])
          extracted_gap <- append(extracted_gap, gap_positions[2])
        }
      }else if(length(gap_diff) > 1){
        if(gap_diff[1] > 1){
          extracted_gap <- append(extracted_gap, gap_positions[1])
        } 
        for (j in 2:(length(gap_positions)-1)) {
          if((gap_diff[j-1] != 1) & (gap_diff[j] != 1)) {
            extracted_gap <- append(extracted_gap, gap_positions[j])
          }
        }
        if (gap_diff[length(gap_diff)] > 1){
          extracted_gap <- append(extracted_gap, gap_positions[length(gap_positions)])
        }
      }
    }
    
    for (k in extracted_gap) {
      if (as.character(seq[[i]][k]) != substr(consensus_string, k, k)) {
        change <- paste0(seq[[i]][k])
        seq[[i]][k] <- substr(consensus_string, k, k)
        log <- append(log, paste0(change, " -> ", substr(consensus_string, k, k)))
        numbers <- append(numbers, k)
        
        tmp_df <- data.frame(sequence_number = i, position = k, change = paste0(change, " -> ", substr(consensus_string, k, k)))
        
        # result <- rbind(result, tmp_df)
      }
    }
  }
  
  return(list(seq = seq, num = i, log = unique(log), numbers = unique(numbers), result = result))
}

first_last_fill <- function(seq, consensus_string) {
  
  first_fills <- c()
  last_fills <- c()
  log <- character()
  numbers <- character()
  
  for (i in 1:length(seq)) {
    cat("\r", paste0("Repairing single errors ", i, "/", length(seq)), sep = "")
    
    gap_positions <- which(strsplit(as.character(seq[[i]]), "")[[1]] == "N")
    extracted_gap <- c()
    gap_diff <- diff(gap_positions)
    
    
    if(length(gap_diff) != 0) {
      firstgap <- find_first_one(gap_diff)
      lastgap <- find_last_one(gap_diff) + 1
      
      gap_list <- Map(base::c, firstgap, lastgap)
      
      for(l in 1:length(gap_list)){
        if(length(gap_list[[l]]) != 0){
          if((gap_positions[as.numeric(gap_list[[l]][1])] %% 3 == 0) & (gap_positions[as.numeric(gap_list[[l]][2])] %% 3 == 1)) {
            first_fills <- append(first_fills, gap_positions[as.numeric(gap_list[[l]][1])])
            last_fills <- append(first_fills, gap_positions[as.numeric(gap_list[[l]][2])])
          } 
          
          
          if((length(first_fills) != 0) & (length(last_fills) != 0)){
            for (k in first_fills){ 
              change <- paste0(k, " (", seq[[i]][k], ")")
              seq[[i]][k] <- substr(consensus_string, k, k) 
              log <- append(log, paste0(change, " -> ", substr(consensus_string, k, k)))
              numbers <- append(numbers, k) }
            
            for (k in last_fills) { 
              change <- paste0(k, " (", seq[[i]][k], ")")
              seq[[i]][k] <- substr(consensus_string, k, k) 
              log <- append(log, paste0(change, " -> ", substr(consensus_string, k, k)))
              numbers <- append(numbers, k) }
            
          }
        }
      }
    }
  }
  return(list(seq = seq, log = unique(log), numbers = unique(numbers)))
  # print(list(log, numbers))
}

find_first_one <- function(x) {
  r <- rle(x)
  starts <- c(1, head(cumsum(r$lengths), -1) + 1)
  ones <- which(r$values == 1)
  starts[ones]
}

find_last_one <- function(x) {
  r <- rle(x)
  ends <- cumsum(r$lengths)
  ones <- which(r$values == 1)
  ends[ones]
}



adjust_sequence <- function(seq, consensus_string) {
  type <- character()
  log <- character()
  numbers <- character()
  
  result <- data.frame(sequence = integer(), aa_position = character(), type = character(), change = character())
  
  pattern <- "[ACGTN]"
  seq <- gsub("N", "-", as.character(seq))
  
  triplets <- strsplit(as.character(seq), "(?<=.{3})", perl = TRUE)
  consensus_triplets <- strsplit(as.character(consensus_string), "(?<=.{3})", perl = TRUE)   # ?<= : positive lookbehind, .{3}: repeated exactly three times
  index_con_gap <- which(consensus_triplets[[1]] == "---")
  
  i <- 1 # Initialize i here
  n <- NA
  
  while(i <= length(triplets) ) {
    j <- 2 # initialize j here
    
    cat("\r", paste0("Processing triplet ", i, "/", length(triplets)), sep = "")
    
    while (j < length(unlist(triplets[i]))) {
      
      if(str_count(triplets[[i]][j], "[^\\-]") == 1){
        
        n <- find_nearest_non_zero(j, triplets[[i]]) %>% as.numeric()
        # print(paste0("n:", n))
        
        
        if(!is.na(n)) {
          count <- str_count(triplets[[i]][n], "[^\\-]")  # Calculate non-hyphen character count
          
          if(count == 1){
            type <-"1-1 -> 0-0"
            change <- paste0(triplets[[i]][j], " ~ ", triplets[[i]][n])
            triplets[[i]][n] <- "---"
            triplets[[i]][j] <- "---"
            log <- paste0(change, " --> ", triplets[[i]][j], " ~ ", triplets[[i]][n])
            numbers <- paste0(j, " ~ ", n)
            tmp_df <- data.frame(sequence = i, aa_position = numbers, type = type, change = log)
            result <- rbind(result, tmp_df)
            # j <- j - 1
          }
          
          else if(count == 2){
            change <- paste0(triplets[[i]][j], " ~ ", triplets[[i]][n])
            nucleotides_front1 <- paste0(str_extract_all(triplets[[i]][j], pattern)[[1]][1])
            nucleotides_back2 <- paste(str_extract_all(triplets[[i]][n], pattern)[[1]][c(1,2)], collapse = "")
            nucleotide_string <- paste0(nucleotides_front1, nucleotides_back2)
            
            if(nucleotide_string %in% c("TAG", "TGA", "TAA")){
              type <- "1-2 -> 0-3"
              triplets[[i]][j] <- "---" #delete front
              triplets[[i]][n] <- consensus_triplets[[1]][n] #fill back            
              log <- paste0(change, " --> ", triplets[[i]][j], " ~ ", triplets[[i]][n])
              numbers <- paste0(j, " ~ ", n)
              tmp_df <- data.frame(sequence = i, aa_position = numbers, type = type, change = log)
              result <- rbind(result, tmp_df)
              
            } else if(n-j > 2){
              type <- "1-2 -> 0-3"
              triplets[[i]][j] <- "---" #delete front
              
              extracted_nucleotides <- sapply(seq, function(seq) as.character(subseq(seq, start = n, end = n)))
              trip_i_n<- strsplit(triplets[[i]][n], "")
              
              if(any(trip_i_n[[1]] %in% extracted_nucleotides)){
                #n자리에 하나라도 일치하는 게 있으면 붙이기 
                triplets[[i]][n] <- nucleotide_string #paste 
                log <- paste0(change, " --> ", triplets[[i]][j], " ~ ", triplets[[i]][n])
                numbers <- paste0(j, " ~ ", n)
                tmp_df <- data.frame(sequence = i, aa_position = numbers, type = type, change = log)
              }else{
                triplets[[i]][n] <- consensus_triplets[[1]][n] #fill back
                log <- paste0(change, " --> ", triplets[[i]][j], " ~ ", triplets[[i]][n])
                numbers <- paste0(j, " ~ ", n)
                tmp_df <- data.frame(sequence = i, aa_position = numbers, type = type, change = log)
                result <- rbind(result, tmp_df)
              }
            } else {
              type <- "1-2 -> 0-3"
              triplets[[i]][j] <- "---"
              triplets[[i]][n] <- nucleotide_string
              log <- paste0(change, " --> ", triplets[[i]][j], " ~ ", triplets[[i]][n])
              numbers <- paste0(j, " ~ ", n)
              tmp_df <- data.frame(sequence = i, aa_position = numbers, type = type, change = log)
              result <- rbind(result, tmp_df)
            }
          } 
          
          else if(count == 3){
            type <- "1-3 -> 0-3"
            change <- paste0(triplets[[i]][j], " ~ ", triplets[[i]][n])
            triplets[[i]][j] <- '---'
            triplets[[i]][n] <- triplets[[i]][n]
            log <- paste0(change, " --> ", triplets[[i]][j], " ~ ", triplets[[i]][n])
            numbers <- paste0(j, " ~ ", n)
            tmp_df <- data.frame(sequence = i, aa_position = numbers, type = type, change = log)
            result <- rbind(result, tmp_df)
          }
          else {
            type <- "1-X"
            change <- paste0(triplets[[i]][j], " ~ ", triplets[[i]][n])
            str_extract_all(triplets[[i]][j], pattern)[[1]] <- str_extract_all(triplets[[i]][j], pattern)[[1]]
            str_extract_all(triplets[[i]][n], pattern)[[1]] <- str_extract_all(triplets[[i]][n], pattern)[[1]]
            log <- paste0(change, " --> ", triplets[[i]][j], " ~ ", triplets[[i]][n])
            numbers <- paste0(j, " ~ ", n)
            tmp_df <- data.frame(sequence = i, aa_position = numbers, type = type, change = log)
            result <- rbind(result, tmp_df)
          }
        } 
      } else if(str_count(triplets[[i]][j], "[^\\-]") == 2) { 
        
        n <- find_nearest_non_zero(j, triplets[[i]])  %>% as.numeric()
        # print(paste0("n:", n))
        
        if(!is.na(n)){
          count <- str_count(triplets[[i]][n], "[^\\-]")
          
          if(count == 1){
            change <- paste0(triplets[[i]][j], " ~ ", triplets[[i]][n])
            nucleotides_front2 <- paste(str_extract_all(triplets[[i]][j], pattern)[[1]][c(1,2)], collapse ="")
            nucleotides_back1 <- paste0(str_extract(triplets[[i]][n], pattern)[[1]][1])
            nucleotide_string <- paste0(nucleotides_front2, nucleotides_back1)
            
            #check nucleotide_string is NOT stop codon (TAG, TGA, TAA)
            if(nucleotide_string %in% c("TAG", "TGA", "TAA")){
              type <- "2-1 -> 3-0"
              if(grepl("N", consensus_triplets[[1]][j])) {
                triplets[[i]][j] <- consensus_triplets[[1]][j] #fill front
              } else {
                triplets[[i]][j] <- nucleotide_string
              }
              
              triplets[[i]][n] <- "---" #delete back
              log <- paste0(change, " --> ", triplets[[i]][j], " ~ ", triplets[[i]][n])
              numbers <- paste0(j, " ~ ", n)
              tmp_df <- data.frame(sequence = i, aa_position = numbers, type = type, change = log)
              result <- rbind(result, tmp_df)
              
            }else{
              if(n-j>2){
                type <- "2-1 -> 3-0"
                
                if(grepl("N", consensus_triplets[[1]][j])) {
                  triplets[[i]][j] <- consensus_triplets[[1]][j] #fill front
                } else {
                  triplets[[i]][j] <- nucleotide_string
                }
                
                # triplets[[i]][j] <- consensus_triplets[[1]][j] #fill front
                nucleotide_string
                triplets[[i]][n] <- "---" #delete back
                log <- paste0(change, " --> ", triplets[[i]][j], " ~ ", triplets[[i]][n])
                numbers <- paste0(j, " ~ ", n)
                tmp_df <- data.frame(sequence = i, aa_position = numbers, type = type, change = log)
                result <- rbind(result, tmp_df)
              }else{
                type <- "2-1 -> 3-0"
                triplets[[i]][j] <- nucleotide_string
                triplets[[i]][n] <- "---"
                log <- paste0(change, " --> ", triplets[[i]][j], " ~ ", triplets[[i]][n])
                numbers <- paste0(j, " ~ ", n)
                tmp_df <- data.frame(sequence = i, aa_position = numbers, type = type, change = log)
                result <- rbind(result, tmp_df)
              }
            }
          }
          
          else if(count== 2){
            type <- "2-3 -> 3-3"
            change <- paste0(triplets[[i]][j], " ~ ", triplets[[i]][n])
            triplets[[i]][j] <- consensus_triplets[[1]][j] 
            triplets[[i]][n] <- consensus_triplets[[1]][n] 
            log <- paste0(change, " --> ", triplets[[i]][j], " ~ ", triplets[[i]][n])
            numbers <- paste0(j, " ~ ", n)
            tmp_df <- data.frame(sequence = i, aa_position = numbers, type = type, change = log)
            result <- rbind(result, tmp_df)
          }
          
          else if(count == 3){
            type <- "2-3 -> 3-3"
            change <- paste0(triplets[[i]][j], " ~ ", triplets[[i]][n])
            triplets[[i]][n] <- triplets[[i]][n]
            triplets[[i]][j] <- consensus_triplets[[1]][j]  
            log <- paste0(change, " --> ", triplets[[i]][j], " ~ ", triplets[[i]][n])
            numbers <- paste0(j, " ~ ", n)
            tmp_df <- data.frame(sequence = i, aa_position = numbers, type = type, change = log)
            result <- rbind(result, tmp_df)
          } else {
            type <- "2-X"
            change <- paste0(triplets[[i]][j], " ~ ", triplets[[i]][n])
            str_extract_all(triplets[[i]][j], pattern)[[1]] <- str_extract_all(triplets[[i]][j], pattern)[[1]]
            str_extract_all(triplets[[i]][n], pattern)[[1]] <- str_extract_all(triplets[[i]][n], pattern)[[1]]
            log <- paste0(change, " --> ", triplets[[i]][j], " ~ ", triplets[[i]][n])
            numbers <- paste0(j, " ~ ", n)
            tmp_df <- data.frame(sequence = i, aa_position = numbers, type = type, change = log)
            result <- rbind(result, tmp_df)
          }
        } else {
          # type <- "3-X"
          # change <- paste0(triplets[[i]][j], " ~ ", triplets[[i]][n])
          # triplets[[i]][j] <- triplets[[i]][j]
          # log <- paste0(change, " --> ", triplets[[i]][j], " ~ ", triplets[[i]][n])
          # numbers <- paste0(j, " ~ ", n)
          # result <- c(result, NA)
          result <- result
          
        }
      } 
      
      j <- j + 1
      
      if(j >= length(unlist(triplets[i]))){break}
    }
    
    i <- i + 1 # Increment i at the end of the loop
    if(i > length(triplets)){ 
      break
    }
  }
  
  sorted_result <- distinct(result, aa_position, .keep_all = TRUE)
  return(list(result = sorted_result, triplets = triplets))
  
}

find_nearest_non_zero <- function(index, triplets) {
  for(n in (index+1):length(triplets)) {
    if(str_count(triplets[n], "[^\\-]") > 0) {     # [^\\-] means NOT '-'
      return(n)
    }
  }
  return(NA)  # Return NA if no non-zero triplet is found
}

### Cleaning

seq <- readDNAStringSet('/Volumes/shared610/H5NX_NEW/fasta/NA_500_modified.fas')
seq <- replaceAmbiguities(seq)
seq <- gsub('-', 'N', seq) %>% Biostrings::DNAStringSet()

d <- DistanceMatrix(seq)
cutoff <- d %>% mean

clusters <- Clusterize(seq, cutoff=cutoff, penalizeGapLetterMatches = F, invertCenters=TRUE)
apply(clusters, 2, function(x) max(abs(x)))
apply(clusters, 2, function(x) (abs(x))) %>% table %>% as.data.frame()

consensus_list <- unique(seq[clusters$cluster < 0])


consd <- DistanceMatrix(consensus_list)
conscutoff <- consd %>% mean


cons_group <- list()
for(i in seq_len(length(table(abs(clusters$cluster))))) {
  matching_indices <- abs(clusters$cluster[clusters$cluster < 0]) == i
  cons_group[[i]] <- c(seq[clusters$cluster < 0][matching_indices], seq[clusters$cluster == i])
}

consensus_list <- lapply(cons_group, cons_matrix)

for(i in length(cons_group):1) {
  if(length(cons_group[[i]]) == 0) {
    cons_group[[i]] <- NULL  
    consensus_list[[i]] <- NULL
  }
}

pro_seq <- list()
pro_seq <- mapply(function(x, y) process_seq(x, DNAStringSet(y)), cons_group, consensus_list, SIMPLIFY = FALSE)
pro_seq_seq <- lapply(pro_seq, function(x) x$seq)

pro_seq2 <- mapply(function(x, y) first_last_fill(x, DNAStringSet(y)), pro_seq_seq, consensus_list, SIMPLIFY = FALSE)
pro_seq2_seq <- lapply(pro_seq2, function(x) x$seq)


adj_seq <- mapply(function(x, y) adjust_sequence(x, y), 
                  pro_seq2_seq, 
                  consensus_list, 
                  SIMPLIFY = FALSE)

triple_auto <- lapply(adj_seq, function(x) lapply(x$triplets, paste, collapse = "")) %>% unlist %>% DNAStringSet()
names(triple_auto) <- lapply(pro_seq_seq, function(x) unlist(names(x))) %>% unlist

triple_auto_fin <- gsub('-', 'N', triple_auto) %>% Biostrings::DNAStringSet()

writeXStringSet(triple_auto_fin, '/Volumes/shared610/H5NX_NEW/fasta/NA_500_modified_finalign.fas')

