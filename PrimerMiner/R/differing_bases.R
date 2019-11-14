differing_bases <- function(files = c("List of files"), predator_name, min_prey = 0.75, max_pred = 0.25,
                            output = "plot") {


  #' Find positions where prey DNA differs from that of its predator, to aid in designing dietary primers. To be used after aligning sequences, as per https://www.youtube.com/watch?v=EcdTxGmcsEM
  #'
  #' @param files a vector of input files
  #' @param predator_name
  #' @param min_prey the minimum proportion of prey taxa which need to have a given base for it to be considered a potentially useful site for a primer
  #' @param max_prey the maximum proportion of predator sequences which need to have a given base for it to be considered a potentially useful site for a primer
  #' @return Either a plot summarising the results, or a list of tibbles with useful values

  if (is.numeric(lines)) {
    if (lines == 0) {
      deletebox <- F
    } else {
      deletebox <- T
    }
  } else {
    deletebox <- F
  }

  # Function calculates the portion of each base at each position within
  # a given file
  entropy <- function(fasta) {
    {
      upac <- read.csv(text = c("ID,comment,A,T,C,G,farbe
A,Adenine,1,0,0,0,F
C,Cytosine,0,0,1,0,F
G,Guanine,0,0,0,1,F
T,Thymine,0,1,0,0,F
R,A or G,0.5,0,0,0.5,T
Y,C or T,0,0.5,0.5,0,T
S,G or C,0,0,0.5,0.5,T
W,A or T,0.5,0.5,0,0,T
K,G or T,0,0.5,0,0.5,T
M,A or C,0.5,0,0.5,0,T
B,C or G or T,0,0.3,0.3,0.3,T
D,A or G or T,0.3,0.3,0,0.3,T
H,A or C or T,0.3,0.3,0.3,0,T
V,A or C or G,0.3,0,0.3,0.3,T
N,any base,0,0,0,0,T
-,gap,0,0,0,0,T"), stringsAsFactors = F)

      alignment <- read.fasta(fasta, seqonly = T)
      alignment <- strsplit(unlist(alignment), split = "")
      alignment <- matrix(unlist(alignment), nrow = length(alignment), ncol = length(alignment[[1]]), byrow = T)
    }

    meep <- c()
    for (i in 1:ncol(alignment)) {
      data <- alignment[, i]

      temp <- match(data, upac$ID)
      colu <- upac[temp, 3:6]

      hey <- colSums(colu)

      # shanon entropy
      p <- hey / sum(hey)
      p[is.na(p)] <- 0
      meep <- rbind(meep, c("ID" = i, p))
    }

    return(meep)
  }


  # Primers in one plot!!
  covcol <- c("Green", "Orange", "Red", "Black")


  # Use the function on all files to make a list of all the base proportions
  # per file


  # Run Vascos entropy function on each input file

  # scoreslist <- lapply(entropy, files)
  scoreslist <- list()
  for (x in 1:length(files)) {
    scoreslist[[x]] <- entropy(files[x])
  }



  # Make a list of the base scores for each base of each file, where the
  # different list items correspond to different input files
  z <- scoreslist

  # Combine the items into a single dataframe where each row has the identifier
  # of the filename
  for (i in 1:length(scoreslist)) {
    z[[i]] <- as.data.frame(scoreslist[[i]])
    z[[i]]$filename <- fastafiles[i]
  }
  all_scores <- do.call(rbind, z)

  colnames(all_scores)[1] <- "nucleotide_position"

  # return(all_scores)


  all_scores <- all_scores %>%
    # Get all the nucleotides into a single column
    gather(A, C, G, `T`, key = "nucleotide", value = "proportion") %>%
    mutate(filename = gsub(".+/| .+|.fasta", "", filename))

  # Make a dataframe without spiders
  no_predators <- all_scores %>%
    filter(!grepl(predator_name, filename))


  predators_only <- all_scores %>%
    filter(grepl(predator_name, filename))

  nosp_summary <- no_predators %>%
    # Use gather by and then summarise to calculate the inter-order mean
    # percentage of that nucleotide at that position
    group_by(nucleotide_position, nucleotide) %>%
    summarise(proportion = mean(proportion)) %>%
    mutate(filename = "all_target_taxa") %>%
    # Reorder the columns so that the two tibbles can be combined
    select(nucleotide_position, filename, nucleotide, proportion)


  min_prey <- 0.75
  max_pred <- 0.25



  # Combine the summary tibble and predators_only tibble
  merged_for_calcs <- left_join(nosp_summary,
    predators_only,
    by = c("nucleotide_position", "nucleotide")
  ) %>%
    # make a new column, for if the minimum prey average threshold AND the maximum
    # predator threshold is met
    mutate(good_base = ifelse(proportion.x >= min_prey &
      proportion.y <= max_pred,
    T, F
    )) %>%
    filter(good_base == T)

  if (output == "list") {
    out_list <- list(merged_for_calcs, nosp_summary, predators_only, all_scores)
    return(out_list)
  }
  if (output == "plot") {


    # More info needs to be added: show the overall consensus bars, and reorder the facet labels
    # Make stacked barplots
    predator_plot <- ggplot(all_scores, aes(x = nucleotide_position, y = proportion)) +
      geom_col(aes(fill = nucleotide)) +
      # Add facets, move facet labels to the left
      facet_wrap(. ~ filename, ncol = 1, strip.position = "left") +
      theme(
        # Rotate the facet labels 180 degrees
        strip.text.y = element_text(angle = 180),
        # remove all the y-axis labelling
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        # format the x axis to increase label size and rotate them by 45 degrees
        axis.text.x = element_text(angle = 45, size = 14)
      ) +
      # Add annotations for bases of interest
      annotate("segment",
        x = merged_for_calcs$nucleotide_position, y = -0.05, xend = merged_for_calcs$nucleotide_position, yend = 0,
        col = "red", arrow = arrow(length = unit(0.3, "cm"))
      ) +
      scale_x_continuous(name = "Nucleotide position", breaks = seq(0, round(max(all_scores$nucleotide_position), 10), 10))


    return(predator_plot)
  }
}
