# Data cleaning functions
# John Godlee (johngodlee@gmail.com)
# 2019_12_05

set.seed(364066)

# Colour palettes
degrad_pal <- c("#369FC2", "#DB2833")
big_pal <- RColorBrewer::brewer.pal(4, "Dark2")

# Create a stem/tree abundance matrix, with plotcode id column
stem_ab <- function(data, id, sp){
  id <- enquo(id)
  sp <- enquo(sp)
  stem_ab_mat <- data %>%
    filter(!is.na(!!sp)) %>%
    group_by(!!id, !!sp, .drop = FALSE) %>%
     tally() %>%
    spread(as_label(sp), n, fill = 0) %>%
    ungroup() %>%
    rename("plotcode" = !!names(.[1]))
  return(stem_ab_mat)
}

# Extract taxa from taxize classification
taxon_extract <- function(x, taxon) {
  ifelse(is.na(x), return(NA), return(x[x$rank == taxon,][1]))
}

# Create tidy dataframe from taxize classification
taxa_df_gen <- function(x) {
  taxa_list <- lapply(x, function(y) {
    data.frame(family = c(as.character(taxon_extract(y, taxon = "family"))),
      species = as.character(c(taxon_extract(y, taxon = "species"))))
  })
  taxa_df <- do.call(rbind, taxa_list)
  rownames(taxa_df) <- NULL
  taxa_df$family <- as.character(taxa_df$family)
  taxa_df$species <- as.character(taxa_df$species)
  taxa_df <- na.omit(taxa_df)
  taxa_df <- taxa_df[!taxa_df$species == "character(0)",]
  return(taxa_df)
}

# Run NMDS with different dimensions 
NMDS.scree <- function(x) {
  plot(rep(1, 10), replicate(10, metaMDS(x, autotransform = F, k = 1)$stress), xlim = c(1, 10),ylim = c(0, 0.30), xlab = "# of Dimensions", ylab = "Stress", main = "NMDS stress plot")
  for (i in 1:10) {
    points(rep(i + 1,10),replicate(10, metaMDS(x, autotransform = F, k = i + 1)$stress))
  }
}

# Format p values for text
p_format <- function(p){
  dplyr::case_when(p < 0.01 ~ "p<0.01",
    p < 0.05 ~ "p<0.05",
    TRUE ~ paste0("p = ", as.character(round(p, digits = 2))))
}

# POM adjustment function

##' @title Stem diameter Point Of Measurement (POM) adjustment
##' @description  Function to estimate stem diameter at 1.3 given measurements 
##'   at other POMs. This is an important adjustment if the stem's biomass is to
##'   be estimated using an allometric eq, which uses diameter at 1.3 m height
##'   as a predictor, e.g. Chave (2014)
##' @author Casey M Ryan
##' @return d130, the estimated diameter at a POM of 1.3 m (in cm). NB: This might
##'   be negative for small trees measured at POMs <1.3m.
##'   
##' @param d_in the diameter measured at the POM (in cm)
##' @param POM the height of the POM (in m)
##'   
##' @details The adjustment is based on a tree taper model developed as part of 
##'   ACES, using data from the miombo of Niassa. The model is a cubic 
##'   polynomial, with three equations for different sized stems. It is far from
##'   perfect! The three size classes are defined by the edges variable, which is
##'   hard coded
##' @section Warning: The model should not be used for POMs above 1.7 m. Extrapolating 
##'   beyond the data will give nonsense. Thus, POMs >1.7 m are not adjusted.
##' @examples
##' POMadj(10, 0.3)
##' POMadj(1, 0.3) # d130 is negative, i.e. the stem probably wasn't 1.3 m tall
##' POMadj(50, 1.9) #Â generates warning, as outside calibration data range
##' \dontrun{
##' POMadj(50, 0) #Â zero or -ve POM is outside range, or nonsense
##' }
##' TODO: not happy that there is still a small correction even when the POM is 130
POMadj <- function(d_in, POM) {
  stopifnot(is.numeric(d_in),
    is.numeric(POM),
    POM >= 0,
    sum(is.na(POM))==0,
    length(POM) == length(d_in))
  if (any(POM > 1.7))
    warning("POMs >1.7 m are outside the calibration data and result in extrapolated nonsense, no correction applied")
  
  NAS <- is.na(d_in)
  d_in_clean <- d_in[!NAS]
  POM_clean <- POM[!NAS]
  # define the size class edges:
  edges <- c(5.0, 15.8, 26.6, 37.4)
  sm <- d_in_clean < edges[2]
  med <- d_in_clean >= edges[2] & d_in_clean < edges[3]
  lg <- d_in_clean >= edges[3]
  
  # compute apredictions for delta_d, for all size classes
  delta_d <- data.frame(
    # if small:
    small =  3.4678+-5.2428 * POM_clean  + 2.9401 * POM_clean ^ 2+-0.7141 * POM_clean ^3,
    # if med
    med =  4.918+-8.819 * POM_clean  + 6.367 * POM_clean ^ 2+-1.871 * POM_clean ^3,
    # if large
    large =  9.474+-18.257 * POM_clean  + 12.873 * POM_clean ^ 2+-3.325 * POM_clean ^3
  )
  # index into the right size class
  dd <- NA_real_
  dd[sm] <- delta_d$small[sm]
  dd[med] <- delta_d$med[med]
  dd[lg] <- delta_d$large[lg]
  dd[POM_clean > 1.7] <- 0 #Â to avoid extrapolation mess
  
  #Â add back in NAs
  d130 <- NA
  d130[NAS] <- NA
  d130[!NAS] <- d_in_clean - dd
  
  
  if (any(d130[!NAS] < 0))
    warning("Negative d130 estimated, repaced with NA")
  d130[d130<=0 & !is.na(d130)] <- NA
  return(d130)
}

# Define a function to draw a convex hull 
find_hull <- function(df){
  df[chull(df$PCoA1, df$PCoA2), ]}

# Format output of linear model
lm_format <- function(x){
  paste0("F(", 
  summary(x)$fstatistic[2],
  ",",
  summary(x)$fstatistic[3],
  ") = ",
  sprintf("%.2f", round(summary(x)$fstatistic[1], 2)),
  ", ",
  p_format(anova(x)$`Pr(>F)`[1]))
}



