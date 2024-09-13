#' Revert antibiotic resistance genes obfuscated by make.names
#' 
#' The built-in function make.names replaces all illegal characters in 
#' R-functions such as parenthesis, apostrophes and so on with dots. This 
#' obfuscates a lot of the gene names in our data set. This function reverts 
#' those names back. Requires a csv containing two columns, one with changed 
#' names labeled Truncated and another with original names labeled Original. 
#' 
#' @param names A character vector. 
#' @return Renamed vector. Strings not found in the list of truncated names
#'  are kept unchanged
#' @examples
#' revert.names("CTX-M-15")
#' revert.names(c("aac.2...IIa", "CTX.M.9", "mismatches are kept"))
revert.names <- function(names) {
  features <- read.csv("original_feature_names.csv")
  
  renamed = c()
  for (name in names) {
    idx <- features$Truncated==name
    if (any(idx)) {
      renamed <- c(renamed, features$Original[idx])
    } else {
      renamed <- c(renamed, name)
    }
  }
  return (renamed)
}

