manage_packages <- function(packages) {
  for (pkg in packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      message(paste("Installing package:", pkg))
      install.packages(pkg, dependencies = TRUE)
    } else {
      message(paste("Package", pkg, "is already installed."))
    }
    library(pkg, character.only = TRUE)
  }
}

# Example usage
packages <- c("tidyverse", "gtsummary", "gt")
manage_packages(packages)

## check out results for funsies
outfile = file.path("/home/nfs/vaithid1/FACETS/RewriteSnpPileup", "outputs", "test1.csv")
output = read.csv(outfile)
output = output[-1,]

## output columns
## chromosome position ref alt  then for each file there is an R E A D

output.long = pivot_longer(output, cols = !c("Chromosome", "Position", "Ref", "Alt"), names_to = c("file", ".value"),  # `file` becomes one column, RAED stay as values
names_pattern = "File(\\d+)([RAED])"  # Regex to separate file number and RAED
) %>%
  mutate(across(
    .cols = -c(Chromosome, Position, Ref, Alt),  # Exclude Ref and Alt columns
    .fns = as.numeric      # Apply as.numeric
  )) %>%
    rename(
      ReferenceCounts = R,
      AlternativeCounts = A,
      Errors = E,
      Deletions = D
    ) 

head(output.long)


write.csv(as.data.frame(summ), file.path("/home/nfs/vaithid1/FACETS/RewriteSnpPileup", "outputs", "test1_pretty.csv"))

  
