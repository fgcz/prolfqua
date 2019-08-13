library(readr)
library(tidyverse)

outdir <- "tmp"
prot <- read_csv("LinearModelP_Values.txt")
pep <- readr::read_csv("PEPTIDE_P_VALUES.txt")


prot <- prot %>% dplyr::rename(fc = estimate)  %>% select(-rhs)
pep <- pep %>% dplyr::rename(fc = estimate)  %>% select(-rhs)

#peptide <- pep %>%filter(lhs == "NR - c") %>% select(-rhs)
pep <- pep %>% tidyr::unite("prot_pep", protein_Id, peptide_Id, remove=FALSE )

file.copy("C:/Users/wolski/prog/LFQService/inst/plotly_reports/ProteinPeptideViewer.Rmd","ProteinPeptideViewer.Rmd")

for(i in unique(pep$lhs)){
  print(i)
  proti <- prot %>%filter(lhs == i)
  pepi <- pep %>% filter(lhs == i)
  name <- make.names(i)
  par <- list(prot = proti, pep=pepi)

  destname<-paste0("ProteinPeptideViewer", name, ".html" )
  rmarkdown::render("ProteinPeptideViewer.Rmd",
                    params=list(par=par),
                    output_file = file.path(outdir, destname))
}
