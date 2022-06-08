# fix example dataset
set.seed(76)
datadir <- file.path(find.package("prolfquadata") , "quantdata")
inputMQfile <-  file.path(datadir, "MAXQuant_IonStar2018_PXD003881.zip")

tmp <- readLines(unz(inputMQfile,"proteinGroups.txt"))
tmpP <- read.csv(unz(inputMQfile,"proteinGroups.txt"),
                 header = TRUE, sep = "\t",
                 stringsAsFactors = FALSE)

proteins <- sample(tmpP$Protein.IDs,100)
length(unique(proteins))

tmptiny <- c(tmp[1], tmp[2:length(tmp)][tmpP$Protein.IDs %in% proteins])
writeLines(tmptiny, "proteinGroups.txt")

# create tiny peptides
tmp <- readLines(unz(inputMQfile,"peptides.txt"))
tmp2 <- read.csv(unz(inputMQfile,"peptides.txt"),
                 header = TRUE, sep = "\t",
                 stringsAsFactors = FALSE)

sum(unique(tmp2$Proteins) %in% proteins)

tinypep <- c(tmp[1], tmp[2:length(tmp)][tmp2$Proteins %in% proteins])
writeLines(tinypep, "peptides.txt")

# create tiny evidence
tmp <- readLines(unz(inputMQfile,"evidence.txt"))
tmp2 <- read.csv(unz(inputMQfile,"evidence.txt"),
                 header = TRUE, sep = "\t",
                 stringsAsFactors = FALSE)
sum(unique(tmp2$Proteins) %in% proteins)

tinyevi <- c(tmp[1],tmp[2:length(tmp)][tmp2$Proteins %in% proteins])
writeLines(tinyevi, "evidence.txt")

zip("tiny2.zip", c("evidence.txt","peptides.txt","proteinGroups.txt"))
file.remove(c("evidence.txt","peptides.txt","proteinGroups.txt"))
