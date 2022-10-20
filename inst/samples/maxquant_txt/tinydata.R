# fix example dataset
set.seed(77)
datadir <- file.path(find.package("prolfquadata") , "quantdata")
inputMQfile <-  file.path(datadir, "MAXQuant_IonStar2018_PXD003881.zip")

tmp <- readLines(unz(inputMQfile,"proteinGroups.txt"))
tmpP <- read.csv(unz(inputMQfile,"proteinGroups.txt"),
                 header = TRUE, sep = "\t",
                 stringsAsFactors = FALSE)

length(tmpP$Protein.IDs)
tmpP <- tmpP[!grepl("CON__|REV__", tmpP$Protein.IDs),]

proteins <- sample(tmpP$id, 400)
length(unique(proteins))

tmptiny <- c(tmp[1], tmp[2:length(tmp)][proteins])
writeLines(tmptiny, "proteinGroups.txt")

# create tiny peptides
tmp <- readLines(unz(inputMQfile,"peptides.txt"))
tmp2 <- read.csv(unz(inputMQfile,"peptides.txt"),
                 header = TRUE, sep = "\t",
                 stringsAsFactors = FALSE)

sum(tmp2$Protein.group.IDs %in% proteins)

tinypep <- c(tmp[1], tmp[2:length(tmp)][tmp2$Protein.group.IDs %in% proteins])
writeLines(tinypep, "peptides.txt")

# create tiny evidence
tmp <- readLines(unz(inputMQfile,"evidence.txt"))
tmp2 <- read.csv(unz(inputMQfile,"evidence.txt"),
                 header = TRUE, sep = "\t",
                 stringsAsFactors = FALSE)
sum(unique(tmp2$Protein.group.IDs) %in% proteins)

tinyevi <- c(tmp[1],tmp[2:length(tmp)][tmp2$Protein.group.IDs %in% proteins])
writeLines(tinyevi, "evidence.txt")

zip("tiny2.zip", c("evidence.txt","peptides.txt","proteinGroups.txt"))
file.copy("tiny2.zip",to = "inst/samples/maxquant_txt/tiny2.zip", overwrite = TRUE)

file.remove(c("tiny2.zip","evidence.txt","peptides.txt","proteinGroups.txt"))
