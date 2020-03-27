
makeClass <-  function(x){

  .ncol <- function(){
    .x = x
    return(ncol(.x))
  }
  .nrow <- function(){
    .x = x
    return(nrow(.x))
  }
  .cells <- function(){
    .x = x
    return(nrow(.x)*ncol(.x))
  }
  return(list(nrow= .nrow, ncol = .ncol, cells = .cells))
}


datF <- data.frame(r = 1:3, b = letters[1:3])
clss <- makeClass(datF)

#print(datF)
#nrow(datF)

clss$nrow()
clss$ncol()
clss$cells()
