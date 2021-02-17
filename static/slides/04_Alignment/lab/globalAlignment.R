# https://github.com/BenLangmead/ads1-slides/blob/master/0450_approx__editdist_friends.pdf
# https://www.youtube.com/watch?v=BGV-hUoHF9k&list=PL2mpR0RYFQsBiCWVJSvVAO3OJ2t7DzoHA&index=38
# https://nbviewer.jupyter.org/github/BenLangmead/comp-genomics-class/blob/master/notebooks/CG_DP_Global.ipynb

rm(list = ls())
alphabet <- c("A", "C", "G", "T")
score <- rbind(c(0, 4, 2, 4, 8),
               c(4, 0, 4, 2, 8),
               c(2, 4, 0, 4, 8),
               c(4, 2, 2, 0, 8),
               c(8, 8, 8, 8, 8))


globalAlignment <- function(x, y){
  D = matrix(data = 0, nrow = nchar(x) + 1, ncol = nchar(y) + 1)
  for (i in 2:(nchar(x) + 1)) {
    D[i, 1] <- D[i - 1, 1] + score[which(substr(x, i - 1, i - 1) == alphabet), ncol(score)]
  }
  for (i in 2:(nchar(y) + 1)) {
    D[1, i] <- D[1, i - 1] + score[nrow(score), which(substr(y, i - 1, i - 1) == alphabet)]
  }
  
  for (i in 2:(nchar(x) + 1)) {
    for (j in 2:(nchar(y) + 1)) {
      distVer <- D[i, j - 1] + score[nrow(score), which(substr(y, j - 1, j - 1) == alphabet)]
      distHor <- D[i - 1, j] + score[which(substr(x, i - 1, i - 1) == alphabet), ncol(score)]
      if (substr(x, i - 1, i - 1) == substr(y, j - 1, j - 1)) {
        distDiag <- D[i - 1, j - 1]
      } else {
        distDiag <- D[i - 1, j - 1] + score[which(substr(x, i - 1, i - 1) == alphabet), which(substr(y, j - 1, j - 1) == alphabet)]
      }
      D[i, j] <- min(distVer, distHor, distDiag)
    }
  }
  
  return(D[nrow(D), ncol(D)])
}

x="TACCAGAATTCGA"
y="TACCAGAATTCGA"

globalAlignment(x, y)

x="TACCAGAATTCGA"
y="TACCAGAATTCA"

globalAlignment(x, y)

x="TACCAGAATTCGA"
y="TACCACAATTCGA"

globalAlignment(x, y)

x="TACCAGAATTCGA"
y="TACCAAAATTCGA"

globalAlignment(x, y)
