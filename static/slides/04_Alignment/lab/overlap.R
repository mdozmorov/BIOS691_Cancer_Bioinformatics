# WIP
# https://www.youtube.com/watch?v=iWwrPC768y0&list=PL2mpR0RYFQsBiCWVJSvVAO3OJ2t7DzoHA&index=43


a = "abcdedddcde"
b = "cde"
overlap(a, b)

overlap <- function(a, b, min_length = 3) {
  start <- 1
  while (TRUE) {
    start <- unlist(regexpr(substr(b, start, min_length), a)) # look for b's suffx in a
    if (start == -1) { # no more occurrences to right
      return(0)
    }
    if (unlist(regexpr(substr(a, start, nchar(a)), b)) > 0) {
      return(nchar(a) - start + 1)
    }
    start <- start + 1
  }
}

overlap("TTACGT", "CGTACCGT")

a = "TTACGT"
b = "CGTACCGT"
