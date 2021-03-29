# http://www.sthda.com/english/wiki/wiki.php?id_contents=7940
install.packages("fpc")
install.packages("dbscan")
# Compute DBSCAN using fpc package
library(fpc)
set.seed(123)
# coords <- c(runif(10, 0, 1), runif(10, 4, 5), runif(10, 6, 7), runif(10, 19, 20)) %>% round

# Small genome
genome_length  <- 100 # Full length of the genome
number_of_TADs <- 10  # Number of TAD boundaries within it
max_TAD_width  <- 10  # Each TAD boundary has certain width
max_TAD_peaks  <- 5   # Within that width, there are peaks with probability 1
# Where points TAD boundaries are located (uniformly across the genome)
TAD_coords     <- round(runif(number_of_TADs, min = 0, max = genome_length))
# Vector to store coordinates of actual TAD boundary regions
coords <- vector(mode = "numeric", length = genome_length)
# Randomly initialize TAD boundary regions
for (i in TAD_coords) {
  # within each region, what are the coordinates with probability 1
  local_TAD_coords <- round(runif(max_TAD_peaks, min = 0, max = max_TAD_width))
  # Add an offset and mark coordinates of TAD boundary peaks with probability 1
  coords[i +  local_TAD_coords] <- min(i + local_TAD_coords, genome_length)
}
# Put things in a data frame and visualize
df <- data.frame(coords, coords)
df <- df[rowSums(df) > 0, ] # Remove coordinates where probabilities of TAD boundaries are <1
plot(df)
# Explore to select eps - inflection point
dbscan::kNNdistplot(df, k =  max_TAD_peaks)
# Need understanding and tweaking
system.time(db <- fpc::dbscan(df, eps = 3, MinPts = 2, method = "raw"))
# Plot DBSCAN results
plot(db, df, main = "DBSCAN", frame = FALSE)
db$cluster %>% table

set.seed(123)
# Medium genome
genome_length  <- 10000 # Full length of the genome
number_of_TADs <- 50  # Number of TAD boundaries within it
max_TAD_width  <- 10  # Each TAD boundary has certain width
max_TAD_peaks  <- 5   # Within that width, there are peaks with probability 1
# Where points TAD boundaries are located (uniformly across the genome)
TAD_coords     <- round(runif(number_of_TADs, min = 0, max = genome_length))
# Vector to store coordinates of actual TAD boundary regions
coords <- vector(mode = "numeric", length = genome_length)
# Randomly initialize TAD boundary regions
for (i in TAD_coords) {
  # within each region, what are the coordinates with probability 1
  local_TAD_coords <- round(runif(max_TAD_peaks, min = 0, max = max_TAD_width))
  # Add an offset and mark coordinates of TAD boundary peaks with probability 1
  coords[i +  local_TAD_coords] <- min(i + local_TAD_coords, genome_length)
}
# Put things in a data frame and visualize
df <- data.frame(coords, coords)
df <- df[rowSums(df) > 0, ] # Remove coordinates where probabilities of TAD boundaries are <1
plot(df)
# Explore to select eps - inflection point
dbscan::kNNdistplot(df, k =  max_TAD_peaks)
# Need understanding and tweaking
system.time(db <- fpc::dbscan(df, eps = 300, MinPts = 2, method = "raw"))
# Plot DBSCAN results
plot(db, df, main = "DBSCAN", frame = FALSE)
db$cluster %>% table


set.seed(123)
# Large genome
genome_length  <- 10000000 # Full length of the genome
number_of_TADs <- 100  # Number of TAD boundaries within it
max_TAD_width  <- 20  # Each TAD boundary has certain width
max_TAD_peaks  <- 10   # Within that width, there are peaks with probability 1
# Where points TAD boundaries are located (uniformly across the genome)
TAD_coords     <- round(runif(number_of_TADs, min = 0, max = genome_length))
# Vector to store coordinates of actual TAD boundary regions
coords <- vector(mode = "numeric", length = genome_length)
# Randomly initialize TAD boundary regions
for (i in TAD_coords) {
  # within each region, what are the coordinates with probability 1
  local_TAD_coords <- round(runif(max_TAD_peaks, min = 0, max = max_TAD_width))
  # Add an offset and mark coordinates of TAD boundary peaks with probability 1
  coords[i +  local_TAD_coords] <- min(i + local_TAD_coords, genome_length)
}
# Put things in a data frame and visualize
df <- data.frame(coords, coords)
df <- df[rowSums(df) > 0, ] # Remove coordinates where probabilities of TAD boundaries are <1
plot(df)
# Explore to select eps - inflection point
dbscan::kNNdistplot(df, k =  max_TAD_peaks)
# Need understanding and tweaking
system.time(db <- fpc::dbscan(df, eps = 100000, MinPts = 2, method = "raw"))
# Plot DBSCAN results
plot(db, df, main = "DBSCAN", frame = FALSE)
db$cluster %>% table
