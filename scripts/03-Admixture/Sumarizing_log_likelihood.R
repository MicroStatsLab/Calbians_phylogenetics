# K values
K <- 6:30

# file names
files <- paste0("data_out/03-Admixture/901_isolates/260109full_901_outfiles_K", K, ".log")

# function to extract best like
get_best_like <- function(f) {
  line <- grep("^best like=", readLines(f), value = TRUE)
  as.numeric(sub("best like=([^ ]+).*", "\\1", line))
}

# extract values
best_like <- sapply(files, get_best_like)

# make table
best_like_table <- data.frame(
  K = K,
  best_like = best_like
)

best_like_table


best_like_table_clean <- best_like_table

# Remove the path, keep only the file name
rownames(best_like_table_clean) <- basename(rownames(best_like_table_clean))

best_like_table_clean

rownames(best_like_table_clean) <- sub("\\.log$", "", basename(rownames(best_like_table_clean)))

best_like_table_clean

best_like_table_clean <- best_like_table_clean
rownames(best_like_table_clean) <- paste0("K", best_like_table_clean$K)
best_like_table_clean$K <- NULL

best_like_table_clean

options(scipen = 999)


plot(
  best_like_table$K,
  best_like_table$best_like,
  type = "b",
  pch = 19,
  xlab = "K",
  ylab = "Log-likelihood",
  main = ""
)

library(ggplot2)
library(scales)  # for nice numeric formatting

ggplot(best_like_table, aes(x = K, y = best_like)) +
  geom_line(color = "black", linewidth = 1) +
  geom_point(size = 3, color = "black") +
  scale_y_continuous(labels = label_number(big.mark = ",")) +
  labs(
    x = "K",
    y = "Log-likelihood",
    title = ""
  ) +
  theme_minimal(base_size = 14)

