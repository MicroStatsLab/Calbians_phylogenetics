# --- Libraries ---
library(tidyverse)
library(Polychrome)
library(colorspace)
library(pophelper)

# --- Ks to plot ---
K_values <- c(2, 3, 4, 6, 9, 15, 27)
K_values <- c(2, 5,10, 16, 27)
# --- File paths ---
qfiles <- sprintf("/Users/abdurahman/Nextcloud/19Abdul-Rahman/C_albicans_manuscript/241101_Current/data_out/03-Admixture/outfiles_K%d.qopt", K_values)

# --- Read Q-matrices into a qlist ---
qlist <- readQ(qfiles)


# --- Align clusters across Ks to fix label switching ---
qlist_aligned <- alignK(qlist, type = "across")
plotQ(qlist_aligned,exportpath="/Users/abdurahman/Nextcloud/19Abdul-Rahman/C_albicans_manuscript/241101_Current/data_out/03-Admixture")
