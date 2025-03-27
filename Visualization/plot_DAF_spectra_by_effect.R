# plot_DAF_spectra_by_effect.R
# Description: Subset the 'variants' data frame by effect type, then plot DAF histograms using ggplot2.
# Author: Alexandra Chapman
# Date: March 2025

library(ggplot2)
library(dplyr)

# 1. Subset  main dataset by effect type
synonymous <- subset(variants, Effect %in% c("synonymous_variant"))
missense   <- subset(variants, Effect %in% c("missense_variant"))
stop_gained <- subset(variants, Effect %in% c("stop_gained"))
# ... etc. for intron, upstream, etc.

# 2. Create histograms
missense_DAF_spectra <- ggplot(missense, aes(x=DAF)) +
  geom_histogram(binwidth = 0.05, fill = "pink", color = "black") +
  scale_x_log10() +
  labs(title = "DAF Spectrum for missense variants", x = "Derived Allele Frequency", y = "Count") +
  theme_bw() +
  theme(
    plot.title = element_text(family = "serif", face = "bold", size = 18),
    axis.title.x = element_text(family = "serif", size = 14),
    axis.title.y = element_text(family = "serif", size = 14),
    axis.text.x  = element_text(size = 12),
    axis.text.y  = element_text(size = 12)
  )
print(missense_DAF_spectra)

# The same steps for synonymous, stop_gained, etc.
# ...
# End of plot_DAF_spectra_by_effect.R
