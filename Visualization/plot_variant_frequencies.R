# plot_variant_frequencies.R
# Description: Plots frequency of variant types or gene categories using ggplot2
# Author: Alexandra Chapman
# Date: March 2025

# Subset the 'variants' data frame to include only the variant types of interest.
# We are only  focusing on certain VEP annotation categories 
Frequencies <- subset(variants, Effect %in% c("intron_variant", "intron_variant,non_coding_transcript_variant", "upstream_gene_variant", "intergenic_variant", "downstream_gene_variant", "non_coding_transcript_exon_variant", "3_prime_UTR_variant", "missense_variant", "synonymous_variant", "5_prime_UTR_variant","intron_variant,NMD_transcript_variant", "stop_gained"))


# Create a bar plot of variant frequencies using ggplot2. with ‘Effect' is on the x-axis
# We set fill' and 'color' to 'Effect' so each variant type has a unique color
ggplot(Frequencies, aes(x = Effect, fill = Effect, color = Effect)) +
  geom_bar()  +
  labs(title = "Distribution of variants by VEP annotation",y = "Frequency", x = "Variant type")+
theme_bw() +
theme(axis.text.x = element_text(angle = 60, hjust = 1, family = "serif", size = 11), 
plot.title = element_text(family = "serif", face = "bold", size = 16), 
axis.title.x = element_text(family = "serif", size = 12, face = "bold"), 
axis.title.y = element_text(family = "serif", size = 10, face = "bold"), 
legend.position = "none", axis.text.y = element_text(family = "serif", size = 9) )

#End of plot_variant_frequencies.R 
