library(SKAT)
library(ggplot2)

data(SKAT.haplotypes)
names(SKAT.haplotypes)
attach(SKAT.haplotypes)  

set.seed(500)
out.c <- Power_Continuous_R(Haplotype, 
                            SNPInfo$CHROM_POS, 
                            SubRegion.Length = 5000, 
                            Causal.Percent = 20, 
                            N.Sim = 100, 
                            MaxBeta = 2,
                            Negative.Percent = 20, 
                            r.corr = 2)

x <- as.vector(rownames(out.c$Power))
y <- as.vector(colnames(out.c$Power))
data <- expand.grid(x=x, y=y)
data$Z <- as.vector(out.c$Power)

# Heatmap 
ggplot(data, aes(x, y, fill= Z)) + 
  geom_tile() +
  scale_fill_gradient(low = "lightblue", high = "Red") +
  theme_bw()

Get_RequiredSampleSize(out.c, Power=0.6)
out.c
