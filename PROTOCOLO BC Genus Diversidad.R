library(ALDEx2)
library(ampvis2)
library(ape)
library(BiocManager)
library(btools)
library(dada2)
library(DECIPHER)
library(DESeq2)
library(doParallel)
library(plyr)
library(dplyr)
library(fantaxtic)
library(ggplot2)
library(ggpubr)
library(grid)
library(gridExtra)
library(kableExtra)
library(lme4)
library(metagenomeSeq)
library(metagMisc)
library(microbiome)
library(miLineage)
library(openssl)
library(phangorn)
library(philr)
library(PhyloMeasures)
library(phyloseq)
library(plyr)
library(plotly)
library(readxl)
library(readr)
library(reshape2)
library(rlang)
library(scales)
library(tidyr)
library(tidytree)
library(tidyverse)
library(vegan)
library(VennDiagram)
library(xtable)


#Crear objeto Phyloseq-------------------------------------------------------------------------------------------------------------
#importar txt de OTU_table Clasification y Maping  
otu_mat = read.delim("D:/Doctorado en Ciencias Biologicas/00)Publicaciones/4) Articulo solo cianos Bacalar/00) Publicar/R GENUS/TABLA GENUS.txt",row.names=1)
tax_mat = read.delim("D:/Doctorado en Ciencias Biologicas/00)Publicaciones/4) Articulo solo cianos Bacalar/00) Publicar/R GENUS/CLASIFICATION GENUS.txt",row.names=1)
samples_df = read.delim("D:/Doctorado en Ciencias Biologicas/00)Publicaciones/4) Articulo solo cianos Bacalar/00) Publicar/R GENUS/METADATA GENUS.txt",row.names=1)

#Transform into matrixes otu and tax tables (sample table can be left as data frame)
otu_mat <- as.matrix(otu_mat)
tax_mat <- as.matrix(tax_mat)

#Transform to phyloseq objects
OTU = otu_table(otu_mat, taxa_are_rows = TRUE)
TAX = tax_table(tax_mat)
samples = sample_data(samples_df)
rooteBC <- read_tree("D:/Doctorado en Ciencias Biologicas/00)Publicaciones/4) Articulo solo cianos Bacalar/00) Publicar/R GENUS/tree.nwk")
BC <- phyloseq (OTU, TAX, samples, rooteBC)
BC <- prune_samples(sample_names(BC) != "Mock", BC) # Remover potenciales muestras sintéticas
BC

#Barplots-------------------------------------------------------------------------------------------------------------------------------------------------
#Set colors
phylum_colors <-c("purple","red","#009E73","green", "#5E738F", "red4","blue","black","yellow","#619CFF","maroon4","#575329", "#00FECF", "#B05B6F",
                  "#8CD0FF", "#3B9700", "#04F757", "#C8A1A1","#1E6E00","yellow", "#320033", "#66E1D3", "#CFCDAC", "#4FC601", "#4A3B53", "#FF2F80",
                  "#61615A", "#BA0900", "#6B7900", "#00C2A0", "#FFAA92", "#FF90C9", "#B903AA", "#D16100", "#575329", "#00FECF", "#B05B6F", "#8CD0FF", 
                  "#3B9700", "#04F757", "#C8A1A1", "#1E6E00", "#7900D7", "#A77500", "#6367A9", "#A05837", "#6B002C", "#772600", "#D790FF", "#9B9700",
                  "#3B5DFF", "#549E79","#7B4F4B", "#A1C299", "purple", "#999999", "#E69F00",  "#009E73","darkorange1","darkgrey","#FC4E07","#D14285",  
                  "#652926", "red4" , "#CC79A7",  "#009E73","#00BA38" , "#252A52" , "brown",    "#D55E00" , "cyan1", "royalblue4","#CBD588", "#5F7FC7", 
                  "#CC79A7","orange","#56B4E9","#DA5724", "#508578","#C84248", "darkorchid", "#AD6F3B", "#673770","#D14285", "#652926", "#C84248", "#8569D5",
                  "darkseagreen", "#5E738F","#D1A33D","#8A7C64", "#599861", "yellow","darkgoldenrod1", "#56B4E9","darkolivegreen1","#F0E442", "#0072B2", "#D55E00")

#Plot genus
P<- ggplot(BC, aes(x = Site, y = Abundance, fill = Genus)) +
  geom_bar(stat = "identity",position = "Fill") + 
   scale_fill_manual(values = phylum_colors) +
  theme(axis.title.x = element_blank()) +
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) +
  ylab("Relative Abundance") +
  ggtitle("Genus")
P <- P+ theme_classic()
P

#Barplot por abundancia-------------------------------------------------------------------------------------------------------------
par(mar = c(10, 4, 4, 2) + 0.1) # make more room on bottom margin
N <- 30
barplot(sort(taxa_sums(BC), TRUE)[1:N]/nsamples(BC), las=2)

#Plot tree por forma y abundancia----------------------------------------------------------------------------------------------------
plot_tree(BC, color="Biome", shape="Class", label.tips="Genus", size="Abundance")

#plot_redes de coocurrecnia----------------------------------------------------------------------------------------------------------
plot_net(BC, maxdist=0.4, color="Site", shape="Biome")

#Alpha diversidad--------------------------------------------------------------------------------------------------------------------
#Cuidado Table deben ser numeros enteros no decimales
alpha_meas = c("Observed", "Chao1", "Shannon", "Simpson")
(p <- plot_richness(BC, "Site", color = "Biome", measures=alpha_meas))

#Phylomesures PD Diversidad Filogenetica----------------------------------------------------------------------------------------------
# Estimate phylogenetic diversity (PD) and mean pairwise phylogenetic distance (MPD) Donde: MNTD = NTI; MPD = NRI 
PD <- phyloseq_phylo_div(BC, measures = c("PD", "MPD", "SES.PD", "SES.MPD", "SES.MNTD"))
PD
write.csv(PD, "Tabla_PD-BIEN.csv")

# Standardized effect size of phylogenetic diversity (PDI)
phyloseq_phylo_div(BC, measures = "PD", standardize = T)
phyloseq_phylo_div(BC, measures = "SES.PD")  # the same



#Ampvis (Heatmap y PCA)--------------------------------------------------------------------------------------------------------------
#Convertir de phyloseq a ampvis
Totu_table =t(otu_table(BC))
otu_table(BC)=Totu_table
av2_otutable <- data.frame(OTU = rownames(t(phyloseq::otu_table(BC)@.Data)),
                           t(phyloseq::otu_table(BC)@.Data),
                           phyloseq::tax_table(BC)@.Data,
                           check.names = F)
#Extract metadata from the phyloseq object:
av2_metadata <- data.frame(phyloseq::sample_data(BC), 
                           check.names = F)
av2_metadata <- cbind(rownames(av2_metadata), av2_metadata)
#Load the data with amp_load:
av2_obj <- amp_load(av2_otutable, av2_metadata)
amp_rankabundance(av2_obj, group_by = "Site")
sqrt
log10

#Heatmap Top--------------------------------------------------------------------------
amp_heatmap(av2_obj,
            group_by = "Site",
            facet_by = "Biome",
            tax_aggregate = "Genus",
            tax_show = 44,
            color_vector = c("white", "red"),
            plot_colorscale = "sqrt",
            plot_values = TRUE) +
  theme(axis.text.x = element_text(angle = 45, size=12, vjust = 1),
        axis.text.y = element_text(size=12),
        legend.position="left")
#PCA--------------------------------------------------------------------------------
amp_ordinate(
  data = av2_obj,
  type = "pca",
    sample_color_by = "Biome",
  sample_colorframe = FALSE,
  sample_colorframe_label = "Site",
  species_plotly = TRUE)             
       
#CA--------------------------------------------------------------------------------
amp_ordinate(
  data = av2_obj,
  type = "ca",
  sample_color_by = "Biome",
  sample_colorframe = FALSE,
  sample_colorframe_label = "Site",
  species_plotly = TRUE)     

#Diagrama de Venn--------------------------------------------------------------------------------
amp_venn(av2_obj, group_by = "Biome", cut_a = 0.001, cut_f = 1, text_size = 5)   



#NMDS----------------------------------------------------------------------------------------------
tab = read.delim("D:/Doctorado en Ciencias Biologicas/00)Publicaciones/4) Articulo solo cianos Bacalar/00) Publicar/R GENUS/TABLA GENUS.txt",row.names=1)
env = read.delim("D:/Doctorado en Ciencias Biologicas/00)Publicaciones/4) Articulo solo cianos Bacalar/00) Publicar/R GENUS/METADATA GENUS.txt",row.names=1)
tab
env
#definir columnas tabla y metadata
com = tab[,1:19]
com
envi = env[,1:10]
envi
#convert com to a matrix
m_com = as.matrix(com)
##Perform the NMDS ordination,nmds code
set.seed(123)
nmds = metaMDS(m_com, distance = "bray")
nmds
#Graficar
plot(nmds, type = "t") 
plot(envi)



