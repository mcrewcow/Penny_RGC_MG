penny <- LoadH5Seurat('C://Bioinf/Penny_project/combined_EKanno.h5Seurat')
library(escape)
gene.sets1 <- getGeneSets(library = "C5", gene.sets = c('GOBP_REGULATION_OF_ANTIGEN_PROCESSING_AND_PRESENTATION',
                                                        'GOBP_REGULATION_OF_ANTIGEN_RECEPTOR_MEDIATED_SIGNALING_PATHWAY',
                                                        'GOBP_REGULATION_OF_APOPTOTIC_CELL_CLEARANCE',
                                                        'GOBP_REGULATION_OF_APOPTOTIC_DNA_FRAGMENTATION',
                                                        'GOBP_REGULATION_OF_APOPTOTIC_PROCESS_INVOLVED_IN_DEVELOPMENT',
                                                        'GOBP_REGULATION_OF_APOPTOTIC_SIGNALING_PATHWAY',
                                                        'GOBP_REGULATION_OF_ATP_BIOSYNTHETIC_PROCESS',
                                                        'GOBP_REGULATION_OF_ATP_DEPENDENT_ACTIVITY',
                                                        'GOBP_REGULATION_OF_ATP_METABOLIC_PROCESS',
                                                        'GOBP_REGULATION_OF_AUTOPHAGY',
                                                        'GOBP_REGULATION_OF_AUTOPHAGY_OF_MITOCHONDRION',
                                                        'GOBP_REGULATION_OF_AUTOPHAGY_OF_MITOCHONDRION_IN_RESPONSE_TO_MITOCHONDRIAL_DEPOLARIZATION',
                                                        'GOBP_ACTIVATION_OF_IMMUNE_RESPONSE',
                                                        'GOBP_ACTIVATION_OF_INNATE_IMMUNE_RESPONSE',
                                                        'GOBP_ACTIVATION_OF_NF_KAPPAB_INDUCING_KINASE_ACTIVITY',
                                                        'GOBP_ACUTE_INFLAMMATORY_RESPONSE',
                                                        'GOBP_REGULATION_OF_BLOOD_BRAIN_BARRIER_PERMEABILITY',
                                                        'GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION',
                                                        'GOBP_REGULATION_OF_CELLULAR_RESPONSE_TO_OXIDATIVE_STRESS',
                                                        'GOBP_REGULATION_OF_CELLULAR_RESPONSE_TO_STRESS',
                                                        'GOBP_ATP_BIOSYNTHETIC_PROCESS',
                                                        'GOBP_ATP_METABOLIC_PROCESS',
                                                        'GOBP_REGULATION_OF_JNK_CASCADE',
                                                        'GOBP_REGULATION_OF_JUN_KINASE_ACTIVITY',
                                                        'GOBP_REGULATION_OF_MACROPHAGE_APOPTOTIC_PROCESS',
                                                        'GOBP_REGULATION_OF_MICROGLIAL_CELL_ACTIVATION',
                                                        'GOBP_REGULATION_OF_NEURON_APOPTOTIC_PROCESS',
                                                        'GOBP_COMPLEMENT_ACTIVATION',
                                                        'GOBP_CYTOKINE_MEDIATED_SIGNALING_PATHWAY',
                                                        'GOBP_CYTOKINE_PRODUCTION',
                                                        'GOBP_CYTOKINE_PRODUCTION_INVOLVED_IN_IMMUNE_RESPONSE',
                                                        'GOBP_CYTOKINE_PRODUCTION_INVOLVED_IN_INFLAMMATORY_RESPONSE',
                                                        'GOBP_CYTOKINESIS',
                                                        'GOBP_CYTOKINETIC_PROCESS',
                                                        'GOBP_REGULATION_OF_PHAGOCYTOSIS',
                                                        'GOBP_ENDOTHELIAL_CELL_ACTIVATION',
                                                        'GOBP_ENDOTHELIAL_CELL_APOPTOTIC_PROCESS',
                                                        'GOBP_EXTRACELLULAR_MATRIX_CELL_SIGNALING',
                                                        'GOBP_EXTRACELLULAR_MATRIX_CONSTITUENT_SECRETION',
                                                        'GOBP_GLIAL_CELL_ACTIVATION',
                                                        'GOBP_GLIAL_CELL_APOPTOTIC_PROCESS',
                                                        'GOBP_GLIAL_CELL_DERIVED_NEUROTROPHIC_FACTOR_RECEPTOR_SIGNALING_PATHWAY',
                                                        'GOBP_INFLAMMATORY_RESPONSE',
                                                        'GOBP_MACROPHAGE_ACTIVATION',
                                                        'GOBP_MICROGLIA_DIFFERENTIATION',
                                                        'GOBP_MICROGLIAL_CELL_ACTIVATION_INVOLVED_IN_IMMUNE_RESPONSE'),
                                                        species = 'Mus musculus')
ES <- enrichIt(obj = penny,
gene.sets = gene.sets1,
groups = 3000)
View(gene.sets1)
penny <- AddMetaData(penny, ES)
penny <- SetIdent(penny, value = 'EK_anno')
ES2 <- data.frame(penny[[]], Idents(penny))
colnames(ES2)[ncol(ES2)] <- "cluster"
head(ES2)
ridgeEnrichment(ES2, gene.set = "GOBP_ACTIVATION_OF_IMMUNE_RESPONSE", group = 'EK_anno',  add.rug = TRUE)

ES2 %>%
ggplot( aes(x=GOBP_ACTIVATION_OF_IMMUNE_RESPONSE, fill=group)) +
geom_histogram( color="#e9ecef", alpha=0.6, position = 'identity', bins = 30) +
scale_fill_manual(values=c("#69b3a2", "#404080", '#a32000')) + facet_wrap(~EK_anno) + 
labs(fill="")

library(ggplot2)

# Ensure ES2 has a 'group' column and columns for each pathway in gene.sets1
# Iterate through the names of the pathways in gene.sets1
for(pathway_name in names(gene.sets1)) {
  
  # Construct the plot
  plot <- ggplot(ES2, aes_string(x=pathway_name, fill="group")) +  # Use aes_string to refer to the column by name
    geom_histogram(color="#e9ecef", alpha=0.6, position='identity', bins=30) +
    scale_fill_manual(values=c("#69b3a2", "#404080", '#a32000')) +
    facet_wrap(~EK_anno) +
    labs(fill="") +
    ggtitle(pathway_name)  # Use the pathway name as the plot title
  
  # Define the PDF file name based on the pathway name
  file_name <- paste0("C:/Users/rodri/Downloads/penny_gsea/", gsub("[^a-zA-Z0-9]", "_", pathway_name), ".pdf")  # Sanitize the pathway name to create a valid file name
  
  # Save the plot as a PDF, 8x8 inches
  ggsave(file_name, plot, width=8, height=4)
}
