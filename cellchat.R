library(CellChat)
table(penny$group)
DefaultAssay(penny) <- 'RNA'
penny_BA <- subset(penny, subset = group == 'BA')
penny_Naive <- subset(penny, subset = group == 'Naive')
penny_Vehicle <- subset(penny, subset = group == 'Vehicle')

cellchat <- createCellChat(object = penny_BA, group.by = "EK_anno")
CellChatDB <- CellChatDB.mouse
CellChatDB.use <- CellChatDB
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.mouse)
cellchat <- computeCommunProb(cellchat, type = "truncatedMean", trim = 0.2, raw.use = FALSE, population.size = TRUE)
cellchat <- filterCommunication(cellchat, min.cells = 3)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
cellchat_BA <- cellchat

cellchat <- createCellChat(object = penny_Naive, group.by = "EK_anno")
CellChatDB <- CellChatDB.mouse
CellChatDB.use <- CellChatDB
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.mouse)
cellchat <- computeCommunProb(cellchat, type = "truncatedMean", trim = 0.2, raw.use = FALSE, population.size = TRUE)
cellchat <- filterCommunication(cellchat, min.cells = 3)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
cellchat_Naive <- cellchat

penny_Vehicle@active.ident <- droplevels(penny_Vehicle@active.ident, exclude = setdiff(levels(penny_Vehicle@active.ident), unique(penny_Vehicle@active.ident)))
penny_Vehicle$EK_anno <- droplevels(penny_Vehicle$EK_anno, exclude = setdiff(levels(penny_Vehicle$EK_anno), unique(penny_Vehicle$EK_anno)))

cellchat <- createCellChat(object = penny_Vehicle, group.by = "EK_anno")
CellChatDB <- CellChatDB.mouse
CellChatDB.use <- CellChatDB
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.mouse)
cellchat <- computeCommunProb(cellchat, type = "truncatedMean", trim = 0.2, raw.use = FALSE, population.size = TRUE)
cellchat <- filterCommunication(cellchat, min.cells = 3)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
cellchat_Vehicle <- cellchat


ht1 <- netAnalysis_signalingRole_heatmap(cellchat_Vehicle, pattern = "outgoing", height = 28)
ht2 <- netAnalysis_signalingRole_heatmap(cellchat_Vehicle, pattern = "incoming", height = 28)
ht1 + ht2
