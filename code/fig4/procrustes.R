
if (!require("vegan")) install.packages("vegan")
library(vegan)
library(readxl)
metabolite_data = metab
microbe_data = merged_data
common_samples <- intersect(rownames(metabolite_data), rownames(microbe_data)) 
metabolite_data <- metabolite_data[common_samples, ]
microbe_data <- microbe_data[common_samples, ]
all(rownames(metabolite_data) == rownames(microbe_data))

metabolite_data[] <- lapply(metabolite_data, function(x) {
  if (is.character(x)) {
    as.numeric(x) 
  } else {
    x 
  }
})
microbe_data[] <- lapply(microbe_data, function(x) {
  if (is.character(x)) {
    as.numeric(x)
  } else {
    x 
  }
})
#---------------------------------------------------------------------------------

group = read.delim('./data/metadata.txt')
group$Timepoint0 <- ifelse(group$Group == "W/ Infection", 
                           paste0("W_", group$Timepoint), 
                           ifelse(group$Group == "W/O Infection", 
                                  paste0("WO_", group$Timepoint), 
                                  group$Timepoint))

wt0 = group$SampleID
merged_data_target = microbe_data[rownames(microbe_data)%in% wt0,]
metab_target = metabolite_data[rownames(metabolite_data)%in% wt0,]


#---------------

metabolite_data = metab_target

microbe_data = merged_data_target

#-----------------------------------------------------------------------------------

microbe_data <- microbe_data[, !(apply(microbe_data, 2, function(x) var(x) == 0))]
#--------------------------------------------------
# metabolite_data = scale(metabolite_data)
# microbe_data = scale(microbe_data)
metabolite_pca <- prcomp(metabolite_data,center = T,scale=T)
metabolite_scores <- metabolite_pca$x  


microbe_pca <- prcomp(microbe_data,center = T,scale=T)
microbe_scores <- microbe_pca$x 

metabolite_pca_reduced <- metabolite_scores[, 1:2]
microbe_pca_reduced <- microbe_scores[, 1:2]

pro <- procrustes(metabolite_pca_reduced, microbe_pca_reduced)

summary(pro)  
#------------------------------------------
x_coords <- pro$X[, 1] 
y_coords <- pro$X[, 2] 


outlier_indices <- which(x_coords >10 | y_coords > 20)  
print(outlier_indices)
print(rownames(pro$X)[outlier_indices])

metabolite_data_clean <- metabolite_data[-outlier_indices, ]
microbe_data_clean <- microbe_data[-outlier_indices, ]

pro_clean <- pro
pro_clean$X <- pro$X[-outlier_indices, ]
pro_clean$Yrot <- pro$Yrot[-outlier_indices, ]
save(pro_clean,mantel_result,file= './figure/suppdata/fig4/procrustes.RData')
plot(pro_clean$X[, 1], pro_clean$X[, 2], 
     xlim = c(-5, 10),
     xlab = "PC 1", ylab = "PC 2",
     pch = 24, col = rgb(0.6118, 0.6902, 0.7647), bg = rgb(0.6118, 0.6902, 0.7647,0.2), cex = 3, lwd = 2)

points(pro_clean$Yrot[, 1], pro_clean$Yrot[, 2], 
       pch = 23, col = rgb(0.956, 0.874, 0.827),cex = 3, lwd = 2)

for (i in 1:nrow(pro_clean$X)) {
  lines(c(pro_clean$X[i, 1], pro_clean$Yrot[i, 1]), 
        c(pro_clean$X[i, 2], pro_clean$Yrot[i, 2]), 
        col = rgb(0.5, 0.5, 0.5, 0.3), lty = 1, lwd = 2)
}

text(0.065, -0.09, expression(italic("P value = 0.006")), 
     cex = 1.2, col = "black", pos = 4) 
pro_test <- protest(metabolite_pca_reduced, microbe_pca_reduced)

summary(pro_test)
print(pro_test$signif)
mantel_result <- mantel(vegdist(metabolite_data), vegdist(microbe_data), method = "spearman")
print(mantel_result)

# library(Rtsne)
# tsne_result <- Rtsne(microbe_data)
# tsne_result1 <- Rtsne(metabolite_data)
# plot(tsne_result$Y, 
#      xlim = c(-10, 10), ylim = c(-10, 10),
#      col = "purple", pch = 17)
# points(tsne_result1$Y[, 1], tsne_result1$Y[, 2], 
#        pch = 1,         
#        # 三角形点
#        col = rgb(0.5, 0, 0.5, 0.5),           # 紫色三角形点，带透明度
#        cex = 1)                             # 调整点大小
# for (i in 1:nrow(tsne_result$Y)) {
#   lines(c(tsne_result$Y[i, 1], tsne_result1$Y[i, 1]), 
#         c(tsne_result$Y[i, 2], tsne_result1$Y[i, 2]), 
#         col = "grey", lty = 1, lwd = 1)     # 灰色细线
# }
# #------------------------------------------------------------------------------
# 
# 
# 
# metabolite_dist <- vegdist(metabolite_data, method = "bray")
# 
# # 确保菌群数据标准化（如 Hellinger 转换）
# microbe_data <- decostand(microbe_data, method = "hellinger")
