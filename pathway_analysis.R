
############################################
#### DEG in NF-kappB signaling pathway ####
############################################

library(Seurat)
library(KEGGREST)
library(ComplexHeatmap)

# identifying DEG 
DEG<-FindMarkers(dat.seurat,ident.1="CD34+ SC",ident.2="Bst1-high",logfc.threshold = 0.5)
DEG<-DEG[DEG$p_val_adj<0.01,]

# mapping DEG to KEGG pathway 
NFkappaB.KEGG.names <- keggGet("mmu04064")[[1]]$GENE
gene.symbols.NFkappaB<-intersect(gene.symbols.NFkappaB,row.names(DEG))

# taking average gene expression
avg.expr=sapply(levels(CC.label.cell.sort),function(x){
  cell.select=names(CC.label.cell.sort[CC.label.cell.sort==x])
  cluster.avg=apply(data.norm[,cell.select],1,mean)
})
expr.NFkappaB.avg<-as.data.frame(avg.expr[gene.symbols.NFkappaB,])

# heatmap
heatmap_afterTrajectory4avgExpr<-function(expr,heatmap_title,heatmap_col,file){
  ht<-Heatmap(expr,column_title=heatmap_title,
              cluster_rows=TRUE, clustering_method_rows="ward.D2",
              cluster_columns = TRUE,
              show_column_names=TRUE,
              row_names_side = c("left"),
              row_names_gp = gpar(fontsize = 5),
              column_names_gp = gpar(fontsize = 6),
              column_names_rot = 45,
              col=heatmap_col,
              heatmap_legend_param = list(
                title = "expression", 
                labels_gp = gpar(fontsize =6),
                title_gp = gpar(fontsize = 7, fontface = "bold")
              )
  )
  pdf(file)
  draw(ht)
  dev.off()
}

heatmap_afterTrajectory4avgExpr(t(scale(t(expr.NFkappaB.avg))),heatmap_title="NF-kappa B signaling pathway",heatmap_col=mycolor,file="heatmap_NFkappaB_avgExpr.pdf") 
  
## similar code apply to Hippo signaling pathway

