
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

################################################
## function for UMAP feature plot for scRNA-seq data from human LN
##### (data download from PMID: 33903764) #####
################################################

marker.plot= function(markers,counts.table,embedding){
  gene.expr.selected=counts.table[row.names(counts.table)==markers,]
  cell.selected=names(gene.expr.selected)[which(gene.expr.selected>0)]
  embedding.selected=embedding[cell.selected,]
  gene.expr.cell.selected=gene.expr.selected[which(gene.expr.selected>0)]
  gene.expr.cell.selected.lg=gene.expr.cell.selected
  gene.expr.cell.selected.lg=as.vector(gene.expr.cell.selected.lg)
  mi=min(gene.expr.cell.selected.lg)
  ma=max(gene.expr.cell.selected.lg)
  v <- round((gene.expr.cell.selected.lg- mi)/(ma - mi)*99 + 1,0)
  gene.selected=data.frame(embedding.selected,expr=gene.expr.cell.selected,expr.modify=v)
  p<-ggplot(embedding, aes(x=UMAP1_orig, y=UMAP2_orig)) + 
    geom_point(size=0.3,color="grey")+
    geom_point( data=gene.selected,mapping=aes( x=UMAP1,y=UMAP2,color=(expr.modify)),size=0.3)+
    scale_colour_gradientn(colours = myPalette(100))+
    theme_bw()+
    ggtitle(markers)+
    theme(
      legend.position="right",
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(),
      panel.background = element_blank(), 
      axis.line =element_blank(),
      panel.border = element_rect(colour = "black", fill=NA, size=1),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      plot.title = element_text(hjust = 0.5,face="italic"),
      text = element_text(size=8)
    )+
    guides(color = guide_colourbar(title="",barwidth = 0.7, barheight =5))
  print(p)
}

marker.plot(markers,counts,emb)
