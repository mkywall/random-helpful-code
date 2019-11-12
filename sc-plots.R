library(Seurat)
library(ggplot2)
library(dplyr)
library(circlize)
library(tidyr)
library(UpSetR)
library(ggrepel)

get_overlap_deg<- function(obj1, obj2, group1name, group2name, colno){
  
  # Upset plots of overlap between gaba 
  obj1$group<- group1name
  obj2$group<- group2name
  #obj3$group<- group3name
  #obj1$X<- rownames(obj1)
  #obj2$X<- rownames(obj2)
  
  df<- rbind(obj1,
             obj2)
  # df<- rbind(df,
  #            obj3)
  # 
  df<- df[, c(1, 3, colno)]
  df_spr = spread(df, group, avg_logFC)
  
  rownames(df_spr)<- df_spr$X
  df_spr$Symbol<- as.character(df_spr$X)
  df_spr$X<- NULL
  df_spr[is.na(df_spr)]<- 0
  df_spr[df_spr!= 0]<- 1
  df_spr$Symbol<- rownames(df_spr)
  rownames(df_spr)<- seq(1:nrow(df_spr))
  #setEPS()
  #postscript(paste0("./upset-",group1name, "--", group2name,".eps"), width = 7, height = 4)
  upset(df_spr, sets = c(group1name,
                         group2name), order.by = ("freq"))
  #dev.off()
}
setup_dg_for_chord<- function(file, fullgroup, posgroupname,neggroupname, negcolor, poscolor, genecat1, genecat2, genecat3){
  df<- read.csv(file,  stringsAsFactors = F)
  df<- df[which(df$p_val_adj < 0.05), c(1, 3,6)]
  df<- df[which(abs(df$avg_logFC) > 0.1),]
  df$group = fullgroup
  df$link_color = "black"
  df$category = "blank"
  df$category[which(df$X %in% genecat1) ]<- 1
  df$category[which(df$X %in% genecat2) ]<- 2
  df$category[which(df$X %in% genecat3) ]<- 3
  #df$category[which(df$category == "blank")]<- 4
  df<- df[which(df$category != "blank"), ]
  df$link_color[which(df$avg_logFC < 0)] = negcolor
  df$link_color[which(df$avg_logFC > 0)] = poscolor
  df$group[which(df$avg_logFC < 0)]<- neggroupname
  df$group[which(df$avg_logFC > 0)]<- posgroupname
  return(df)
}
setup_dg_for_chord2<- function(file, fullgroup, posgroupname,neggroupname, negcolor, poscolor, genecat){
  df<- read.csv(file,  stringsAsFactors = F)
  df<- df[which(df$p_val_adj < 0.05), c(1, 3,6)]
  df<- df[which(abs(df$avg_logFC) > 0.1), ]
  df$group = fullgroup
  df$link_color = "black"
  df$category = "blank"
  for ( i in 1:length(genecat)){
    df$category[which(df$X %in% genecat[[i]]) ]<- i }
  #df$category[which(df$category == "blank")]<- 4
  df<- df[which(df$category != "blank"), ]
  df$link_color[which(df$avg_logFC < 0)] = negcolor
  df$link_color[which(df$avg_logFC > 0)] = poscolor
  df$group[which(df$avg_logFC < 0)]<- neggroupname
  df$group[which(df$avg_logFC > 0)]<- posgroupname
  return(df)
}
make_volcano<- function(file, logfcthresh, genes, lim = NULL){
  allgenes<- read.csv(file, stringsAsFactors = F)
  allgenes$neglog10<- -log10(allgenes$p_val_adj)
  if(length(which(allgenes$neglog10 == Inf) > 0)){
    allgenes$neglog10[which(allgenes$neglog10 == Inf)]<- max(allgenes$neglog10[-which(allgenes$neglog10 == Inf)])
  }
  labeldf<- allgenes %>%
    dplyr::filter(p_val_adj < 0.01) %>%
    dplyr::filter(abs(avg_logFC) > logfcthresh) %>%
    dplyr::filter(X %in% genes) %>%
    dplyr::arrange(desc(abs(avg_logFC)))
  labeltheseup<- labeldf$X[which(labeldf$avg_logFC > 0)]
  labelthesedown<- labeldf$X[which(labeldf$avg_logFC < 0)]
  if(is.null(lim)){
    lowerlim<- min(allgenes$avg_logFC)
    upperlim<- max(allgenes$avg_logFC)
    
    lim<- max(abs(lowerlim), abs(upperlim))
    lim<- lim + 0.75}
  
  ggplot(allgenes, aes(avg_logFC, neglog10)) +
    geom_point(data = allgenes[which(allgenes$neglog10 < 1.3), ], colour = "grey", size = 1) +
    geom_point(data=allgenes[which(allgenes$avg_logFC < 0 & allgenes$neglog10 > 1.3), ], colour="blue", size=1) +
    geom_point(data=allgenes[which(allgenes$avg_logFC > 0 & allgenes$neglog10 > 1.3), ], colour="red", size=1) +
    
    
    
    
    geom_text_repel(data= allgenes[which(allgenes$X %in% labeltheseup), ],
                    aes(label=allgenes$X[which(allgenes$X %in% labeltheseup)]), 
                    force = 1,
                    direction = "y",
                    min.segment.length = 0,
                    segment.color = 'grey10',
                    segment.size = 0.5,
                    arrow = arrow(length = unit(0.01, 'npc'),
                                  type = 'closed', ends = 'first'),
                    size = 4, 
                    hjust = 0, 
                    nudge_x = lim - allgenes$avg_logFC[which(allgenes$X %in% labeltheseup)]) + 
    
    geom_text_repel(data= allgenes[which(allgenes$X %in% labelthesedown), ],
                    aes(label=allgenes$X[which(allgenes$X %in% labelthesedown)]), 
                    force = 1,
                    direction = "y",
                    min.segment.length = 0,
                    segment.color = 'grey10',
                    segment.size = 0.5,
                    arrow = arrow(length = unit(0.01, 'npc'),
                                  type = 'closed', ends = 'first'),
                    size = 4, 
                    hjust = 1, 
                    nudge_x = (-1*lim) + allgenes$avg_logFC[which(allgenes$X %in% labelthesedown)]) + 
    
    theme_classic()+
    theme(text = element_text(size = 20))+
    xlab('avg. LFC')+
    ylab('-log10(adj. p)')+
    xlim((-1*lim), lim) }
