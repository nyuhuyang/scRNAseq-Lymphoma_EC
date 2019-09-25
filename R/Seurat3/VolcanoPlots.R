########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################

library(ggplot2)
library(ggrepel)
rm(list=ls());gc()
path <- paste0("./output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)

# 3.1.1 load data ==================
dates <- c("2018-10-18","2018-12-30")
cluster <- c("EC+RO2", "EC+3119")
dataset <- list()
for(i in seq_along(dates)){
        dataset[[i]] = read.csv(file=paste0("output/20190924/EC_markers_",dates[i],".csv"),
                            row.names = 1,stringsAsFactors=F)
        dataset[[i]]$cluster = cluster[i]
        dataset[[i]]$gene = rownames(dataset[[i]])
        jpeg(paste0(path,"Volcano_plot_EC_",clusters[i],".jpeg"), units="in", width=10, height=7,res=600)
        print(VolcanoPlots(dataset[[i]],alpha = 0.9,size = 1.5,cut_off_logFC = 0.5)+
                      ggtitle(paste(cluster[i]," vs. naive EC"))+
                      theme(plot.title = element_text(size=15, hjust = 0.5,face="plain")))
        dev.off()
}


VolcanoPlots <- function(data, cut_off_pvalue = 0.0000001, cut_off_logFC = 0.25,
                         cols = c("#0000ff","#d2dae2","#ff0000"),alpha=0.8, size=2,
                         legend.size = 12){
        data = data[data$p_val<0.05,]
        rownames(data) = data$gene
        data$change = ifelse(data$p_val_adj < cut_off_pvalue &
                                        abs(data$avg_logFC) >= cut_off_logFC, 
                                ifelse(data$avg_logFC > cut_off_logFC ,'Up','Down'),
                                'Stable')
        colnames(data)[grep("cluster",colnames(data))]="cluster"
        # 将需要标记的基因放置在单独的数组
        Up <- data[data$change %in% "Up",]
        Down <- data[data$change %in% "Down",]
        Up_gene_index <- rownames(Up)[Up$p_val_adj <= head(sort(Up$p_val_adj,decreasing = F),15) %>% tail(1)]
        Down_gene_index <- rownames(Down)[Down$p_val_adj <= head(sort(Down$p_val_adj,decreasing = F),15) %>% tail(1)]
        p<-ggplot(
                #设置数据
                data, 
                aes(x = avg_logFC, 
                    y = -log10(p_val_adj), 
                    colour=change))+
                
                # 辅助线
                geom_vline(xintercept=c(-cut_off_logFC,cut_off_logFC),lty=4,col="black",lwd=0.8) +
                geom_hline(yintercept = -log10(cut_off_pvalue),lty=4,col="black",lwd=0.8) +
                
                # 坐标轴
                labs(x="log2(fold change)",
                     y="-log10 (p-value)")+
                theme_bw()+
                
                # 图例
                theme(plot.title = element_text(hjust = 0.5), 
                      axis.title=element_text(size=12),
                      legend.position="right", 
                      legend.title = element_blank(),
                      legend.text = element_text(size = legend.size),
                )
        
        p = p + ggrepel::geom_text_repel(data = data[c(Down_gene_index[1:15],Up_gene_index[1:15]),], 
                              aes(label = gene),
                          size = 3,box.padding = unit(0.5, "lines"),
                          point.padding = unit(0.8, "lines"), 
                          segment.color = "black", 
                          show.legend = FALSE)+
                geom_point(alpha=alpha, size=size) +
                scale_color_manual(values=cols)
        return(p)
}


