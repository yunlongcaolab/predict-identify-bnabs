library(ggplot2)
library(tidyverse)
library(ggrastr)
library(ggrepel)

# Fig. 2a lines
data <- read.csv(paste("../processed_source_data/Fig2/site-calc-sum-WT.csv",sep=""))

use_src <- "WT"
use_weight <- "WT"

pdf(paste0("plots/Fig2/Fig2a-Rplot-calc-site-scores-",use_src,".pdf"), width=5, height=2)
p <- ggplot(data %>% filter(absrc==use_src & weight == use_weight), aes(site, mut_escape_adj)) + 
    geom_line(color="#A03429", size=0.8, alpha=0.8) + geom_point(color="#A03429", shape=21)+ theme_classic() + theme(
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        # axis.text.x=element_blank(),
        # axis.ticks.x=element_blank(),
        axis.title=element_blank(),
        axis.text.x=element_text(angle=60, hjust=0.5, vjust=0.5)
    )+scale_x_continuous(breaks=seq(331,519,5),limits=c(331,519))+scale_color_manual(values=colors)+
geom_label_repel(data=data %>% filter(absrc == use_src & weight == use_weight & (mut_escape_adj > 0.13)) %>% group_by(site) %>% summarise(mut_escape_adj=max(mut_escape_adj)), 
                    aes(label=site), min.segment.length = 0, direction="both", fill = alpha(c("white"),0.5), max.overlaps=50)
print(p)
dev.off()


library(ComplexHeatmap)
library(circlize)

# Fig. 2c and Extend
show_cols <- c("D614G_IC50", "D614G+E484K", "D614G-S1", "D614G-S2", "D614G-S3", "D614G-S4", "D614G-S5", "BA1_IC50", "BA2_IC50", "BA5_IC50", "BQ1_1_IC50", "XBB1_5_IC50", 'HK3_1_IC50','JN1_IC50', "filter_worst")
data_raw <- read.csv("../source_data/Fig2-WT-NAbs-neut.csv", check.names=F)
rownames(data_raw) <- data_raw$id
data <- data_raw %>% 
    pivot_longer(!c(id, source,v_gene_H))

data$value[data$value > 10] = 10
data$value[data$value < 0.0005] = 0.0005

wtmat <- as.data.frame(data %>% pivot_wider(id_cols=name, names_from = id, values_from = value))

rownames(wtmat) <- wtmat$name
wtmat <- log10(wtmat[,-1])[show_cols,]
x <- hclust(dist(t(wtmat[c(show_cols, "JN1_IC50", "filter_worst"),]<log10(0.05))))

rownames(wtmat) <- c('B.1 (D614G)', 'B.1-E484K', 'B.1-S1', 'B.1-S2', 'B.1-S3', "B.1-S4", "B.1-S5", 'BA.1', 'BA.2', 'BA.5', 'BQ.1.1', 'XBB.1.5', 'HK.3.1', 'JN.1', 'S1-S5')

col_fun = colorRamp2(c(-3,log10(0.05),0), c("#4575b4","white","#d73027"))
col_fun_bin = colorRamp2(c(0, 0.3,1), c("#7585c4","white","#cccccc"))
col_fun_fc = colorRamp2(c(-1, -log10(3), log10(3),2), c("#4575b4","white","white","#d73027"))

wtmat <- wtmat[,x$order]

sources <- data_raw[colnames(wtmat), 'source']

pdf("plots/Fig2/fig2c-heatmap.pdf", width=7, height=2.7)
Heatmap(
    wtmat, cluster_rows = F, cluster_columns = F, na_col = "#666666", show_column_names=F,
        column_split=sources, border=T,col=col_fun
)

wtbin <- apply(wtmat > log10(0.05), 2, as.numeric)
rownames(wtbin) <- rownames(wtmat)

Heatmap(
    wtbin, cluster_rows = F, cluster_columns = F, show_column_dend=F, na_col = "#666666", show_column_names=F,
        column_split=sources,border=T, col=col_fun_bin
)
wtmat_fc <- t(t(wtmat) - t(wtmat)[,1])
wtmat_fc <- wtmat_fc[-1,]
wtmat_fc[wtmat_fc > 5] <- 5
Heatmap(
    wtmat_fc, cluster_rows = F, cluster_columns = F, show_column_dend=F, na_col = "#666666", show_column_names=F,
        column_split=sources,border=T, col=col_fun_fc
)
dev.off()

# Fig. 2d

pdf("plots/Fig2/fig2d-proportion.pdf", width=2.8, height=2.5)
res <- read.csv("../processed_source_data/Fig2/fig2d_WT_only_True.csv")
ggplot(res, aes(x=indicator, y=value)) + geom_bar(data=res %>% filter(target == "all"), stat="identity", fill="#CDCDCD", width=0.8) +
     geom_bar(data=res %>% filter(target == "BA1_IC50_BA2_IC50"), stat="identity", fill="#F4B183", width=0.8) +
     geom_bar(data=res %>% filter(target == "BA1_IC50_BA2_IC50_BA5_IC50"), stat="identity", fill="#9BBB59", width=0.8)+
     geom_bar(data=res %>% filter(target == "BA1_IC50_BA2_IC50_BA5_IC50_BQ1_1_IC50_XBB1_5_IC50"), stat="identity", fill="#CF6B50", width=0.8)+
     geom_line(aes(group=target), linetype="dashed") + 
     geom_point(shape=21, fill="white") + 
    theme_classic()

res <- read.csv("../processed_source_data/Fig2/fig2d_WT_only_False.csv")
ggplot(res, aes(x=indicator, y=value)) + geom_bar(data=res %>% filter(target == "all"), stat="identity", fill="#CDCDCD", width=0.8) +
     geom_bar(data=res %>% filter(target == "BA1_IC50_BA2_IC50"), stat="identity", fill="#F4B183", width=0.8) +
     geom_bar(data=res %>% filter(target == "BA1_IC50_BA2_IC50_BA5_IC50"), stat="identity", fill="#9BBB59", width=0.8)+
     geom_bar(data=res %>% filter(target == "BA1_IC50_BA2_IC50_BA5_IC50_BQ1_1_IC50_XBB1_5_IC50"), stat="identity", fill="#CF6B50", width=0.8)+
     geom_line(aes(group=target), linetype="dashed") + 
     geom_point(shape=21, fill="white") + 
    theme_classic()
dev.off()

# Fig. 2e and Ext - p values
cbp = c('#3CB371', '#4682B4', '#5F9EA0', '#A0522D', '#F4A460')
pdf("plots/Fig2/fig2e-p-values-stat.pdf", width=6.5, height=3)
res <- read.csv("../processed_source_data/Fig2/fig2e_WT_only_False.csv")
res$filter <- factor(res$filter, levels=c('D614G+E484K', 'D614G-S1', 'D614G-S2', 'D614G-S3', 'D614G-S4', 'D614G-S5', 'D614G-S1_D614G-S2_D614G-S3_D614G-S4_D614G-S5'))
ggplot(res, aes(x=filter, y=-log10(pval)))+
    geom_bar(stat="identity", aes(fill=target), position=position_dodge(width = 0.7), width=0.6)+
    geom_hline(yintercept = -log10(0.01), linetype="dashed")+
    scale_fill_manual(values=cbp)+
    theme_classic()+coord_flip()

res <- read.csv("../processed_source_data/Fig2/fig2e_WT_only_True.csv")
res$filter <- factor(res$filter, levels=c('D614G+E484K', 'D614G-S1', 'D614G-S2', 'D614G-S3', 'D614G-S4', 'D614G-S5', 'D614G-S1_D614G-S2_D614G-S3_D614G-S4_D614G-S5'))
ggplot(res, aes(x=filter, y=-log10(pval)))+
    geom_bar(stat="identity", aes(fill=target), position=position_dodge(width = 0.7), width=0.6)+
    geom_hline(yintercept = -log10(0.01), linetype="dashed")+
    scale_fill_manual(values=cbp)+
    theme_classic()+coord_flip()
dev.off()

# Fig. 2f and Ext - enrichment curves
cbp = c('#8CB371', '#F8766D', '#00BFC4', '#3F6E88', '#A0522D', '#F4A460', '#EE0000')
pdf("plots/Fig2/enrichment_curve_wt.pdf", width=4, height=2)
res <- read.csv("../processed_source_data/Fig2/enrichment_curve_data_wt.csv")
for (tg in unique(res$target)) {
    p <- ggplot(res %>% filter(target == tg), aes(id,ratio))+
        geom_line(aes(color=metric))+
        scale_color_manual(values=cbp)+
        scale_fill_manual(values=cbp)+
        geom_ribbon(aes(ymin = low, ymax = high,fill=metric), alpha = 0.1)+
        theme_classic()+ylim(0, 1)+xlab('k')+ylab(paste('ratio of',tg,'-effective \n NAbs in top-k'))
    
    print(p)
}

dev.off()
pdf("plots/Fig2/Fig2f_enrichment_curve.pdf", width=3.5, height=2)
res <- read.csv("../processed_source_data/Fig2/enrichment_curve_data_wt_sars.csv")
for (tg in c("BA5_IC50", "XBB1_5_IC50")) {
    p <- ggplot(res %>% filter(target == tg & metric %in% c("D614G_IC50", "filter_worst")), aes(id,ratio))+
        geom_line(aes(color=metric))+
        scale_color_manual(values=cbp)+
        scale_fill_manual(values=cbp)+
        geom_ribbon(aes(ymin = low, ymax = high,fill=metric), alpha = 0.1)+
        theme_classic()+ylim(0, 1)+xlab('k')+ylab(paste('ratio of',tg,'-effective \n NAbs in top-k'))
    
    print(p)
}
dev.off()