library(tidyverse)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)

# ExtFig2
show_cols <- c(
    "D614G_IC50", 
    "B.1+E484K", 
    "B.1+V483G",
    "B.1+F486I",
    "B.1+F490L",
    "B.1+F490S",
    "B.1+K378E",
    "B.1+K378N",
    "B.1+K417N",
    "B.1+L452R",
    "B.1+N450K",
    "B.1+Q493K",
    "B.1+R346S",
    "B.1+R346T",
    "B.1+K444E",
    "B.1+K444N",
    "B.1+V445D",
    "B.1+G446V",
    "BA1_IC50",
    "BA2_IC50",
    "BA5_IC50",
    "filter_worst_single"
)
data_raw <- read.csv("../processed_source_data/ExtFig2/neutralization_single.csv", check.names=F)
rownames(data_raw) <- data_raw$id
data <- data_raw %>% 
    pivot_longer(!c(id, source,v_gene_H)) %>% filter(name %in% show_cols)

data$value[data$value < 0.0005] = 0.0005

wtmat <- as.data.frame(data %>% pivot_wider(id_cols=name, names_from = id, values_from = value))

rownames(wtmat) <- wtmat$name
wtmat <- log10(wtmat[,-1])[show_cols,]
x <- hclust(dist(t(wtmat[c('filter_worst_single', 'BA1_IC50','BA2_IC50', 'BA5_IC50'),]<log10(0.05))))

col_fun = colorRamp2(c(-3,log10(0.05),0), c("#4575b4","white","#d73027"))
col_fun_bin = colorRamp2(c(0, 0.3,1), c("#7585c4","white","#cccccc"))
col_fun_fc = colorRamp2(c(-1, -log10(3), log10(3),2), c("#4575b4","white","white","#d73027"))

wtmat <- wtmat[,x$order]

sources <- (data_raw[colnames(wtmat), 'B.1+E484K'] < 0.05)

pdf("plots/ExtFig2/single-heatmap.pdf", width=7, height=4.2)
Heatmap(
    wtmat, cluster_rows = F, cluster_columns = F, na_col = "#666666", show_column_names=F,
        column_split=sources, border=T,col=col_fun
)

dev.off()
pdf("plots/ExtFig2/enrich_bars.pdf", width=2, height=2.5)
res <- read.csv("../processed_source_data/ExtFig2/enrich_bars_data.csv")

res$indicator <- factor(res$indicator, levels=c("D614G_IC50", "B.1+E484K", "filter_worst_single"))
ggplot(res, aes(x=indicator, y=value)) + geom_bar(data=res %>% filter(target == "all"), stat="identity", fill="#CDCDCD", width=0.8) +
     geom_bar(data=res %>% filter(target == "BA1_IC50_BA2_IC50"), stat="identity", fill="#F4B183", width=0.8) +
     geom_bar(data=res %>% filter(target == "BA1_IC50_BA2_IC50_BA5_IC50"), stat="identity", fill="#9BBB59", width=0.8)+
     geom_bar(data=res %>% filter(target == "BA1_IC50_BA2_IC50_BA5_IC50_BQ1_1_IC50_XBB1_5_IC50"), stat="identity", fill="#CF6B50", width=0.8)+
     geom_line(aes(group=target), linetype="dashed") + 
     geom_point(shape=21, fill="white") + 
    theme_classic()

dev.off()
# ExtFig5
plot <- function (fit_data, fit_params, fit_params_ind, data, name) {
    cb <- c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#e6ab02", "#a6761d", "#666666","#006d2c","#54278f","#810f7c","#08519c","#807dba","#6baed6", "#2caddc", "#9e5ac8", "#f33272")

    pdf(paste0(name, "_curve.pdf"), width=7, height=4)
    print(ggplot(fit_data, aes(x=conc, y=inhibition))+geom_line(aes(color=variant), alpha=0.7)+
        geom_point(data=data, aes(fill=variant), shape=21, size=3, alpha=0.5)+
        scale_x_log10()+
        scale_color_manual(values=cb)+
        scale_fill_manual(values=cb)+
        theme_classic())
    dev.off()

    pdf(paste0(name, "_EC50_bar.pdf"), width=round(nrow(fit_params)/2)+2, height=4)
    print(ggplot(fit_params, aes(variant, EC50*1000))+geom_bar(aes(fill=variant), stat="identity", color='black', width=.8)+
        geom_point(data=fit_params_ind, aes(variant, EC50*1000), shape=21, fill="white", size=2)+scale_y_log10()+
        geom_text(data=fit_params, aes(y=EC50*1000*1.8, label=signif(EC50,2)))+
        scale_fill_manual(values=cb)+
        theme_classic()+
        theme(axis.text.x=element_text(angle=45, hjust=1),axis.title.x=element_blank()))
    dev.off()
}

dir.create('plots/ExtFig5')
exclude <- c("BA.1", "BA.2", "HV.1", "XBB.1.5+FLip+D420Y+Y489H","JN.4", "JN.5.1", 'KP.3', "JN.1+R346T","BA.2.86+FLip", "XBB.1.5+S490Y", "JN.1+R346T+F456L+S490Y")

fit_data <- read.csv("../processed_source_data/ExtFig5/ACE2_neut_fit_data.csv") %>% filter(!(variant %in% exclude))
data <- read.csv("../source_data/ExtFig5-ACE2_neutralization.csv") %>% filter(!(variant %in% exclude))
fit_params <- read.csv("../processed_source_data/ExtFig5/ACE2_neut_fit_params.csv") %>% filter(!(variant %in% exclude))
fit_params_ind <- read.csv("../processed_source_data/ExtFig5/ACE2_neut_fit_params_ind.csv") %>% filter(!(variant %in% exclude))

plot(fit_data, fit_params,fit_params_ind, data, 'plots/ExtFig5/ACE2_neut')

# ExtFig6e
use_v <- c("XBB.1.5", "XBB.1.5+S490Y", "JN.1+R346T+F456L+S490Y")
dir.create('plots/ExtFig6')

fit_data <- read.csv("../processed_source_data/ExtFig5/ACE2_neut_fit_data.csv") %>% filter(variant %in% use_v)
data <- read.csv("../source_data/ExtFig5-ACE2_neutralization.csv") %>% filter(variant %in% use_v)
fit_params <- read.csv("../processed_source_data/ExtFig5/ACE2_neut_fit_params.csv") %>% filter(variant %in% use_v)
fit_params_ind <- read.csv("../processed_source_data/ExtFig5/ACE2_neut_fit_params_ind.csv") %>% filter(variant %in% use_v)

plot(fit_data, fit_params, fit_params_ind, data, 'plots/ExtFig6/ExtFig6e-ACE2_neut')

# ExtFig5c - SPR
library(circlize)
library(ggrepel)

data <- read.csv("../source_data/SPR_data_RBD_1205.csv") %>% filter(antibody == "BD55-1205")
data$KD[data$KD < 1e-12] <- 1e-12
use_ag <- (data %>% group_by(antibody,antigen) %>% summarise(count=n()) %>% filter(count > 1))$antigen
data <- data %>% filter(antigen %in% use_ag)
data$lgKD <- log10(data$KD)

cgcolors <- colorRamp2(c(1,5,12,15,20),c('#113C80', "#C9E7F1",'#F6e9dd', '#FBcEc3','#Ec684F'))(1:20)

data_gmean <- as.data.frame(data %>% group_by(antigen) %>% summarise(ka_gmean=exp(mean(log(ka))), kd_gmean=exp(mean(log(kd))), KD_gmean=exp(mean(log(KD)))))
use_levels <- data_gmean[order(data_gmean$KD_gmean),'antigen']
data$antigen <- factor(data$antigen, levels=use_levels)

rownames(data_gmean) <- data_gmean$antigen

data_gmean$label <- signif(1e9*data_gmean$KD_gmean, 2)
antigens <- data_gmean$antigen

data_plot <- data %>% group_by(antigen) %>% mutate(
    rank = rank(runif(n())),
    jitter_x = (rank / (n() + 1)) * 0.5 - 0.25
)

pdf("plots/ExtFig5/ExtFig5c-SPR_bar_1205.pdf",width=8, height=3.5)
ggplot(data_plot, aes(x=antigen, y=lgKD+14, fill=antigen))+stat_summary(fun='mean', geom='bar',show.legend = F, color='#999999',size=0.2)+
    geom_point(aes(x=as.numeric(antigen)+jitter_x),size=1.5, alpha=0.5,shape=21, fill='white', show.legend = F)+
    scale_y_continuous(breaks=c(2, 5), labels=c(expression('10'^{-12}),expression('10'^{-9})))+ylab('Apparent KD (M)')+
    coord_cartesian(ylim=c(0.5,7))+
    scale_fill_manual(values = cgcolors)+
    theme_classic()+theme(axis.text.x=element_text(angle=45, hjust=1,vjust=1),axis.title.x=element_blank())+
    geom_text(data=data_gmean, aes(x=antigen, label=label, y=log10(KD_gmean)+14.8), size=3)
dev.off()

# ExtFig. 5d - sensorgram
home <- '../source_data/ExtFig5d-SPR-raw/'
antigens <- c('WT', 'BA.5', 'XBB.1.5', 'XBB.1.5.70', 'BA.2.86', 'JN.1')

data <- data.frame()
concs <- c("1.5625nM", "6.25nM", "12.5nM", "25nM", "50nM")
pdf('plots/ExtFig5/ExtFig5d-SPR_sensorgram.pdf', width=4.5,height=2.5)
for (ag in antigens) {
    file_d <- read.delim(paste(home, '/', ag, '-RBD.txt',sep=''))[,1:20]
    file <- data.frame()
    for (i in 1:length(concs)) {
        f1 <- file_d[,(i*4-3):(i*4)]
        colnames(f1) <- c("X", "Y", "Fitted_X", "Fitted_Y")
        f1$conc <- concs[i]
        file <- rbind(f1, file)
    }
    
    file$conc <- factor(file$conc, levels=concs)
    
    print(
        ggplot()+
            geom_line(data=file, aes(Fitted_X, Fitted_Y, group=conc), alpha=1.0, size=1.0,color="#DDDDDD")+
            geom_line(data=file, aes(X,Y,color=conc), alpha=0.9,size=0.5)+
            scale_color_manual(values=c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e"))+
            theme_classic()+
            annotate("text", size=3.5, x=250, y=max(file$Fitted_Y,na.rm=T)*0.85,hjust=0, label=as.expression(bquote(k[a]~"="~.(format(data_gmean[ag,"ka_gmean"]/1e5, digits=2))~phantom()%*% 10^5~M^-1~s^-1)))+
            annotate("text", size=3.5, x=250, y=max(file$Fitted_Y,na.rm=T)*0.75,hjust=0, label=as.expression(bquote(k[d]~"="~.(format(data_gmean[ag,"kd_gmean"]*1e4, digits=2))~phantom()%*% 10^-4~s^-1)))+
            annotate("text", size=3.5, x=250, y=max(file$Fitted_Y,na.rm=T)*0.45, hjust=0, label=as.expression(bquote(K[D]~"="~.(format(data_gmean[ag,"KD_gmean"]*1e9, digits=2))~"nM")))+
            ggtitle(paste(ag," RBD",sep=''))+xlab("Time (s)")+ylab("Relative response (RU)")
    )

    file$antigen <- ag
    data <- rbind(data)
}
dev.off()

# ExtFig 6a
data <- read.csv("../source_data/ExtFig6a-authentic-virus-neut.csv")

data$strain <- factor(data$strain, levels=c("WT", "BA.5.2.1", "FL.8", "XBB.1.5.6", "JN.3"))

summary_data <- data %>%
  group_by(strain, source) %>%
  summarize(geo_mean_IC50 = exp(mean(log(IC50))),  # Geometric mean
            geo_sd_IC50 = exp(sd(log(IC50))),      # Geometric SD, using log transformation
            .groups = 'drop')

pdf("plots/ExtFig6/ExtFig6a-authentic-neut.pdf", width=6,height=2.5)

# Create the dot plot with symmetric error bars on log scale
ggplot(data, aes(x=strain, y=IC50, color=strain)) +
  geom_errorbar(data=summary_data, aes(x=strain, y=geo_mean_IC50, ymin=geo_mean_IC50/geo_sd_IC50, ymax=geo_mean_IC50*geo_sd_IC50), width=0.2, color='black') +  # Add symmetric error bars
  geom_errorbar(data=summary_data, aes(x=strain, y=geo_mean_IC50, ymin=geo_mean_IC50, ymax=geo_mean_IC50), width=0.4, color='black') +  # Add symmetric error bars
  scale_y_log10(limits=c(0.001,0.1), breaks=c(seq(0.001,0.01,0.001), seq(0.01,0.1,0.01)), ) +  # Log scale for the y-axis
  geom_point(position=position_jitter(width=0.2), size=3, alpha=.8, shape=21, show.legend=F) +  # Draw the individual points
  labs(x="", y="Neutralization IC50 (µg/mL)") +  # Axis labels
    geom_text(data=summary_data, aes(x=strain,label=sprintf("%.3f",geo_mean_IC50)), y=log10(0.08), color='black')+
    scale_color_manual(values=c("#4D80E4", "#E44D7B", "#3AA361", "#994021", "#7E57C2"))+
  theme_classic()  # Minimal theme
dev.off()

# ExtFig 6d
fit_data <- read.csv("../processed_source_data/ExtFig6/490_neut_fit_data.csv")
data <- read.csv("../processed_source_data/ExtFig6/490_neut_data_clean.csv")
fit_params <- read.csv("../processed_source_data/ExtFig6/490_neut_fit_params.csv")
fit_params_ind <- read.csv("../processed_source_data/ExtFig6/490_neut_fit_params_ind.csv")

cb <- c("#1b9e77",  "#7570b3")

pdf(paste0("plots/ExtFig6/ExtFig6d-fit_curve.pdf"), width=7, height=4)
print(ggplot(fit_data, aes(x=conc, y=inhibition))+geom_line(aes(color=variant), alpha=0.7)+
    geom_point(data=data, aes(fill=variant), shape=21, size=3, alpha=0.5)+
    scale_x_log10()+
    scale_color_manual(values=cb)+
    scale_fill_manual(values=cb)+
    theme_classic())
dev.off()

# ExtFig 9e
data <- read.csv("../source_data/SPR_data_RBD_1205.csv")
data$KD[data$KD < 1e-12] <- 1e-12

data$lgKD <- log10(data$KD)

data_gmean <- as.data.frame(data %>% group_by(antibody, antigen) %>% summarise(ka_gmean=exp(mean(log(ka))), kd_gmean=exp(mean(log(kd))), KD_gmean=exp(mean(log(KD)))))
use_ag <- (data_gmean %>% pivot_wider(id_cols=antigen, names_from=antibody, values_from="KD_gmean") %>% na.omit())$antigen

data <- data %>% filter(antigen %in% use_ag)
data_gmean <- data_gmean %>% filter(antigen %in% use_ag)

data$antigen <- factor(data$antigen, levels=c('WT', 'BA.5', 'XBB.1.5', 'XBB.1.5.70', 'JN.1'))
data_gmean$antigen <- factor(data_gmean$antigen, levels=c('WT', 'BA.5', 'XBB.1.5', 'XBB.1.5.70', 'JN.1'))
data_gmean$label <- signif(1e12*data_gmean$KD_gmean, 2)

dir.create("plots/ExtFig9")
pdf("plots/ExtFig9/ExtFig9e-SPR_bar_1205_rev.pdf",width=6, height=2)
ggplot(data, aes(x=antigen, y=lgKD+14))+
    geom_bar(data=data_gmean, aes(y=log10(KD_gmean)+14, fill=antibody), stat="identity", width=.8, position=position_dodge(.9))+
    geom_hline(yintercept = 2, linetype='dashed')+
    geom_point(aes(fill=antibody), position=position_jitterdodge(jitter.width = .1, dodge.width = .9, seed=42), size=1.5, alpha=0.7,shape=21, show.legend = F)+
    scale_y_continuous(breaks=c(2, 5), labels=c(expression('10'^{-12}),expression('10'^{-9})))+ylab('Apparent KD (M)')+
    coord_cartesian(ylim=c(1.5,7))+
    scale_fill_manual(values = c("#264653", "#2a9d8f", "#e76f51", "#e9c46a", "#f4a261"))+
    theme_classic()+theme(axis.text.x=element_text(angle=45, hjust=1,vjust=1),axis.title.x=element_blank())

dev.off()
