library(ggplot2)
library(tidyverse)
library(ggrastr)

data <- read.csv("../processed_source_data/Fig1/selected_Ab_clean_IC50.csv")

data$source[data$source == "BA.1 BTI + BA.5/BF.7" | data$source == "BA.2 BTI + BA.5/BF.7"] <- "BA.1/BA.2 BTI + BA.5/BF.7"
data$IC50[data$IC50 > 10] = 10
data$IC50[data$IC50 < 0.0005] = 0.0005

# show_variants <- c('D614G', 'BA.1', 'BA.2', 'BA.5', 'BQ.1.1', 'XBB.1.5')
# show_variants <- c('D614G', 'BA.1', 'BA.2', 'BA.5', 'BQ.1.1', 'XBB.1.5', 'SARS')
show_variants <- c('SARS', 'D614G', 'BA.1', 'BA.2', 'BA.5', 'BQ.1.1', 'XBB.1.5', 'HK.3.1', 'JN.1', 'KP.3')

# show_variants <- c('SARS', 'D614G', 'BA.1', 'BA.2', 'BA.5', 'BQ.1.1', 'XBB.1.5', 'HK.3.1', 'JD.1.1', 'BA.2.86', 'JN.1')
data <- data %>% filter(variant %in% show_variants)
data$variant <- factor(data$variant, levels=show_variants)

data$is_neut <- data$IC50 < 0.05
stat <- data %>% group_by(source, variant) %>% na.omit() %>% summarise(count=n(), neut=sum(is_neut))

# stat$label <- paste(stat$neut, '/', stat$count, '\n',round(stat$neut/stat$count*100), '%', sep='')
stat$label <- paste0(stat$neut, '\n',round(stat$neut/stat$count*100), '%')



colors <- c("SARS"="#E69F00", "D614G"="#56B4E9", 
            "BA.1"="#009E73", "BA.2"="#F0E442", "BA.5"="#0072B2", "BQ.1.1"="#D55E00", 
            "XBB.1.5"="#CC79A7", "HK.3.1"="#BB99AC", "JD.1.1"="#CCBD99", "BA.2.86"="#3CA091", "JN.1"="#1B9E77")

dir.create("plots/Fig1", showWarnings = FALSE)
pdf("plots/Fig1/fig1c-h-sars.pdf", width=4.4, height=1.8)
for (src in unique(data$source)) {
p <- ggplot(data %>% filter(source == src & variant %in% c('SARS', 'D614G', 'BA.1', 'BA.2', 'BA.5', 'BQ.1.1', 'XBB.1.5', 'HK.3.1', 'JN.1','KP.3')), aes(variant, log10(IC50)))+geom_hline(yintercept=1, linetype='dashed')+
    geom_hline(yintercept=log10(0.0005), linetype='dashed')+geom_hline(yintercept=log10(0.05), color='red', linetype='dashed')+
    geom_line(aes(group=id), alpha=0.01)+
    geom_point_rast(aes(color=variant), alpha=0.05, size=1, show.legend = F)+
    stat_summary(fun=mean, geom="point", shape=21, size=3, color="black") + labs(y = "Pseudovirus IC50 (\u00b5g/mL)")+
    scale_color_manual(values=colors)+
    scale_y_reverse(limits=c(1,-5), breaks=c(1,0,-1,-2,-3), labels=c(expression('10'^{1}),expression('10'^{0}),expression('10'^{-1}),expression('10'^{-2}),expression('10'^{-3})))+
    theme_classic()+ggtitle(
        paste0(src, ' (n=', (stat%>%filter(source == src))$count[1], ')')
    )+theme(title = element_text(size=8), axis.title.x=element_blank(), axis.text.x=element_text(angle=45, hjust=1))

p <- p+geom_text(data=stat%>%filter(source == src), aes(x=variant, label=label), y=4.4, size=3)
print(p)
}
dev.off()
pdf("plots/Fig1/fig1c-h.pdf", width=4, height=1.8)
for (src in unique(data$source)) {
p <- ggplot(data %>% filter(source == src & variant %in% c('D614G', 'BA.1', 'BA.2', 'BA.5', 'BQ.1.1', 'XBB.1.5', 'HK.3.1', 'JN.1','KP.3')), aes(variant, log10(IC50)))+geom_hline(yintercept=1, linetype='dashed')+
    geom_hline(yintercept=log10(0.0005), linetype='dashed')+geom_hline(yintercept=log10(0.05), color='red', linetype='dashed')+
    geom_line(aes(group=id), alpha=0.01)+
    geom_point_rast(aes(color=variant), alpha=0.05, size=1, show.legend = F)+
    stat_summary(fun=mean, geom="point", shape=21, size=3, color="black") + labs(y = "Pseudovirus IC50 (\u00b5g/mL)")+
    scale_color_manual(values=colors)+
    scale_y_reverse(limits=c(1,-5), breaks=c(1,0,-1,-2,-3), labels=c(expression('10'^{1}),expression('10'^{0}),expression('10'^{-1}),expression('10'^{-2}),expression('10'^{-3})))+
    theme_classic()+ggtitle(
        paste0(src, ' (n=', (stat%>%filter(source == src))$count[1], ')')
    )+theme(title = element_text(size=8), axis.title.x=element_blank(), axis.text.x=element_text(angle=45, hjust=1))

p <- p+geom_text(data=stat%>%filter(source == src & variant != "SARS"), aes(x=variant, label=label), y=4.4, size=3)
print(p)
}
dev.off()

use_src <- c('WT', 'SARS+WT', 'BA.1 BTI', 'BA.2 BTI', 'BA.5/BF.7 BTI', 'BA.1 BTI + BA.5/BF.7', 'BA.2 BTI + BA.5/BF.7')
data_all <- read.csv("../processed_source_data/Fig1/all_mAbs_IC50.csv") %>% filter(source %in% use_src)

for (col in grep("_IC50", x= colnames(data_all), value=T)) {
    data_all[,col][data_all[,col] > 10] = 10
    data_all[,col][data_all[,col] < 0.0005] = 0.0005    
}
data_all$source <- factor(data_all$source, levels=use_src)
data_all %>% group_by(source) %>% summarise(count=n())
cbp = c('#2F4F4F', '#791351', '#4682B4', '#5F9EA0', '#A0522D', '#F4A460', '#C93636')
# names(cbp) <- c(use_src, 'SARS')
names(cbp) <- use_src

pdf("plots/Fig1/fig1b.pdf", width=5, height=3)
ggplot(data_all %>% filter(auto_IC50 < 0.05), aes(log10(auto_IC50), log10(XBB1_5_IC50)))+
    geom_point_rast(data=data_all %>% filter(auto_IC50 > 0.05), color="#AAAAAA", alpha=0.3)+
    scale_color_manual(values=cbp)+
    geom_point_rast(aes(color=source),alpha=0.6)+
    geom_hline(yintercept = log10(0.05), linetype="dashed")+geom_vline(xintercept = log10(0.05), linetype="dashed")+theme_classic()+
    scale_x_reverse(limits=c(1,-4), breaks=c(1,0,-1,-2,-3), labels=c(expression('10'^{1}),expression('10'^{0}),expression('10'^{-1}),expression('10'^{-2}),expression('10'^{-3})))+
    scale_y_reverse(limits=c(1,-4), breaks=c(1,0,-1,-2,-3), labels=c(expression('10'^{1}),expression('10'^{0}),expression('10'^{-1}),expression('10'^{-2}),expression('10'^{-3})))+
    theme(aspect.ratio = 1.0)
dev.off()
pdf("plots/Fig1/extfig1a-b.pdf", width=7, height=3)
ggplot(data_all %>% filter(auto_IC50 < 0.05), aes(log10(XBB1_5_IC50)))+
    geom_density(aes(fill=source), color="black", alpha=0.5)+
    scale_fill_manual(values=cbp)+
    geom_vline(xintercept = log10(0.05), linetype="dashed")+theme_classic()+
    scale_x_reverse(limits=c(1,-4), breaks=c(1,0,-1,-2,-3), labels=c(expression('10'^{1}),expression('10'^{0}),expression('10'^{-1}),expression('10'^{-2}),expression('10'^{-3})))


ggplot(data_all %>% filter(auto_IC50 < 0.05), aes(log10(JN1_IC50)))+
    geom_density(aes(fill=source), color="black", alpha=0.5)+
    scale_fill_manual(values=cbp)+
    geom_vline(xintercept = log10(0.05), linetype="dashed")+theme_classic()+
    scale_x_reverse(limits=c(1,-4), breaks=c(1,0,-1,-2,-3), labels=c(expression('10'^{1}),expression('10'^{0}),expression('10'^{-1}),expression('10'^{-2}),expression('10'^{-3})))

ggplot(data_all %>% filter(auto_IC50 < 0.05), aes(log10(KP3_IC50)))+
    geom_density(aes(fill=source), color="black", alpha=0.5)+
    scale_fill_manual(values=cbp)+
    geom_vline(xintercept = log10(0.05), linetype="dashed")+theme_classic()+
    scale_x_reverse(limits=c(1,-4), breaks=c(1,0,-1,-2,-3), labels=c(expression('10'^{1}),expression('10'^{0}),expression('10'^{-1}),expression('10'^{-2}),expression('10'^{-3})))
ggplot(data_all %>% filter(auto_IC50 < 0.05), aes(log10(auto_IC50), log10(JN1_IC50)))+
    geom_point_rast(data=data_all %>% filter(auto_IC50 > 0.05), color="#AAAAAA", alpha=0.3)+
    scale_color_manual(values=cbp)+
    geom_point_rast(aes(color=source),alpha=0.6)+
    geom_hline(yintercept = log10(0.05), linetype="dashed")+geom_vline(xintercept = log10(0.05), linetype="dashed")+theme_classic()+
    scale_x_reverse(limits=c(1,-4), breaks=c(1,0,-1,-2,-3), labels=c(expression('10'^{1}),expression('10'^{0}),expression('10'^{-1}),expression('10'^{-2}),expression('10'^{-3})))+
    scale_y_reverse(limits=c(1,-4), breaks=c(1,0,-1,-2,-3), labels=c(expression('10'^{1}),expression('10'^{0}),expression('10'^{-1}),expression('10'^{-2}),expression('10'^{-3})))+
    theme(aspect.ratio = 1.0)
ggplot(data_all %>% filter(auto_IC50 < 0.05), aes(log10(auto_IC50), log10(KP3_IC50)))+
    geom_point_rast(data=data_all %>% filter(auto_IC50 > 0.05), color="#AAAAAA", alpha=0.3)+
    scale_color_manual(values=cbp)+
    geom_point_rast(aes(color=source),alpha=0.6)+
    geom_hline(yintercept = log10(0.05), linetype="dashed")+geom_vline(xintercept = log10(0.05), linetype="dashed")+theme_classic()+
    scale_x_reverse(limits=c(1,-4), breaks=c(1,0,-1,-2,-3), labels=c(expression('10'^{1}),expression('10'^{0}),expression('10'^{-1}),expression('10'^{-2}),expression('10'^{-3})))+
    scale_y_reverse(limits=c(1,-4), breaks=c(1,0,-1,-2,-3), labels=c(expression('10'^{1}),expression('10'^{0}),expression('10'^{-1}),expression('10'^{-2}),expression('10'^{-3})))+
    theme(aspect.ratio = 1.0)

dev.off()

for (var in c('XBB1_5_IC50', 'JN1_IC50', 'KP3_IC50')){
    data_all[is.na(data_all[,var]),var] <- Inf
    data_all$auto_pass <- data_all$auto_IC50 < 0.05
    data_all$var_pass <- data_all[,var] < 0.05
    print(var)
    print(data_all %>% group_by(auto_pass, var_pass) %>% summarise(count=n()))
}

pdf("plots/Fig1/Fig1a.pdf", width=4.5, height=1.5)
ggplot(data_all, aes(source, log10(auto_IC50)))+
    geom_dotplot(aes(color=source), binaxis = "y", alpha=0.3, stackdir = "center", binwidth = 0.05, show.legend = F)+
    stat_summary(fun=mean, geom="point", shape=21, size=3, color="black") + labs(y = "Pseudovirus IC50 (\u00b5g/mL)")+
    scale_color_manual(values=cbp)+geom_hline(yintercept=1, linetype='dashed')+
    geom_hline(yintercept=log10(0.0005), linetype='dashed')+geom_hline(yintercept=log10(0.05), color='red', linetype='dashed')+
    scale_y_reverse(limits=c(1,-4), breaks=c(1,0,-1,-2,-3), labels=c(expression('10'^{1}),expression('10'^{0}),expression('10'^{-1}),expression('10'^{-2}),expression('10'^{-3})))+
    theme_classic()
dev.off()
