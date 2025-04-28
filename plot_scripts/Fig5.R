library(ggplot2)
library(tidyverse)
library(scales)

gmean <- function(x){exp(mean(log(x)))}
gsd <- function(x){exp(sd(log(x)))}
gci_lower <- function(x){
    a <- log(x)
    logmean <- mean(a)
    logsd <- sd(a)
    
    exp(logmean - 1.96*(logsd / sqrt(length(x))))
}

gci_higher <- function(x){
    a <- log(x)
    logmean <- mean(a)
    logsd <- sd(a)
    
    exp(logmean + 1.96*(logsd / sqrt(length(x))))
}

plot_bar <- function(data, fit_params, pdfname) {
    samples <- unique(data$mouse)

    fit_params$tm <- log(2)/(fit_params$fit_k+fit_params$fit_k3)/24
    fit_params$ta <- log(2)/fit_params$fit_k3/24
    groups <- unique(fit_params$group)

    for (g in groups) {
        print(g)
        print(mean((fit_params %>% filter(group == g))$ta))
        print(sd((fit_params %>% filter(group == g))$ta))
    }

    print("24h")
    print(data %>% filter(hours == 24) %>% group_by(group) %>% summarise(gmt=gmean(conc*0.15)))
    print("48h")
    print(data %>% filter(hours == 48) %>% group_by(group) %>% summarise(gmt=gmean(conc*0.15)))

    use_breaks <- c(seq(10,100,10), seq(200, 1000, 100)) 

    data_summary <- data %>% filter(hours == 48) %>% group_by(group) %>% summarise(gmean_conc=gmean(conc*0.15), gsd_conc=gsd(conc*0.15), gcilow=gci_lower(conc*0.15), gcihigh=gci_higher(conc*0.15))
    pdf(pdfname, width=1.5, height=2)
    print(
        ggplot(data_summary, aes(x=group, y=gmean_conc))+
            geom_bar(stat='identity', aes(fill=group),color='black', width=0.8, show.legend = F)+
        scale_y_log10(limits=c(10,2000), breaks=use_breaks, oob = rescale_none)+
            geom_point(data=data %>% filter(hours == 48), aes(y=conc*0.15), shape=21, fill='white', position=position_jitter(0.3,seed=42))+scale_fill_manual(values=c("#1b9e77", "#A0492F", "#64A1EC"))+
            geom_errorbar(aes(ymin=gcilow, ymax=gcihigh), width=.5)+
            theme_classic()+
        theme(axis.title=element_blank(), axis.text=element_blank())
    )
    dev.off()
}

plot_fit <- function(data, fit_data, pdfname) {
    data_summary <- data  %>% group_by(hours, group) %>% summarise(gmean_conc=gmean(conc*0.15), gsd_conc=gsd(conc*0.15), gcilow=gci_lower(conc*0.15), gcihigh=gci_higher(conc*0.15))
    samples <- unique(data$group)
    use_breaks <- c(1:9, seq(10,100,10), seq(200, 1000, 100)) 

    pdf(pdfname, width=4.3, height=2)
    print(ggplot(fit_data, aes(hours, conc*0.15))+
        geom_point(data=data, aes(hours, conc*0.15,color=group), alpha=0.3)+
        geom_line(data=data, aes(hours, conc*0.15, group=mouse,color=group), alpha=0.1)+
        geom_line(aes(color=group))+
        scale_y_log10(limits=c(1, 2100), breaks=use_breaks)+
        scale_color_manual(values=c("#1b9e77", "#7570b3", "#1b9e77", "#A0492F", "#64A1EC"))+
        geom_point(data=data_summary, aes(hours, gmean_conc, color=group), shape=21, size=2)+
        geom_errorbar(data=data_summary, aes(x=hours, y=gmean_conc, ymin=gcilow, ymax=gcihigh,color=group))+
    theme_classic()+
        theme(axis.title=element_blank()))
    dev.off()
}

dir.create("plots/Fig5", recursive = T)
data <- read.csv("../source_data/Fig5-mice-data-elisa.csv", check.names=F, colClasses = c("numeric","character","character","numeric","character"))
fit_params <- read.csv("../processed_source_data/Fig5/ELISA_fit_params.csv", check.names=F)

data$group <- paste(data$type, data$sex, sep="-")

plot_bar(data %>% filter(group != "CHIK-24-F"), fit_params %>% filter(group!= "CHIK-24-F"), 'plots/Fig5/Ext-compare-Tg32.pdf')
plot_bar(data %>% filter(type != "Tg32"), fit_params %>% filter(group != "Tg32-F" & group!= "Tg32-M"), 'plots/Fig5/Fig5c_bar.pdf')

# Fig.5b fit
fit_data <- read.csv("../processed_source_data/Fig5/ELISA_fit_data_agg.csv", check.names=F) %>% filter(group != "Tg32-F" & group!= "Tg32-M")
plot_fit(data %>% filter(type != "Tg32"), fit_data %>% filter(group != "Tg32-F" & group!= "Tg32-M"), 'plots/Fig5/Fig5b_fit.pdf')

# Fig.5d-e neutralization
neut <- read.csv("../source_data/Fig5-mice-neutralization.csv", check.names=F, colClasses = c("numeric","character","character","numeric","character","character"))
neut$group <- paste(neut$type, neut$sex, sep="-")
neut$variant <- factor(neut$variant, levels=c("XBB.1.5", "HK.3.1", "JN.1"))
use_breaks <- c(seq(10,100,10), seq(200, 1000, 100)) 
neut <- neut  %>% filter(group == "Tg32-SCID-F")

pdf("plots/Fig5/Fig5e-neut_bar.pdf", width=2.2, height=2.5)
neut_stat <- neut %>% filter(hours==48) %>% group_by(group, variant) %>% summarise(gmean_titer=gmean(ID50), gsd_titer=gsd(ID50), gcilow=gci_lower(ID50), gcihigh=gci_higher(ID50))
print(neut_stat)

ggplot(neut_stat, aes(x=variant, y=gmean_titer, fill=variant))+
    geom_bar(stat="identity", color='black',width=.8,show.legend = F)+
    scale_fill_manual(values=c("#1b9e77", "#7570b3", "#d95f02"))+
    scale_y_log10(limits=c(100,10000), breaks=use_breaks*10, oob = rescale_none)+
    ylab('48h Serum ID50')+
geom_point(data=neut %>% filter(hours==48), aes(x=variant, y=ID50), shape=21, fill='white', position=position_jitter(0.3, seed=42))+
    geom_errorbar(data=neut_stat, aes(x=variant, y=gmean_titer, ymin=gcilow, ymax=gcihigh), width=.5)+
theme_classic()+
theme(axis.title=element_blank(), axis.text.x=element_blank())
dev.off()

neut_stat <- neut %>% group_by(group, variant, hours) %>% summarise(GMT=exp(mean(log(ID50))))

pdf("plots/Fig5/Fig5d-neut_line.pdf", width=4.5, height=2.3)
ggplot(neut_stat, aes(hours, GMT))+
    scale_y_log10(limits=c(100, 10000))+scale_color_manual(values=c("#1b9e77", "#7570b3", "#d95f02"))+
geom_line(data=neut, aes(x=hours, y=ID50, group=interaction(variant, mouse), color=variant), alpha=.1)+
geom_line(aes(color=variant))+
geom_point(aes(color=variant), shape=21, size=2, fill="white")+
# geom_hline(yintercept=20, linetype='dashed')+
geom_point(data=neut, aes(x=hours, y=ID50, shape=variant, color=variant), size=1, alpha=.6)+
theme_classic()
dev.off()

# Fig.5f corr
data_merge <- merge(data, neut %>% pivot_wider(id_cols=c(hours, sex, mouse), names_from=variant, values_from=ID50), by=c('hours', 'sex', 'mouse'), all=T) %>% na.omit()
pdf("plots/Fig5/Fig5f-conc_neut_corr.pdf", width=3, height=3)
Min <- min(neut$ID50, na.rm = T)
Max <- max(neut$ID50, na.rm = T)
ggplot(data_merge, aes(conc*0.146, XBB.1.5))+geom_point(aes(color=hours), size=2, alpha=.6)+scale_colour_distiller(palette = "Spectral")+scale_x_log10()+scale_y_log10(limits=c(Min, Max))+theme_classic()+theme(aspect.ratio = 1.0,axis.title=element_blank())
ggplot(data_merge, aes(conc*0.146, HK.3.1))+geom_point(aes(color=hours), size=2, alpha=.6)+scale_colour_distiller(palette = "Spectral")+scale_x_log10()+scale_y_log10(limits=c(Min, Max))+theme_classic()+theme(aspect.ratio = 1.0,axis.title=element_blank())
ggplot(data_merge, aes(conc*0.146, JN.1))+geom_point(aes(color=hours), size=2, alpha=.6)+scale_colour_distiller(palette = "Spectral")+scale_x_log10()+scale_y_log10(limits=c(Min, Max))+theme_classic()+theme(aspect.ratio = 1.0,axis.title=element_blank())
dev.off()

neut <- read.csv("../source_data/Fig5-mice-neutralization.csv", check.names=F, colClasses = c("numeric","character","character","numeric","character","character"))
neut$group <- paste(neut$type, neut$sex, sep="-")
neut$variant <- factor(neut$variant, levels=c("XBB.1.5", "HK.3.1", "JN.1"))
pdf("plots/Fig5/Ext-Tg32-conc_neut_corr.pdf", width=3, height=3)
data_merge <- merge(data%>%select(!group), neut %>% pivot_wider(id_cols=c(hours, sex, mouse, group), names_from=variant, values_from=ID50), by=c('hours', 'sex', 'mouse'), all=T) %>% na.omit()
Min <- min(neut$ID50, na.rm = T)
Max <- max(neut$ID50, na.rm = T)
ggplot(data_merge, aes(conc*0.146, XBB.1.5))+geom_point(aes(color=group, shape=group), size=2, alpha=.6)+scale_color_manual(values=c("#1b9e77", "#7570b3", "#d95f02"))+scale_x_log10()+scale_y_log10(limits=c(Min, Max))+theme_classic()+theme(aspect.ratio = 1.0,axis.title=element_blank())+ggtitle('XBB.1.5')
ggplot(data_merge, aes(conc*0.146, HK.3.1))+geom_point(aes(color=group, shape=group), size=2, alpha=.6)+scale_color_manual(values=c("#1b9e77", "#7570b3", "#d95f02"))+scale_x_log10()+scale_y_log10(limits=c(Min, Max))+theme_classic()+theme(aspect.ratio = 1.0,axis.title=element_blank())+ggtitle('HK.3.1')
ggplot(data_merge, aes(conc*0.146, JN.1))+geom_point(aes(color=group, shape=group), size=2, alpha=.6)+scale_color_manual(values=c("#1b9e77", "#7570b3", "#d95f02"))+scale_x_log10()+scale_y_log10(limits=c(Min, Max))+theme_classic()+theme(aspect.ratio = 1.0,axis.title=element_blank())+ggtitle('JN.1')
dev.off()
