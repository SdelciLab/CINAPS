library(readxl)
library(cowplot)
library(dplyr)
library(gganimate)
library(ggbreak)
library(ggforce)
library(ggpattern)
library(ggplot2)
library(ggpointdensity)
library(ggpubr)
library(ggridges)
library(ggsankey)
library(grid)
library(gridExtra)
library(magick)
library(MASS)
library(ptinpoly)
library(readxl)
library(rstatix)
library(tidyverse)
library(viridis)
library(zoo)


jco = c("#0073C2FF","#EFC000FF")
jco2 = c("#EFC000FF", "#0073C2FF")

# ATP sensor experiments ----------------------------------

combined_atp <- read_xlsx("supplementary_data_1.xlsx", sheet = "atp_sensor") 

# Nuclear ATP DF 
# Filter main data frame for nuclear ATP only.
atpsensor_nuc <- combined_atp %>%  
    dplyr::filter(Staining == "nuc") %>%  
    dplyr::select(-Mean, -Slice) %>%
    unique()

# Count cell numbers for each
atpsensor_nuc_n <- atpsensor_nuc %>% 
    group_by(Condition, Treatment, Cell_line) %>% 
    tally() 

# Add cell numbers to DF
atpsensor_nuc <- atpsensor_nuc %>% left_join(atpsensor_nuc_n) %>% 
    mutate(labelwithN=paste0(Treatment,"\nn=",n))

#Normalizing values
lm_meanATP <- lm(meanX ~ as.factor(Date), data = atpsensor_nuc)
noeffrep_ATP <- residuals(lm_meanATP) + coef(lm_meanATP)[1]

atpsensor_nuc_corrected <- cbind(atpsensor_nuc ,  
                                 meanATP_corrected = noeffrep_ATP) %>% 
    as_tibble()

# HeLa Nuclear ATP 
#Stat tests
#Batch correction not used because it looks like there is not a lot of variability across replicates
stat.test.nuc.hela <- compare_means(data=atpsensor_nuc_corrected %>% dplyr::filter(Cell_line=='hela'), 
                                    meanX~labelwithN,
                                    # group.by = "Treatment",
                                    method='wilcox.test') %>% 
    slice(1,7,8,10,15)


ggplot(atpsensor_nuc_corrected %>% 
           dplyr::filter(Cell_line == "hela",
                         !Date %in% c() ),
       aes(x=factor(labelwithN,  
                    level=c("control\nn=78",
                            "control\nn=125",
                            "latrunculin a (500 nm)\nn=26",
                            "latrunculin a (500 nm)\nn=35",
                            "oligomycin\nn=46",
                            "oligomycin\nn=53")),
           y=meanX))+
    geom_boxplot(aes(fill=Condition),
                 colour='black',
                 width=0.8, 
                 outlier.shape=NA, notch = F)+
    scale_fill_manual(values=jco)+
    scale_color_manual(values=jco)+
    theme_bw()+
    theme(panel.grid = element_blank(),
          axis.title.x = element_blank(),
          legend.position = "none")+
    ylab("Nuclear ATP")+
    xlab("Treatment")+
    geom_jitter(shape=21,size=3, aes(fill=Condition), width = 0.1)+
    stat_pvalue_manual(stat.test.nuc.hela,
                       label = "p.signif",
                       y.position=c(7,7.5,8,7,7))+
    ggtitle("HeLa")



# U2OS Nuclear ATP 
#Stats
#Batch correction used in case of u2os because independent experiments show high variability unlike hela
stat.test.nuc.u2os <- compare_means(data=atpsensor_nuc_corrected %>% dplyr::filter(Cell_line=='u2os'), 
                                    meanATP_corrected~labelwithN,
                                    # group.by = "Treatment",
                                    method='wilcox.test') %>% 
    slice(1,5,6)

#Plot U2OS
ggplot(atpsensor_nuc_corrected %>% 
           dplyr::filter(Cell_line == "u2os",
                         !Date %in% c() ),
       aes(x=factor(labelwithN,  
                    level=c("control\nn=43",
                            "control\nn=66",
                            "oligomycin\nn=50",
                            "oligomycin\nn=79")),
           y=meanATP_corrected))+
    geom_boxplot(aes(fill=Condition),
                 colour='black',
                 width=0.8, 
                 outlier.shape=NA, notch = F)+
    scale_fill_manual(values=jco)+
    scale_color_manual(values=jco)+
    theme_bw()+
    theme(panel.grid = element_blank(),
          axis.title.x = element_blank(),
          legend.position = "none")+
    ylab("Nuclear ATP")+
    xlab("Treatment")+
    geom_jitter(shape=21,size=2, aes(fill=Condition), width = 0.1)+
    stat_pvalue_manual(stat.test.nuc.u2os,
                       label = "p.signif",
                       y.position=c(6,6.2,6))+
    ggtitle("U2OS")


# Mitochondrial ATP -----------------------------------------------
atpsensor_mito <- combined_atp %>%  
    dplyr::filter(Staining %in%  c("nam", "mito_nam", "mito", "mito_all")) %>% 
    mutate(Category = case_when(Staining %in% c("nam", "mito_nam") ~ "NAM",
                                Staining %in% c("mito", "mito_all") ~ "Total Mito")) %>% 
    dplyr::select(-Mean, -Slice) %>% 
    unique()

atpsensor_mito_n <- atpsensor_mito %>% group_by(Condition, Treatment, Category) %>% 
    tally()

atpsensor_mito <- atpsensor_mito %>% left_join(atpsensor_mito_n) %>% 
    mutate(labelwithN = paste0(Category, "\n", Condition,"\nn=",n))

stat.test.mito <- compare_means(data = atpsensor_mito, 
                                meanX ~ labelwithN,
                                method = 'wilcox.test') %>% 
    dplyr::filter(row_number() %in% c(2,5))

ggplot(atpsensor_mito , aes(x=factor(labelwithN,  
                                     level=c("Total Mito\nsusp\nn=11",
                                             "Total Mito\nconf\nn=47",
                                             "NAM\nsusp\nn=11",
                                             "NAM\nconf\nn=44")),
                            y=meanX))+
    geom_boxplot(aes(fill=Condition),
                 colour='black',
                 width=0.8, 
                 outlier.shape=NA, notch = F)+
    geom_jitter(shape=21,size=3, aes(fill=Condition), width = 0.1)+
    stat_pvalue_manual(stat.test.mito,  
                       label = "p.signif", 
                       y.position=5)+
    scale_fill_manual(values=jco)+
    scale_color_manual(values=jco)+
    theme_bw()+
    theme(panel.grid = element_blank(),axis.title.x = element_blank(), legend.position="none")+
    ylab("Mitochondrial ATP levels")

# Time and Nuclear ATP levels ------------------------------
#Plotting only done for control confined hela cells stained for nuclear ATP

nuc_df_time_stamp <- read_xlsx("supplementary_data_1.xlsx", sheet = "time_nucATP")  %>%
    mutate(time = as.numeric(format(as.POSIXct(Time_rel), format = "%M")))

ggplot(nuc_df_time_stamp, 
       aes(x=time, y=as.numeric(meanATP_corrected)))+
    geom_point()+
    # geom_line()+
    geom_smooth(method = lm)+
    theme_bw()+
    # facet_wrap(~Date, scales = "free")+
    theme(panel.grid = element_blank())+
    # stat_regline_equation(label.y = 1.5, aes(label = ..eq.label..)) +
    stat_regline_equation(label.y = 3.1,label.x = 20, aes(label = ..rr.label..))+
    xlab("Time (min)")+
    ylab("Nuclear ATP")





# NAM ------------------------------
mitoloc <- read_xlsx("supplementary_data_1.xlsx", sheet = "mitoloc") 

# Compute means
meanNAM <- mitoloc %>%
    dplyr::filter(Treatment=="Control",
                  Concentration_stain %in% c("500 nM", "200 nM", "100 nM", "50 nM")) %>% 
    group_by(Treatment, Condition, labs, Exp, Cell_line, CellID, Concentration_stain) %>% 
    summarise(meanRID = mean(RawIntDen)) %>% 
    ungroup() %>% 
    group_by(Treatment, Condition, Cell_line) %>% 
    add_tally() %>% 
    ungroup() %>% 
    mutate(labelwithN=paste0(Condition,"\nn=",n),
           across(Condition, factor, levels=c("Susp.", "Conf.")))

stat.meanNAM <- compare_means(data = meanNAM,
                              meanRID ~ labelwithN,
                              method = "wilcox.test") %>% 
    slice(4,11,17,22) %>% 
    cbind(Cell_line = c("HeLa","MiaPaca2","MDA-MB-231","U2OS"),
          ypos = c(3e+6,8e+5,6e+5,3e+6))


#plot all controls across cell lines
ggplot(meanNAM,aes(x=factor(labelwithN,
                            level=c("Susp.\nn=167","Conf.\nn=300",
                                    "Susp.\nn=29","Conf.\nn=23",
                                    "Susp.\nn=19","Conf.\nn=24",
                                    "Susp.\nn=38","Conf.\nn=53")),
                   y=meanRID))+
    geom_boxplot(aes(fill=Condition),lwd=0.8, fatten=0.8, outlier.shape = NA)+
    geom_jitter(shape=21,size=3, aes(fill=Condition), width = 0.1, alpha=1)+
    scale_fill_manual(values=c(jco[2],jco[1]))+
    theme_bw()+
    theme(
        legend.position = "none",
        axis.title.x=element_blank(),
        axis.text.x=element_text(size=9))+
    ylab("NAM\n(Total mitotracker intensity)")+
    facet_wrap(~Cell_line, scales="free", nrow=1)+
    stat_pvalue_manual(stat.meanNAM,
                       label = "p.signif",
                       y.position= "ypos")


# Only experiments with BAPTA, LAT (500) CK666 SMIFh2 (100) and 50 nM Mitotracker intensity 
meanNAM_drugs <- mitoloc %>% 
    dplyr::filter(Cell_line == "HeLa", 
                  Exp %in% c("Exp008",
                             "Exp009",
                             "Exp018",
                             "Exp019",
                             "Exp020",
                             "Exp021",
                             "Exp022",
                             "Exp023"),
                  Treatment %in% c("Control", 
                                   "BAPTA-AM", 
                                   "Latrunculin A (500 nM)", 
                                   "CK666", 
                                   "SMIFH2"
                  ),
                  Concentration_stain %in% c("500 nM", "200 nM", "100 nM", "50 nM"))

meanNAM_drugs_separate <- meanNAM_drugs %>% 
    group_by(Treatment, Condition, labs, Exp, CellID, Concentration_drug) %>% 
    summarise(meanRID = mean(RawIntDen)) %>% 
    ungroup() %>% 
    group_by(Treatment, Condition, Concentration_drug) %>% 
    add_tally() %>% 
    ungroup() %>% 
    mutate(labelwithN = paste0(Treatment, "\n", Condition,"\nn=",n),
           across(Condition, factor, levels=c("Susp.", "Conf.")))


# Perform statistics for both the categories
stat.meanNAM.separate <- compare_means(data = meanNAM_drugs_separate,
                                       meanRID ~ labelwithN,
                                       method = "wilcox.test") %>%
    slice(1,18,31,40,45) %>% 
    cbind(Treatment = c("BAPTA-AM", "CK666", "Control", "Latrunculin A (500 nM)", "SMIFH2"),
          y.position = c(1.5e+6, 1.0e+6, 1.0e+6, 4.2e+5, 5.0e+6))



# Plot the NAM for all various drug treatments, but independently, with their own respective controls

ggplot(meanNAM_drugs_separate %>% mutate(Filt = paste0(Treatment,meanRID)) %>% 
           # Not plotting extreme outlier for nicer plot. Note, they were used for stats.
           dplyr::filter(Filt != "BAPTA-AM2311570.69230769",
                         Filt != "Latrunculin A (500 nM)852615.75"),
       aes(
           x=factor(
               labelwithN,level=c(
                   "Control\nSusp.\nn=62",
                   "Control\nConf.\nn=124",
                   "BAPTA-AM\nSusp.\nn=25",
                   "BAPTA-AM\nConf.\nn=47",
                   "Latrunculin A (500 nM)\nSusp.\nn=26",
                   "Latrunculin A (500 nM)\nConf.\nn=36",
                   "CK666\nSusp.\nn=32",
                   "CK666\nConf.\nn=64",
                   "SMIFH2\nSusp.\nn=34",
                   "SMIFH2\nConf.\nn=65"
               )),
           y=meanRID
           
       ))+
    geom_boxplot(lwd=0.8, fatten=0.8, outlier.shape = NA, aes(fill=Condition))+
    scale_fill_manual(values=c(jco[2],jco[1]))+
    theme_bw()+
    theme(
        legend.position = "none",
        axis.title.x=element_blank(),
        axis.text.x=element_text(size=12))+
    geom_jitter(shape=21,size=3, aes(fill=Condition), width = 0.1)+
    facet_wrap(~factor(Treatment, levels = c("Control",
                                             "BAPTA-AM",
                                             "Latrunculin A (500 nM)",
                                             "Cytochalasin D",
                                             "CK666",
                                             "SMIFH2")), nrow=2, scales="free")+
    ylab("NAM\n(Total mitotracker intensity)")+
    stat_pvalue_manual(stat.meanNAM.separate,
                       label = "p.signif")



#Normalizing values to present data with single control suspension
meanNAM_drugs_control_susp = dplyr::filter(meanNAM_drugs_separate, Treatment=="Control", Condition=="Susp.") %>% 
    mutate(meanRID_mean = mean(meanRID),sdRID = sd(meanRID),
           newRID = ((meanRID - meanRID_mean)/sdRID)+1
    )
meanNAM_drugs_control_conf = dplyr::filter(meanNAM_drugs_separate, Treatment=="Control", Condition=="Conf.") %>% 
    mutate(
        newRID = 
            ((meanRID-(meanNAM_drugs_control_susp$meanRID_mean %>% unique()))/(meanNAM_drugs_control_susp$sdRID %>% unique())+1)
    )

meanNAM_drugs_bapta_susp = dplyr::filter(meanNAM_drugs_separate, Treatment=="BAPTA-AM", Condition=="Susp.") %>% 
    mutate(meanRID_mean = mean(meanRID),sdRID = sd(meanRID),
           newRID = ((meanRID - meanRID_mean)/sdRID)+1
    )
meanNAM_drugs_bapta_conf = dplyr::filter(meanNAM_drugs_separate, Treatment=="BAPTA-AM", Condition=="Conf.") %>% 
    mutate(
        newRID = 
            ((meanRID-(meanNAM_drugs_bapta_susp$meanRID_mean %>% unique()))/(meanNAM_drugs_bapta_susp$sdRID %>% unique()))+1
    )

meanNAM_drugs_lat_susp = dplyr::filter(meanNAM_drugs_separate, Treatment=="Latrunculin A (500 nM)", Condition=="Susp.") %>% 
    mutate(meanRID_mean = mean(meanRID),sdRID = sd(meanRID),
           newRID = ((meanRID - meanRID_mean)/sdRID)+1
    )
meanNAM_drugs_lat_conf = dplyr::filter(meanNAM_drugs_separate, Treatment=="Latrunculin A (500 nM)", Condition=="Conf.") %>% 
    mutate(
        newRID = 
            ((meanRID-(meanNAM_drugs_lat_susp$meanRID_mean %>% unique()))/(meanNAM_drugs_lat_susp$sdRID %>% unique()))+1
    )

meanNAM_drugs_ck_susp = dplyr::filter(meanNAM_drugs_separate, Treatment=="CK666", Condition=="Susp.") %>% 
    mutate(meanRID_mean = mean(meanRID),sdRID = sd(meanRID),
           newRID = ((meanRID - meanRID_mean)/sdRID)+1
    )
meanNAM_drugs_ck_conf = dplyr::filter(meanNAM_drugs_separate, Treatment=="CK666", Condition=="Conf.") %>% 
    mutate(
        newRID = 
            ((meanRID-(meanNAM_drugs_ck_susp$meanRID_mean %>% unique()))/(meanNAM_drugs_ck_susp$sdRID %>% unique()))+1
    )

meanNAM_drugs_smifh2_susp_100 = dplyr::filter(meanNAM_drugs_separate, Treatment=="SMIFH2", Condition=="Susp.", labelwithN == "SMIFH2\nSusp.\nn=34") %>% 
    mutate(meanRID_mean = mean(meanRID), sdRID = sd(meanRID),
           newRID = ((meanRID - meanRID_mean)/sdRID)+1
    )
meanNAM_drugs_smifh2_conf_100 = dplyr::filter(meanNAM_drugs_separate, Treatment=="SMIFH2", Condition=="Conf.", labelwithN == "SMIFH2\nConf.\nn=65") %>% 
    mutate(
        newRID = 
            ((meanRID-(meanNAM_drugs_smifh2_susp_100$meanRID_mean %>% unique()))/(meanNAM_drugs_smifh2_susp_100$sdRID %>% unique()))+1
    )


meanNAM_drugs_corrected = rbind(meanNAM_drugs_control_susp %>% dplyr::select(-meanRID_mean, -sdRID),
                                meanNAM_drugs_control_conf,
                                meanNAM_drugs_bapta_susp %>% dplyr::select(-meanRID_mean, -sdRID),
                                meanNAM_drugs_bapta_conf,
                                meanNAM_drugs_lat_susp %>% dplyr::select(-meanRID_mean, -sdRID),
                                meanNAM_drugs_lat_conf,
                                meanNAM_drugs_ck_susp %>% dplyr::select(-meanRID_mean, -sdRID),
                                meanNAM_drugs_ck_conf,
                                meanNAM_drugs_smifh2_susp_100 %>% dplyr::select(-meanRID_mean, -sdRID),
                                meanNAM_drugs_smifh2_conf_100) 

meanNAM_drugs_corrected = meanNAM_drugs_corrected %>% 
    mutate(
        Treatment = case_when(grepl("Latrunculin",labelwithN)==TRUE ~ "LatA",
                              !grepl("Latrunculin",labelwithN)==TRUE ~ Treatment),
        labelwithN_comb=case_when(grepl("Susp.",labelwithN)==TRUE ~ "Susp.\nn=179", #Sum of all suspensions in this plot
                                  grepl("Conf.",labelwithN)==TRUE ~ labelwithN),
        across(Condition, factor, levels=c("Susp.", "Conf."))
    )



# Plot single Susp ctrl - corrected
stat.meanNAM.corrected <- compare_means(data = meanNAM_drugs_corrected,
                                        newRID ~ labelwithN_comb,
                                        method = "wilcox.test") %>% 
    slice(1,6:9,14)


# Plot all drugs standardised with controls, and hence plot single control.
ggplot(meanNAM_drugs_corrected, aes(x=factor(labelwithN_comb, 
                                             level=c("Susp.\nn=179",
                                                     "Control\nConf.\nn=124",
                                                     "BAPTA-AM\nConf.\nn=47",
                                                     "Latrunculin A (500 nM)\nConf.\nn=36",
                                                     "CK666\nConf.\nn=64",
                                                     "SMIFH2\nConf.\nn=65")),
                                    y=newRID))+
    geom_boxplot(lwd=0.8, fatten=0.8, outlier.shape = NA, aes(fill=Condition))+
    geom_jitter(shape=21,size=2, aes(fill=Condition), width = 0.1)+
    scale_fill_manual(values=c(jco[2],jco[1]))+
    theme_bw()+
    theme(
        panel.grid = element_blank(),
        legend.position = "none",
        axis.title.x=element_blank(),
        axis.text.x=element_text(size=12))+
    ylab("NAM\n(Total mitotracker intensity)")+
    stat_pvalue_manual(stat.meanNAM.corrected,  
                       label = "p.signif", 
                       y.position=c(20,21,22,23,19,20.5))

# Time and NAM levels ------------------------------
time_NAM_df <- meanNAM_drugs_corrected %>% 
    dplyr::filter(Treatment == "Control", Exp!= "Exp019")

time_stamp <- read_xlsx("supplementary_data_1.xlsx", sheet = "time_nam")  %>% 
    dplyr::select(Treatment, Condition, labs, Exp, CellID, Time_rel)

time_NAM_df_stamp <- left_join(time_NAM_df, time_stamp) %>% 
    mutate(time = as.numeric(format(as.POSIXct(Time_rel), format = "%M"))) %>% 
    dplyr::filter(Condition == "Conf.")

ggplot(time_NAM_df_stamp, 
       aes(x=time, y=newRID))+
    geom_point()+
    geom_smooth(method = lm, )+
    theme_bw()+
    facet_wrap(~Condition, scales="free")+
    theme(panel.grid = element_blank())+
    stat_regline_equation(label.y = 12,label.x = 5, aes(label = ..rr.label..))+
    xlab("Time (h)")+
    ylab("NAM\n(Total mitotracker intensity)")



# Coefficient of Variation ----------------------------------------------------------------------

cv_data <- read_xlsx("supplementary_data_1.xlsx", sheet = "CoV") 

cv_process <- cv_data %>% mutate(coeffvar = StdDev/Mean,
                                 across(Condition, factor, level=c("Susp.", "Conf.")),
                                 across(Treatment, factor, level=c("Control","Oligomycin"))
) %>%
    mutate(Treat_Cond = paste0(Treatment, Condition)) %>%  
    group_by(Treatment, Condition, Bin) %>%
    ungroup()

stat.cov <- data.frame()

for(i in 1:8){
    stat <- cv_process %>% dplyr::filter(Bin==i) %>% 
        compare_means(data =.,  coeffvar ~ Treat_Cond, method="t.test")  %>% 
        slice(1,2,5,6) %>% 
        mutate(Bin = i)
    
    stat.cov <- rbind(stat.cov, stat)
}

print(stat.cov)

cv_process %>% ggplot(aes(x = Bin, y = coeffvar,
                          fill=Condition, 
                          colour=Condition,
                          linetype=Treatment))+
    geom_smooth()+
    theme_bw()+
    scale_fill_manual(values= jco2)+
    scale_colour_manual(values = jco2)+
    ylab("Coefficient of Variation\n(StdDev/Mean)")


# 53BP1 under confinement ------------------------------------------------
damage_data <- read_xlsx("supplementary_data_1.xlsx", sheet = "dna_dam_conf") %>% 
    dplyr::select(Repetition, Condition, Treatment, Area, CellID) %>% 
    dplyr::filter(Area < 3) %>%
    group_by(Repetition, Condition, Treatment, CellID) %>% 
    add_tally() %>% 
    group_by(Repetition, Condition, Treatment, CellID) %>% 
    mutate(mean_Area = mean(Area)) %>% 
    ungroup() %>% 
    dplyr::select(-Area,-CellID) %>% 
    unique() 

#Adding the zeros
df1<-data.frame(Repetition=c("Artificial"),
                Condition=c("Susp."),
                Treatment=c("Control"),
                n=c(0),
                mean_Area=c(0))
df1<-df1[rep(seq_len(nrow(df1)), each = 13), ]



df2<-data.frame(Repetition=c("Artificial"),
                Condition=c("Conf."),
                Treatment=c("Control"),
                n=c(0),
                mean_Area=c(0))
df2<-df2[rep(seq_len(nrow(df1)), each = 6), ]



df3<-data.frame(Repetition=c("Artificial"),
                Condition=c("Susp."),
                Treatment=c("Oligomycin"),
                n=c(0),
                mean_Area=c(0))
df3<-df3[rep(seq_len(nrow(df3)), each = 9), ]



df4<-data.frame(Repetition=c("Artificial"),
                Condition=c("Conf."),
                Treatment=c("Oligomycin"),
                n=c(0),
                mean_Area=c(0))
df4<-df4[rep(seq_len(nrow(df4)), each = 16), ]

damage_data <- damage_data %>% 
    rbind(df1, df2, df3, df4) %>% 
    na.omit()

damage_data_n <- damage_data %>% 
    group_by(Condition, Treatment) %>% 
    tally() %>% 
    set_names("Condition", "Treatment", "count")

damage_data <- left_join(damage_data, damage_data_n) %>% 
    mutate(labelwithN = paste0(Treatment, "\n", Condition,"\nn=", count),
           across(labelwithN, factor, levels=c("Control\nSusp.\nn=34",
                                               "Control\nConf.\nn=55",
                                               "Oligomycin\nSusp.\nn=37",
                                               "Oligomycin\nConf.\nn=65")))

stat.test.damage <- compare_means(data=damage_data, 
                                  n~labelwithN,
                                  method='wilcox.test') %>% 
    dplyr::filter(row_number() %in% c(1,2,5,6))


ggplot(damage_data, aes(x=labelwithN, y=n))+
    geom_boxplot(aes(fill=Condition),
                 colour='black',
                 width=0.8, 
                 outlier.shape=NA)+
    geom_jitter(shape=21,size=3, aes(fill=Condition), width = 0.1)+
    theme_bw()+
    scale_fill_manual(values=jco)+
    scale_color_manual(values=jco)+
    stat_pvalue_manual(stat.test.damage,  
                       label = "p.signif", 
                       y.position=c(42,44,46,42))+
    theme(panel.grid = element_blank(),
          axis.title.x = element_blank(),
          legend.position="none")

# 53BP1 post-confinement -------------------------------------

post_conf_53bp1 <- read.csv("supplementary_data_2.csv")

fp_53bp1 = post_conf_53bp1 %>% 
    mutate(dif488 = nuc_int_488 - ring_int_488,
           difTurq = nuc_int_Turq - ring_int_Turq,
           dif546 = nuc_int_546 - ring_int_546,
           dif488_comp = dif488-(0.571*difTurq-26.59)) %>%
    dplyr::filter(dif488_comp > 0,
                  difTurq > 0,
                  dif546 > 0) %>%
    group_by(Replicate, Treatments) %>%
    mutate(difTurq = log10(difTurq),
           dif488_comp = log10(dif488_comp),
           dif546 = log10(dif546),
           zdifTurq = scale(difTurq),
           zdif488_comp = scale(dif488_comp),
           zdif546 = scale(dif546),
    ) %>%
    filter(between(zdifTurq ,-3,+3),
           between(zdif488_comp ,-3,+3),
           between(zdif546 ,-3,+3)) %>%
    mutate(
        n488 = (dif488_comp-min(dif488_comp))/(max(dif488_comp)-min(dif488_comp)),
        nTurq = (difTurq-min(difTurq))/(max(difTurq)-min(difTurq)),
        n546 = (dif546-min(dif546))/(max(dif546)-min(dif546))
    ) %>%
    ungroup() %>%
    dplyr::filter(n488 >= 0,
                  nTurq >= 0,
                  n546 >= 0) %>% 
    mutate(across(Treatments, factor, levels = c("Control\nUnconfined",
                                                 "Control\nConfined",
                                                 "Oligomycin\nUnconfined",
                                                 "Oligomycin\nConfined")))


df_1_6 <- fp_53bp1 %>% 
    dplyr::filter(Replicate %in% c("Rep2","Rep3","Rep4"),
                  spot_count < 40,
                  Timepoint %in% c(3:7)) %>% 
    mutate(Timepoint_addon = case_when(Treatments == "Control\nUnconfined" ~ 0,
                                       Treatments == "Control\nConfined" ~ 6,
                                       Treatments == "Oligomycin\nUnconfined" ~ 12,
                                       Treatments == "Oligomycin\nConfined" ~ 18),
           Time = Timepoint+Timepoint_addon)

stats <- data.frame()
for(i in df_1_6$Timepoint %>% unique()){
    stat.test.damage <- compare_means(data = df_1_6 %>% dplyr::filter(Timepoint == i),
                                      spot_count ~ Treatments,
                                      method='wilcox.test') %>% 
        slice(1,6)
    stats <- stats %>% rbind(stat.test.damage %>% mutate(Timepoint = i))
}

ggplot(df_1_6, aes(x = factor(Time), y=spot_count, fill=Treatments))+
    geom_boxplot(outlier.shape=NA)+
    geom_jitter(height=0.5, width=0.2, alpha=0.02, size=1)+
    ylab("53BP1 foci")+
    theme_bw()+
    ylim(0,20)+
    scale_fill_manual(values=c(jco2,jco2))

# Cell cycle post confinement FUCCI -------------------------------------------------------------------

fuccidata <- read.csv("supplementary_data_3.csv") 

# Intensity corrections and FUCCI calculations 
fp = fuccidata %>%
    mutate(dif488 = nuc_int_488 - ring_int_488,
           difTurq = nuc_int_Turq - ring_int_Turq,
           dif546 = nuc_int_546 - ring_int_546,
           dif488_comp = dif488-(0.571*difTurq-26.59)) %>%
    dplyr::filter(dif488_comp > 0,
                  difTurq > 0,
                  dif546 > 0) %>%
    group_by(Replicate, Treatments) %>%
    mutate(difTurq = log10(difTurq),
           dif488_comp = log10(dif488_comp),
           dif546 = log10(dif546),
           zdifTurq = scale(difTurq),
           zdif488_comp = scale(dif488_comp),
           zdif546 = scale(dif546),
    ) %>%
    filter(between(zdifTurq ,-3,+3),
           between(zdif488_comp ,-3,+3),
           between(zdif546 ,-3,+3)) %>%
    mutate(
        n488 = (dif488_comp-min(dif488_comp))/(max(dif488_comp)-min(dif488_comp)),
        nTurq = (difTurq-min(difTurq))/(max(difTurq)-min(difTurq)),
        n546 = (dif546-min(dif546))/(max(dif546)-min(dif546))
    ) %>%
    ungroup() %>%
    dplyr::filter(n488 >= 0,
                  nTurq >= 0,
                  n546 >= 0) 

# Proliferation
proliferation_df <- fp %>% 
    dplyr::filter(Timepoint %in% c(1:36)) %>% 
    group_by(Replicate, Treatments, Timepoint) %>% 
    unique() %>% 
    tally() %>% 
    group_by(Replicate, Treatments) %>% 
    mutate(n_norm = n/n[Timepoint==1]) %>% 
    ungroup() %>% 
    mutate(across(Treatments, levels = c("Control\nUnconfined",
                                         "Control\nConfined",
                                         "Oligomycin\nUnconfined",
                                         "Oligomycin\nConfined")))
# Wilcox Test on averaged data
test_data = proliferation_df %>% 
    group_by(Treatments, Timepoint) %>% 
    mutate(SD=sd(n_norm),
           Mean = mean(n_norm)) %>% 
    ungroup() %>% 
    dplyr::select(Timepoint, Treatments, Mean) %>% 
    unique() %>% 
    pivot_wider(names_from = Treatments, id_cols=Timepoint, values_from = Mean)

stat.ctrl <- wilcox.test(test_data$`Control\nUnconfined`,
                         test_data$`Control\nConfined`)
# 
stat.oligo <- wilcox.test(test_data$`Oligomycin\nUnconfined`,
                          test_data$`Oligomycin\nConfined`)

# Cell Proliferation
ggplot(proliferation_df, aes(x=Timepoint, y=n_norm))+
    geom_smooth(aes(color=Treatments), size=1)+
    theme_bw() +
    ggtitle("Combined plot Rep3, Rep4, Rep6: NO FILTER") +
    annotate(geom="text", x=1, y=2, col="black",
             label=paste("Mann Whitney U test",
                         "\nCtrl (Conf. vs Unconf) p = ",signif(stat.ctrl$p.value,3),
                         "\nOligo (Conf. vs Unconf) p = ",signif(stat.oligo$p.value,3)),
             hjust = 0)+
    scale_y_continuous(breaks=c(1,1.5,2))



# Establish cell cycle gating
shiny_gates = read_xlsx("supplementary_data_1.xlsx", sheet="fucci_gates") %>%
    mutate(Treatments = case_when(Treat=='UC'~"Control\nUnconfined",
                                  Treat=='CC'~"Control\nConfined",
                                  Treat=='UT'~"Oligomycin\nUnconfined",
                                  Treat=='CT'~"Oligomycin\nConfined",
                                  TRUE~as.character(Treat)),
           phase=name) %>%
    dplyr::select(-Treat,-name)

phased_fucci <- data.frame()

for(rep in fp$Replicate %>% unique()){
    for(treat in fp$Treatments %>% unique()){
        
        gates = dplyr::filter(shiny_gates, Replicate == rep, Treatments == treat)
        
        phased_temp <- fp %>%
            dplyr::filter(Replicate == rep, Treatments == treat) %>%
            mutate(G1_check = pip2d(as.matrix(gates %>% dplyr::filter(phase=="G1") %>% dplyr::select(1,2)), as.matrix(data.frame(nTurq, n488))),
                   S_check = pip2d(as.matrix(gates %>% dplyr::filter(phase=="S") %>% dplyr::select(1,2)), as.matrix(data.frame(nTurq, n488))),
                   G2_check = pip2d(as.matrix(gates %>% dplyr::filter(phase=="G2") %>% dplyr::select(1,2)), as.matrix(data.frame(nTurq, n488))),
                   M_check = pip2d(as.matrix(gates %>% dplyr::filter(phase=="M") %>% dplyr::select(1,2)), as.matrix(data.frame(nTurq, n488))),
                   Sen_check = pip2d(as.matrix(gates %>% dplyr::filter(phase=="Sen") %>% dplyr::select(1,2)), as.matrix(data.frame(nTurq, n488))),
                   phase = case_when(G1_check==1~"G1",
                                     S_check==1~"S",
                                     G2_check==1~"G2",
                                     M_check==1~"M",
                                     Sen_check==1~"Sen"
                   ))
        
        phased_temp$phase[is.na(phased_temp$phase)] <- "None"
        
        phased_fucci <- rbind(phased_fucci, phased_temp)
    }
}

startPhase <- phased_fucci %>%
    dplyr::filter(Timepoint == 1) %>%
    dplyr::select(cellID, phase) %>%
    unique()


# Phase-wise_abundance 

#Assign time point based on one replicate. All replicates have same time points.
tp <- fp$Timepoint[fp$Replicate == "Rep3"] %>% unique() %>% sort()

phased_fp = data.frame()
for(treat in c("Control\nUnconfined",
               "Control\nConfined",
               "Oligomycin\nUnconfined",
               "Oligomycin\nConfined")){
    
    for(rep in c("Rep3",
                 "Rep4",
                 "Rep6")){
        
        gates <- shiny_gates %>% 
            dplyr::filter(Treatments == treat,
                          Replicate == rep)
        
        temp <- fp %>% 
            dplyr::filter(Treatments == treat,
                          Replicate == rep,
                          Timepoint %in% tp) %>%
            arrange(n546) %>% 
            mutate(across(Treatments, factor, levels=c("Control\nUnconfined",
                                                       "Control\nConfined",
                                                       "Oligomycin\nUnconfined",
                                                       "Oligomycin\nConfined")),
                   G1_check = pip2d(as.matrix(gates %>% dplyr::filter(phase=="G1") %>% dplyr::select(1,2)), as.matrix(data.frame(nTurq, n488))),
                   S_check = pip2d(as.matrix(gates %>% dplyr::filter(phase=="S") %>% dplyr::select(1,2)), as.matrix(data.frame(nTurq, n488))),
                   G2_check = pip2d(as.matrix(gates %>% dplyr::filter(phase=="G2") %>% dplyr::select(1,2)), as.matrix(data.frame(nTurq, n488))),
                   M_check = pip2d(as.matrix(gates %>% dplyr::filter(phase=="M") %>% dplyr::select(1,2)), as.matrix(data.frame(nTurq, n488))),
                   Sen_check = pip2d(as.matrix(gates %>% dplyr::filter(phase=="Sen") %>% dplyr::select(1,2)), as.matrix(data.frame(nTurq, n488))),
                   phase = case_when(G1_check==1~"G1",
                                     S_check==1~"S",
                                     G2_check==1~"G2",
                                     M_check==1~"M",
                                     Sen_check==1~"Sen"
                   ))
        
        
        phased_fp <- phased_fp %>% rbind(temp)
    }
}

# Snapshot of phase proliferation 

sampled_S <- phased_fucci %>%
    dplyr::filter(Timepoint == 1,
                  phase == "S") %>%
    group_by(Treatments) %>%
    # sample_n(664) %>%
    unique()

phased_fp %>%
    dplyr::filter(Timepoint %in% c(1),
                  Replicate != "Rep6") %>% 
    na.omit() %>% 
    arrange(phase) %>% 
    ggplot(aes(x=nTurq, y=n488))+
    geom_point(aes(fill=phase, colour=phase), alpha=0.8,size=4, shape=21)+
    
    scale_colour_manual(values = c("gray","gray","gray","black","gray"))+
    scale_fill_manual(values = c("gray","gray","gray","#29B473","gray"))+
    facet_wrap(Timepoint~Treatments, nrow=2)+
    theme_bw()

# S phase progression 
phase = c("G1","G2","S","M")
Timepoint = tp
Replicate = c("Rep3","Rep4","Rep6")
Condition = c("Control\nUnconfined",
              "Control\nConfined",
              "Oligomycin\nUnconfined",
              "Oligomycin\nConfined")

temp_phase = expand.grid(Replicate, Condition, Timepoint, phase) %>% 
    set_names("Replicate", "Treatments", "Timepoint", "phase")

S_phased_fp <- phased_fp %>%
    dplyr::filter(!phase %in% c("Sen",NA),
                  cellID %in% startPhase$cellID[startPhase$phase=="S"]) %>% 
    arrange(Timepoint)

insertion_phase_1 <- S_phased_fp %>% 
    dplyr::select(Treatments, cellID) %>%
    unique() %>% 
    mutate(Treatments = as.character(Treatments))

insertion_phase_2 <-  expand.grid(cellID = insertion_phase_1$cellID, 
                                  Timepoint = S_phased_fp %>%
                                      dplyr::filter(Timepoint>=1) %>% 
                                      .$Timepoint %>% 
                                      unique()) %>% 
    as_tibble() %>% 
    mutate(cellID = as.character(cellID)) 



insertion_phase <-left_join(insertion_phase_1, insertion_phase_2) %>% 
    mutate(phase="Tracking\nend")


div_df_s <- bind_rows(S_phased_fp, insertion_phase) %>% 
    dplyr::select(Treatments, Timepoint, cellID, phase) %>% 
    dplyr::filter(Timepoint>=1) %>% 
    group_by(cellID, Timepoint) %>% 
    arrange(Timepoint, phase) %>% 
    slice(1) %>% 
    ungroup() %>% 
    group_by(Treatments, Timepoint, phase) %>% 
    tally() %>% 
    ungroup()


unique_Timepoint = S_phased_fp$Timepoint %>% unique()
unique_phase = c("S", "G2", "M", "Tracking\nend")

temp_phase = expand.grid(unique_Timepoint, unique_phase) %>% 
    set_names("Timepoint", "phase")

abundance_fp_mod <- left_join(temp_phase, div_df_s)  %>% 
    mutate(n = replace_na(n, 0)) %>% 
    group_by(Treatments, Timepoint) %>% 
    mutate(total_n = sum(n),
           percentage_n = signif(100*(n/total_n), 6)) %>% 
    ungroup() %>% 
    arrange(Timepoint) %>% 
    mutate(across(Treatments, factor, levels=c("Control\nUnconfined",
                                               "Control\nConfined",
                                               "Oligomycin\nUnconfined",
                                               "Oligomycin\nConfined")),
           across(phase, factor, levels = c("M", "G2","S"))) %>% 
    na.omit()

test_data = abundance_fp_mod %>% 
    dplyr::filter(phase=="S") %>% 
    group_by(Treatments, Timepoint) %>% 
    dplyr::select(-phase) %>% 
    unique() 

stat.ctrl.s <- wilcox.test(percentage_n ~ Treatments, data=test_data %>% 
                               dplyr::filter(Treatments %in% c("Control\nUnconfined", "Control\nConfined")))
# 
stat.oligo.s <- wilcox.test(percentage_n ~ Treatments, data=test_data %>% 
                                dplyr::filter(Treatments %in% c("Oligomycin\nUnconfined",
                                                                "Oligomycin\nConfined")))

stat.ctrl_oligo.s <- wilcox.test(percentage_n ~ Treatments, data=test_data %>% 
                                     dplyr::filter(Treatments %in% c("Control\nUnconfined",
                                                                     "Oligomycin\nUnconfined")))

ggplot(abundance_fp_mod %>% dplyr::filter(phase=="S"), aes(x=Timepoint, y=log2(percentage_n))) + 
    geom_line(aes(colour=Treatments)) +
    theme_bw() +
    ylab("Log2 percentage population")

# S phase ridge plot progression over time 
S_phased_fp %>% 
    dplyr::filter(Timepoint %in% c(1,3,6,9,12,15,18,21,24,27,30,33,36),
                  phase == "S") %>%
    na.omit() %>% 
    ggplot(aes(x=n488, y=factor(-Timepoint),fill = stat(x),height = ..count..))+
    geom_density_ridges_gradient(stat = "density")+
    scale_fill_viridis_c(name = "n488", option = "D") +
    facet_wrap(~Treatments, nrow=1)+
    theme_bw() +
    xlim(0,1)


# G2 progression over time

G2_phased_fp <- phased_fp %>%
    dplyr::filter(!phase %in% c("Sen",NA),
                  cellID %in% startPhase$cellID[startPhase$phase=="G2"]) %>% 
    arrange(Timepoint)

insertion_phase_1 <- G2_phased_fp %>% 
    dplyr::select(Treatments, cellID) %>%
    unique() %>% 
    mutate(Treatments = as.character(Treatments)) 

insertion_phase_2 <-  expand.grid(cellID = insertion_phase_1$cellID, 
                                  Timepoint = G2_phased_fp %>%
                                      dplyr::filter(Timepoint>=1) %>% 
                                      .$Timepoint %>% 
                                      unique()) %>% 
    as_tibble() %>% 
    mutate(cellID = as.character(cellID)) 



insertion_phase <-left_join(insertion_phase_1, insertion_phase_2) %>% 
    mutate(phase="Tracking\nend")


div_df_g2 <- bind_rows(G2_phased_fp, insertion_phase) %>% 
    dplyr::select(Treatments, Timepoint, cellID, phase) %>% 
    group_by(cellID, Timepoint) %>% 
    arrange(Timepoint, phase) %>% 
    slice(1) %>% 
    ungroup() %>% 
    group_by(Treatments, Timepoint, phase) %>% 
    tally() %>% 
    ungroup()


unique_Timepoint = G2_phased_fp$Timepoint %>% unique()
unique_phase = c("G2", "M", "Tracking\nend")

temp_phase = expand.grid(unique_Timepoint, unique_phase) %>% 
    set_names("Timepoint", "phase")

abundance_fp_mod <- left_join(temp_phase, div_df_g2)  %>% 
    mutate(n = replace_na(n, 0)) %>%
    dplyr::filter(Timepoint>=1) %>% 
    group_by(Treatments, Timepoint) %>% 
    mutate(total_n = sum(n),
           percentage_n = signif(100*(n/total_n), 6)) %>% 
    ungroup() %>% 
    arrange(Timepoint) %>% 
    mutate(across(Treatments, factor, levels=c("Control\nUnconfined",
                                               "Control\nConfined",
                                               "Oligomycin\nUnconfined",
                                               "Oligomycin\nConfined")),
           across(phase, factor, levels = c("M", "G2"))) %>% 
    na.omit()

test_data = abundance_fp_mod %>% 
    dplyr::filter(phase=="G2") %>% 
    group_by(Treatments, Timepoint) %>% 
    dplyr::select(-phase) %>% 
    unique() 

stat.ctrl.g2 <- wilcox.test(percentage_n ~ Treatments, data=test_data %>% 
                                dplyr::filter(Treatments %in% c("Control\nUnconfined", "Control\nConfined")))
# 
stat.oligo.g2 <- wilcox.test(percentage_n ~ Treatments, data=test_data %>% 
                                 dplyr::filter(Treatments %in% c("Oligomycin\nUnconfined",
                                                                 "Oligomycin\nConfined")))

stat.ctrl_oligo.g2 <- wilcox.test(percentage_n ~ Treatments, data=test_data %>% 
                                      dplyr::filter(Treatments %in% c("Control\nUnconfined",
                                                                      "Oligomycin\nUnconfined")))

ggplot(abundance_fp_mod %>% dplyr::filter(phase=="G2"), aes(x=Timepoint, y=log2(percentage_n))) + 
    geom_line(aes(colour=Treatments)) +
    theme_bw() +
    ylab("Log2 percentage population")



# G1 progression
G1_phased_fp <- phased_fp %>%
    dplyr::filter(!phase %in% c("Sen",NA),
                  cellID %in% startPhase$cellID[startPhase$phase=="G1"]) %>% 
    arrange(Timepoint)

insertion_phase_1 <- G1_phased_fp %>% 
    dplyr::select(Treatments, cellID) %>%
    unique() %>% 
    mutate(Treatments = as.character(Treatments)) 

insertion_phase_2 <-  expand.grid(cellID = insertion_phase_1$cellID, 
                                  Timepoint = G1_phased_fp %>% 
                                      .$Timepoint %>% 
                                      unique()) %>% 
    as_tibble() %>% 
    mutate(cellID = as.character(cellID)) 



insertion_phase <-left_join(insertion_phase_1, insertion_phase_2) %>% 
    mutate(phase="Tracking\nend")


div_df_g1 <- bind_rows(G1_phased_fp, insertion_phase) %>% 
    dplyr::select(Treatments, Timepoint, cellID, phase) %>% 
    group_by(cellID, Timepoint) %>% 
    arrange(Timepoint, phase) %>% 
    slice(1) %>% 
    ungroup() %>% 
    group_by(Treatments, Timepoint, phase) %>% 
    tally() %>% 
    ungroup()


unique_Timepoint = G1_phased_fp$Timepoint %>% unique()
unique_phase = c("G2", "M","S","G1","Tracking\nend")

temp_phase = expand.grid(unique_Timepoint, unique_phase) %>% 
    set_names("Timepoint", "phase")

abundance_fp_mod <- left_join(temp_phase, div_df_g1)  %>% 
    mutate(n = replace_na(n, 0))  %>%
    dplyr::filter(Timepoint >= 1) %>% 
    group_by(Treatments, Timepoint) %>% 
    mutate(total_n = sum(n),
           percentage_n = signif(100*(n/total_n), 6)) %>% 
    ungroup() %>% 
    arrange(Timepoint) %>% 
    mutate(across(Treatments, factor, levels=c("Control\nUnconfined",
                                               "Control\nConfined",
                                               "Oligomycin\nUnconfined",
                                               "Oligomycin\nConfined")),
           across(phase, factor, levels = c("M","G2","S","G1"))) %>% 
    na.omit()

test_data = abundance_fp_mod %>% 
    dplyr::filter(phase=="G1") %>% 
    group_by(Treatments, Timepoint) %>% 
    dplyr::select(-phase) %>% 
    unique() 

stat.ctrl.g1 <- wilcox.test(percentage_n ~ Treatments, data=test_data %>% 
                                dplyr::filter(Treatments %in% c("Control\nUnconfined", "Control\nConfined")))

stat.oligo.g1 <- wilcox.test(percentage_n ~ Treatments, data=test_data %>% 
                                 dplyr::filter(Treatments %in% c("Oligomycin\nUnconfined", "Oligomycin\nConfined")))

stat.ctrl_oligo.g1 <- wilcox.test(percentage_n ~ Treatments, data=test_data %>% 
                                      dplyr::filter(Treatments %in% c("Control\nUnconfined",
                                                                      "Oligomycin\nUnconfined")))

ggplot(abundance_fp_mod %>% dplyr::filter(phase=="G1"), aes(x=Timepoint, y=log2(percentage_n))) + 
    geom_line(aes(colour=Treatments)) +
    theme_bw() +
    ylab("Log2 percentage population")


# Time spent in phase - G1, S, G2
G1_time <- phased_fp %>% 
    dplyr::filter(!phase %in% c("Sen",NA),
                  cellID %in% startPhase$cellID[startPhase$phase=="G1"]) %>% 
    dplyr::filter(phase=="G1") %>% 
    group_by(Treatments, cellID) %>% 
    arrange(Timepoint) %>% 
    tally() %>% 
    mutate(n=n/3) # data was acquired every 20 min. Divide by three for hourly representation.


G1_time_stat <- compare_means(data = G1_time %>% as.data.frame(), 
                              n ~ Treatments,
                              method='wilcox.test') %>% 
    slice(1,2,5,6)

ggplot(G1_time, aes(x=Treatments, y=n))+
    geom_boxplot(aes(fill=Treatments), outlier.shape=NA)+
    geom_jitter(alpha=0.2,width=0.2, size=0.3)+
    theme_bw()+
    ylab("Duration detected in phase (h)")+
    ggtitle("Cells confined in G1 phase")+
    scale_fill_manual(values = c(jco2, jco2))+
    stat_pvalue_manual(G1_time_stat,
                       label = "p.signif",
                       y.position=c(40,41,44,40))


S_time <- phased_fp %>% 
    dplyr::filter(!phase %in% c("Sen",NA),
                  cellID %in% startPhase$cellID[startPhase$phase=="S"]) %>% 
    dplyr::filter(phase=="S") %>% 
    group_by(Treatments, cellID) %>% 
    arrange(Timepoint) %>% 
    tally() %>% 
    mutate(n = n/3)

S_time_stat = compare_means(data = S_time %>% as.data.frame(), 
                            n ~ Treatments,
                            method='wilcox.test') %>% 
    slice(1,2,5,6)

ggplot(S_time, aes(x=Treatments, y=n))+
    geom_boxplot(aes(fill=Treatments), outlier.shape=NA)+
    geom_jitter(alpha=0.2,width=0.2, size=0.3)+
    theme_bw()+
    ylab("Duration detected in phase (h)")+
    ggtitle("Cells confined in S phase")+
    scale_fill_manual(values = c(jco2, jco2))+
    stat_pvalue_manual(S_time_stat,
                       label = "p.signif",
                       y.position=c(40,41,44,40))


G2_time <- phased_fp %>% 
    dplyr::filter(!phase %in% c("Sen",NA),
                  cellID %in% startPhase$cellID[startPhase$phase=="G2"]) %>% 
    dplyr::filter(phase=="G2") %>% 
    group_by(Treatments, cellID) %>% 
    arrange(Timepoint) %>% 
    tally() %>% 
    mutate(n = n/3)

G2_time_stat = compare_means(data = G2_time %>% as.data.frame(), 
                             n ~ Treatments,
                             method='wilcox.test') %>% 
    slice(1,2,5,6)

ggplot(G2_time, aes(x=Treatments, y=n))+
    geom_boxplot(aes(fill=Treatments), outlier.shape=NA)+
    geom_jitter(alpha=0.2,width=0.2, size=0.3)+
    theme_bw()+
    ylab("Duration detected in phase (h)")+
    ggtitle("Cells confined in G2 phase")+
    scale_fill_manual(values = c(jco2, jco2))+
    stat_pvalue_manual(G2_time_stat,
                       label = "p.signif",
                       y.position=c(40,41,44,40))



