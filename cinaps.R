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
library(sjmisc)
library(circlize)
library(ComplexHeatmap)
library(data.table)
library(ggrepel)
library(biomaRt)
library(org.Hs.eg.db)
library(clusterProfiler)

jco = c("#EFC000FF", "#0073C2FF", "#CD534CFF", "#868686FF")
jco2 = c("#EFC000FF", "#0073C2FF")


# Nuc ATP Latrunculin  --------------------
atp_raw <- read_excel("supplementary_data_1.xlsx", sheet = "atp_fret") 

atp_pre <- atp_raw %>% 
    filter(Exp %in% c(
        # "Exp001", #DMEM
        # "Exp002", #DMEM-Glu
        # "Exp003", #DMEM-Glu
        # "Exp004", #DMEM-Glu
        # "Exp005", #DMEM-Glu
        "Exp006", #DMEM
        "Exp015", #DMEM
        "Exp016", #DMEM
        "Exp017" #DMEM
        
    ),
    Media == "DMEM",
    Treatment %in% c("Control", "Latrunculin A (500 nM)")) 

atp_pre_nuc <- atp_pre %>% 
    filter(Staining == "nuc", Cell_line == "hela", Treatment!= "CGP", Condition!= "release") %>% 
    group_by(Date, Exp, Condition, Treatment, Cell_line, CellID) %>% 
    summarise(meanATP = mean(Mean)) %>%
    ungroup() %>% 
    group_by(Condition, Treatment, Cell_line) %>% 
    add_tally() %>% 
    mutate(labelwithN = paste0(Treatment,"\n", Condition, "\nn=", n)) %>% 
    ungroup() %>% 
    mutate(across(Condition, factor, levels = c("susp", "conf"))) 
# across(labelwithN, factor, levels = c("Control\nsusp\nn=61", 
#                                       "Control\nconf\nn=90",
#                                       "Latrunculin A (500 nM)\nsusp\nn=26",
#                                       "Latrunculin A (500 nM)\nconf\nn=36")))



stat.lata <- compare_means(data = atp_pre_nuc, meanATP ~ labelwithN, test = "wilcox")

ggplot(atp_pre_nuc, aes(x=labelwithN, y=meanATP))+
    geom_boxplot(aes(fill = Condition), outlier.shape = NA) +
    geom_jitter(width=0.1, aes(fill = Exp), shape = 21) +
    theme_bw() +
    # facet_wrap(~Exp)
    # scale_fill_manual(values = jco) +
    stat_pvalue_manual(stat.lata, label = "p.format", y.position = 5, step.increase = 0.08)

# U2OS Nuclear ATP 
#Stats
#Batch correction used in case of u2os because independent experiments show high variability unlike hela
stat.test.nuc.u2os <- compare_means(data=atpsensor_nuc_corrected %>% dplyr::filter(Cell_line=='u2os'), 
                                    meanATP_corrected~labelwithN,
                                    # group.by = "Treatment",
                                    method='wilcox.test') %>% 
    slice(1,5,6)


lm_meanATP <- lm(meanX ~ as.factor(Date), data = atpsensor_nuc)
noeffrep_ATP <- residuals(lm_meanATP) + coef(lm_meanATP)[1]

atpsensor_nuc_corrected <- cbind(atpsensor_nuc ,  
                                 meanATP_corrected = noeffrep_ATP) %>% 
    as_tibble()

#Plot U2OS Nuc ATP
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




# Nuc ATP Release --------------------
# These experiments are one prior to revision (001), and one during revision (017)
atp_revision_release_ <- atp_raw %>% 
    filter(Exp %in% c("Exp001", "Exp017")) 

atp_revision_release <- atp_revision_release_ %>% 
    filter(Staining == "nuc", 
           Cell_line == "hela") %>% 
    group_by(Date, Exp, Condition, Treatment, Cell_line, CellID) %>% 
    summarise(meanATP = mean(Mean)) %>%
    ungroup() %>% 
    group_by(Condition, Treatment, Cell_line) %>% 
    add_tally() %>% 
    mutate(labelwithN = paste0(Condition, "\nn=", n)) %>% 
    ungroup() %>% 
    mutate(across(Condition = factor, levels = c("susp", "conf", "release")),
           across(labelwithN, factor, levels = c("susp\nn=28", 
                                                 "conf\nn=36",
                                                 "release\nn=25"))) 

#Normalizing values due to very different experimental values; however indepedent experiments show the same result
glm_meanATP <- glm(as.numeric(meanATP) ~ Exp, data = atp_revision_release)
noeffrep_ATP <- residuals(glm_meanATP) + coef(glm_meanATP)[1]

#get time of imaging for each cell.
releaseTime <- read_excel("supplementary_data_1.xlsx", sheet = "nucATP_releaseTime") %>% 
    select(Date, Exp, Condition, Treatment, Cell_line, CellID, Time_rel) %>% 
    filter(Exp %in% c("Exp001", "Exp017")) %>% 
    mutate(time = (as.numeric(format(as.POSIXct(Time_rel), format = "%M")))+1)

atp_revision_release <- cbind(atp_revision_release, meanATP_corrected = noeffrep_ATP) %>% 
    as_tibble() %>% 
    mutate(across(Condition = factor, levels = c("susp", "conf", "release")),
           across(labelwithN, factor, levels = c("susp\nn=28", 
                                                 "conf\nn=36",
                                                 "release\nn=25"))) %>% 
    left_join(releaseTime)
# mutate(label = paste(labelwithN, CellID,Exp, sep="_")) %>% 
# filter(label != "release\nn=26_011")

stats.release <- compare_means(data = atp_revision_release, meanATP_corrected ~ labelwithN, method = "wilcox")

ggplot(atp_revision_release, aes(x=labelwithN, y=meanATP_corrected))+
    geom_boxplot(aes(fill = Condition), outlier.shape = NA) +
    geom_jitter(width=0.1, aes(fill = Condition), shape = 21, size=2) +
    theme_bw() +
    scale_fill_manual(values = jco) +
    stat_pvalue_manual(stats.release, label = "p.format", y.position = 6.1, step.increase = 0.08)

# release.lm <- lm(meanATP_corrected ~ time, atp_revision_release %>% filter(Condition == "release"))
# 
# eq <- paste0("y = ", round(release.lm$coefficients[2], 3), "x + ", 
#              round(release.lm$coefficients[1],3), "\nR2 = ", round(summary(release.lm)$r.squared, 5))


atp_revision_release %>% 
    filter(Condition == "release") %>% 
    ggplot(aes(x = time, y = meanATP_corrected))+
    geom_point(size = 5, shape =19 )+
    geom_smooth(method = lm) +
    # geom_abline(slope = coef(release.lm)[["time"]], 
    #             intercept = coef(release.lm)[["(Intercept)"]]) +
    # annotate("text", x = 1, y = release.lm$coefficients[1] + 0.1*release.lm$coefficients[1], 
    #          hjust = 0, label = eq) +
    theme_bw()+
    stat_regline_equation(label.y = 3,label.x = 5, aes(label = ..rr.label..))



# Nuc ATP Oligo  DMEM +/- Glucose ---------------------------

atp_revision_glu <- atp_raw %>%
    filter(Staining == "nuc", 
           Cell_line == "hela",
           Condition %in% c("susp", "conf"),
           Treatment %in% c("Control", "Oligomycin")) %>% # Control misbehaving
    group_by(Date, Exp, Condition, Treatment, Cell_line, CellID, Media) %>% 
    summarise(meanATP = mean(Mean)) %>%
    ungroup() %>% 
    mutate(Media = case_when(Media == "DMEM-Glu" ~ "-Glu",
                             Media == "DMEM" ~ "+Glu")) %>% 
    group_by(Condition, Treatment, Cell_line, Media) %>% 
    add_tally() %>% 
    mutate(labelwithN = paste(Media, Treatment, Condition, "n=", n, sep="\n")) %>% 
    ungroup() %>% 
    mutate(across(Condition, factor, levels = c("susp", "conf", "release")),
           across(Media, factor, levels = c("+Glu", "-Glu")),
           across(labelwithN, factor, levels = c("+Glu\nControl\nsusp\nn=\n142", 
                                                 "+Glu\nControl\nconf\nn=\n175",
                                                 "+Glu\nOligomycin\nsusp\nn=\n60",
                                                 "+Glu\nOligomycin\nconf\nn=\n83",
                                                 "-Glu\nControl\nsusp\nn=\n35",
                                                 "-Glu\nControl\nconf\nn=\n46",
                                                 "-Glu\nOligomycin\nsusp\nn=\n46",
                                                 "-Glu\nOligomycin\nconf\nn=\n53")))

stats.glucose <- compare_means(data = atp_revision_glu, meanATP ~ labelwithN, method = "wilcox")


ggplot(atp_revision_glu, aes(x = labelwithN, y = meanATP))+
    geom_boxplot(aes(fill = Condition), outlier.shape = NA) +
    geom_jitter(width=0.1, aes(fill = Condition), shape = 21, size = 2) +
    theme_bw() +
    scale_fill_manual(values = jco)+
    stat_pvalue_manual(stats.glucose, label = "p.format", y.position = 6.1, step.increase = 0.08, size = 2)

# +
#     geom_text(data = image_files_atp_glu, aes(label = paste0(Exp, "_", CellID)), check_overlap = F)


# # Image files for FP..
# image_files_atp_glu = atp_revision_glu %>%
#     group_by(labelwithN) %>%
#     # mutate(median_NAM = median(meanATP)) %>%
#     mutate(first=quantile(meanATP,probs=0.25),
#            second=quantile(meanATP,probs=0.5),
#            third=quantile(meanATP,probs=0.75)) %>%
#     ungroup() %>%
#     mutate(image_file_meanATP = case_when(labelwithN == "+Glu\nControl\nsusp\nn=\n116" &
#                                               meanATP < 1.1*second &
#                                               meanATP > 0.95*second ~ 1,
# 
#                                           labelwithN == "+Glu\nControl\nconf\nn=\n157" &
#                                               meanATP < 1.0*third &
#                                               meanATP > 0.95*third ~ 1,
# 
#                                           labelwithN == "+Glu\nOligomycin\nsusp\nn=\n30" &
#                                               meanATP < 1.1*second &
#                                               meanATP > 0.95*second ~ 1,
# 
#                                           labelwithN == "+Glu\nOligomycin\nconf\nn=\n51" &
#                                               meanATP < 1.05*second &
#                                               meanATP > 0.9*second ~ 1,
# 
#                                           labelwithN == "-Glu\nControl\nsusp\nn=\n35" &
#                                               meanATP < 1.3*second &
#                                               meanATP > 1.15*second ~ 1,
# 
#                                           labelwithN == "-Glu\nControl\nconf\nn=\n46" &
#                                               meanATP < 1.1*third &
#                                               meanATP > 0.95*third ~ 1,
# 
#                                           labelwithN == "-Glu\nOligomycin\nsusp\nn=\n46" &
#                                               meanATP < 1.1*second &
#                                               meanATP > 0.95*second ~ 1,
# 
#                                           labelwithN == "-Glu\nOligomycin\nconf\nn=\n53" &
#                                               meanATP < 1.05*second &
#                                               meanATP > 0.9*second ~ 1)) %>%
#     filter(image_file_meanATP == 1) %>%
#     arrange(Cell_line, Condition) %>%
#     mutate(Media = paste0("DMEM",Media),
#            labelwithN = paste0("DMEM",labelwithN))
# 
# write_csv(image_files_atp_glu, "image_files_atp_glu.csv")


# Nuc ATP Oligo Pyruvate ------------------------------------------------

atp_oligo_pyruvate_ <- atp_raw %>% 
    # filter(Exp %in% c("Exp030", "Exp031","Exp033","Exp034"),
    filter(Exp %in% c("Exp001", "Exp006", "Exp013", "Exp014", "Exp015", "Exp016", "Exp017", "Exp018", "Exp019", "Exp020", "Exp030", "Exp031",  "Exp033", "Exp034"),
           Treatment %in% c("Control", "Oligomycin", "Pyruvate", "Pyruvate_Oligomycin")) 


atp_oligo_pyruvate <- atp_oligo_pyruvate_ %>%
    filter(Staining == "nuc", Cell_line == "hela", Condition!= "release", Media == "DMEM") %>%
    group_by(Date, Exp, Condition, Treatment, Cell_line, CellID, Concentration_drug) %>%
    summarise(meanATP = mean(Mean)) %>%
    ungroup() %>%
    group_by(Condition, Treatment, Cell_line, Concentration_drug) %>%
    add_tally() %>%
    mutate(labelwithN = paste0(Concentration_drug, "\n", Treatment,"\n", Condition, "\nn=", n)) %>%
    ungroup() %>%
    mutate(across(Condition, factor, levels = c("susp", "conf")),
           across(labelwithN, factor, levels = c("0?M\nControl\nsusp\nn=142",
                                                 "0?M\nControl\nconf\nn=175",
                                                 "1?M\nOligomycin\nsusp\nn=60",
                                                 "1?M\nOligomycin\nconf\nn=83",
                                                 "3mM\nPyruvate\nsusp\nn=37",
                                                 "3mM\nPyruvate\nconf\nn=57",
                                                 "3mM_1?M\nPyruvate_Oligomycin\nsusp\nn=55",
                                                 "3mM_1?M\nPyruvate_Oligomycin\nconf\nn=59")))

stat.oligo.pyruvate <- compare_means(data = atp_oligo_pyruvate, meanATP ~ labelwithN, test = "wilcox")

ggplot(atp_oligo_pyruvate, aes(x=labelwithN, y=meanATP))+
    geom_boxplot(aes(fill = Condition), outlier.shape = NA) +
    geom_jitter(width=0.1, aes(fill = Condition), shape = 21) +
    theme_bw() +
    scale_fill_manual(values = jco) +
    stat_pvalue_manual(stat.oligo.pyruvate, label = "p.format", y.position = 5, step.increase = 0.04, size =3)

# +
#     geom_text(data = image_files_atp_bam_pyr, aes(label = paste0(Exp, "_", CellID)), check_overlap = F)


# # Image files for FP..
# image_files_atp_bam_pyr = atp_oligo_bam15 %>%
#     group_by(labelwithN) %>%
#     # mutate(median_NAM = median(meanATP)) %>%
#     mutate(first=quantile(meanATP,probs=0.25),
#            second=quantile(meanATP,probs=0.5),
#            third=quantile(meanATP,probs=0.75)) %>%
#     ungroup() %>%
#     mutate(image_file_meanATP = case_when(
#                                           
#                                           labelwithN == "3mM\nPyruvate\nsusp\nn=51" &
#                                               meanATP < 1.1*second &
#                                               meanATP > 0.9*second ~ 1,
#                                           
#                                           labelwithN == "3mM\nPyruvate\nconf\nn=52" &
#                                               meanATP < 1.05*second &
#                                               meanATP > 0.85*second ~ 1,
#                                           
#                                           labelwithN == "3mM_1?M\nPyruvate_Oligomycin\nsusp\nn=51" &
#                                               meanATP < 1.2*second &
#                                               meanATP > 1.0*second ~ 1,
#                                           
#                                           labelwithN == "3mM_1?M\nPyruvate_Oligomycin\nconf\nn=50" &
#                                               meanATP < 1.05*second &
#                                               meanATP > 0.85*second ~ 1)) %>%
#     filter(image_file_meanATP == 1) %>%
#     arrange(Cell_line, Condition) 
# 
# write_csv(image_files_atp_bam_pyr, "image_files_atp_bam_pyr.csv")


# Nuc ATP KO Cell ----------------------------------------------------

atp_ko_ <- atp_raw %>% 
    filter(Exp %in% c("Exp001", "Exp002", "Exp008", "Exp009", "Exp010", "Exp011", "Exp012",
                      "Exp013", "Exp014", "Exp015", "Exp016", 
                      "Exp023", "Exp024", "Exp025", "Exp027", "Exp028", "Exp029", "Exp031"),
           Treatment %in% c("Control")) 

atp_ko <- atp_ko_ %>% 
    filter(Staining == "nuc", Treatment == "Control", Condition != "release",
           Cell_line %in% c("hela", "hela_drp1KO", "hela_fis1KO", "hela_mfn1KO")) %>% 
    group_by(Date, Exp, Condition, Treatment, Cell_line, CellID) %>% 
    summarise(meanATP = mean(Mean)) %>%
    ungroup() %>% 
    group_by(Condition, Treatment, Cell_line) %>% 
    add_tally() %>% 
    mutate(labelwithN = paste0(Cell_line,"\n", Condition, "\nn=", n)) %>% 
    ungroup() %>% 
    mutate(across(Condition, factor, levels = c("susp", "conf")),
           across(labelwithN, factor, levels = c("hela\nsusp\nn=110", 
                                                 "hela\nconf\nn=148",
                                                 "hela_drp1KO\nsusp\nn=39",
                                                 "hela_drp1KO\nconf\nn=62",
                                                 "hela_fis1KO\nsusp\nn=37",
                                                 "hela_fis1KO\nconf\nn=49",
                                                 "hela_mfn1KO\nsusp\nn=39",
                                                 "hela_mfn1KO\nconf\nn=40")))

stat.ko.atp <- compare_means(data = atp_ko, meanATP ~ labelwithN, method = "wilcox")

ggplot(atp_ko, aes(x=labelwithN, y=meanATP))+
    geom_boxplot(aes(fill = Condition), outlier.shape = NA) +
    geom_jitter(width=0.1, aes(fill = Condition), shape = 21, size= 2) +
    theme_bw() +
    scale_fill_manual(values = jco) +
    stat_pvalue_manual(stat.ko.atp, label = "p.format", y.position = 5.1, step.increase = 0.08, size =2) 
# +
#     geom_text(data = image_files_atp_ko, aes(label = paste0(Exp, "_", CellID)), check_overlap = F)


# # Image files for FP..
# image_files_atp_ko = atp_ko %>%
#     group_by(labelwithN) %>%
#     # mutate(median_NAM = median(meanATP)) %>%
#     mutate(first=quantile(meanATP,probs=0.25),
#            second=quantile(meanATP,probs=0.5),
#            third=quantile(meanATP,probs=0.75)) %>%
#     ungroup() %>%
#     mutate(image_file_meanATP = case_when(labelwithN == "hela\nsusp\nn=111" &
#                                               meanATP < 1.05*second &
#                                               meanATP > 0.95*second ~ 1,
#                                           
#                                           labelwithN == "hela\nconf\nn=148" &
#                                               meanATP < 1.0*third &
#                                               meanATP > 0.9*third ~ 1,
#                                           
#                                           labelwithN == "hela_drp1KO\nsusp\nn=39" &
#                                               meanATP < 1.05*second &
#                                               meanATP > 0.95*second ~ 1,
#                                           
#                                           labelwithN == "hela_drp1KO\nconf\nn=62" &
#                                               meanATP < 1.00*third &
#                                               meanATP > 0.9*third ~ 1,
#                                           
#                                           labelwithN == "hela_fis1KO\nsusp\nn=37" &
#                                               meanATP < 1.15*first &
#                                               meanATP > 1.05*first ~ 1,
#                                           
#                                           labelwithN == "hela_fis1KO\nconf\nn=49" &
#                                               meanATP < 1.0*third &
#                                               meanATP > 0.9*third ~ 1,
#                                           
#                                           labelwithN == "hela_mfn1KO\nsusp\nn=39" &
#                                               meanATP < 1.0*second &
#                                               meanATP > 0.9*second ~ 1,
#                                           
#                                           labelwithN == "hela_mfn1KO\nconf\nn=40" &
#                                               meanATP < 1.05*second &
#                                               meanATP > 0.95*second ~ 1)) %>%
#     filter(image_file_meanATP == 1) %>%
#     arrange(Cell_line, Condition) 
# 
# write_csv(image_files_atp_ko, "image_files_atp_ko.csv")

#Normalised ATP
atp_ko_norm <- list()
for(j in atp_ko$Cell_line %>% unique){
    
    # Calculate mean of susp per cell line per experiment
    norm_mean = atp_ko %>% filter(
        Condition == "susp",
        Cell_line == j,
    ) %>%
        .$meanATP %>% mean(na.rm = TRUE)
    
    # Calculate sd of susp per cell line per experiment
    norm_sd = atp_ko %>% filter(
        Condition == "susp",
        Cell_line == j,
    )  %>%
        .$meanATP %>% sd(na.rm = TRUE)
    
    # Normalise
    # nam_ko_norm[[paste(j, k, l, sep="_")]] <- meanNAM_KO %>% 
    atp_ko_norm[[paste(j)]] <- atp_ko %>% 
        filter(
            Cell_line == j,
        ) %>%
        mutate(meanATP_normalised = (meanATP-norm_mean)/norm_sd)
}

# }
atp_ko_norm_df <- bind_rows(atp_ko_norm)     

stat.ko.atp.normal <- compare_means(data = atp_ko_norm_df, meanATP_normalised ~ labelwithN, method = "wilcox")

ggplot(atp_ko_norm_df, aes(x=labelwithN, y=meanATP_normalised))+
    geom_boxplot(aes(fill = Condition), outlier.shape = NA) +
    geom_jitter(width=0.1, aes(fill = Condition), shape = 21, size= 2) +
    theme_bw() +
    scale_fill_manual(values = jco)
stat_pvalue_manual(stat.ko.atp.normal, label = "p.format", y.position = 5.1, step.increase = 0.08, size =2) 


# NucATP Bam15 -----------------------------------------------

atp_bam15_ <- atp_raw %>% 
    filter(Treatment %in% c("Control", "Bam15")) 

atp_bam15 <- atp_bam15_ %>% 
    filter(Staining == "nuc", 
           Condition %in% c("susp", "conf"),
           Media == "DMEM",
           Cell_line %in% c("hela"),
           Concentration_drug %in% c("0?M", "10?M"),
           Exp %in% c("Exp001", "Exp006", "Exp030", "Exp031")) %>% 
    group_by(Date, Exp, Condition, Treatment, Cell_line, CellID) %>% 
    summarise(meanATP = mean(Mean)) %>%
    ungroup() %>% 
    group_by(Condition, Treatment, Cell_line) %>% 
    add_tally() %>% 
    mutate(labelwithN = paste0(Cell_line,"\n", Condition, "\nn=", n)) %>% 
    ungroup() %>% 
    mutate(across(Condition, factor, levels = c("susp", "conf")),
           across(labelwithN, factor, levels = c("hela\nsusp\nn=35", 
                                                 "hela\nconf\nn=40",
                                                 "hela\nsusp\nn=33", 
                                                 "hela\nconf\nn=34")))

stat.bam15 <- compare_means(data = atp_bam15, meanATP ~ labelwithN, test = "wilcox")

ggplot(atp_bam15, aes(x=labelwithN, y=meanATP))+
    geom_boxplot(aes(fill = Condition), outlier.shape = NA) +
    geom_jitter(width=0.1, aes(fill = Condition), shape = 21, size = 2) +
    theme_bw() +
    scale_fill_manual(values = jco) +
    stat_pvalue_manual(stat.bam15, label = "p.format", y.position = 5.1, step.increase = 0.08) 


# Mitochondrial ATP -----------------------------------------------
atpsensor_mito <- atp_raw %>%  
    mutate(Staining = tolower(Staining)) %>% 
    dplyr::filter(Staining %in%  c("nam", "mito_nam", "mito", "mito_all")) %>% 
    mutate(Category = case_when(Staining %in% c("nam", "mito_nam") ~ "NAM",
                                Staining %in% c("mito", "mito_all") ~ "Total Mito")) %>% 
    unique() %>% 
    group_by(Date, Exp, Condition, Treatment, Cell_line, CellID, Category) %>% 
    summarise(meanATP = mean(Mean)) %>%
    ungroup() %>% 
    group_by(Condition, Treatment, Category, Cell_line) %>% 
    add_tally() %>%
    ungroup() %>% 
    mutate(labelwithN = paste0(Category,"\n", Treatment, "\n",Condition, "\nn=", n), 
           across(Condition, factor, levels = c("susp", "conf")),
           across(labelwithN, factor, levels = c("Total Mito\nControl\nsusp\nn=11", 
                                                 "Total Mito\nControl\nconf\nn=47",
                                                 # "Total Mito\nOligomycin\nsusp\nn=8",
                                                 # "Total Mito\nOligomycin\nconf\nn=28",
                                                 
                                                 "NAM\nControl\nsusp\nn=11",
                                                 "NAM\nControl\nconf\nn=44"
                                                 # "NAM\nOligomycin\nsusp\nn=8",
                                                 # "NAM\nOligomycin\nconf\nn=28"
           )))

stat.test.mito <- compare_means(data = atpsensor_mito, 
                                meanATP ~ labelwithN,
                                method = 'wilcox.test') 

ggplot(atpsensor_mito , aes(x=labelwithN, y=meanATP))+
    geom_boxplot(aes(fill=Condition), colour='black', width=0.8, outlier.shape=NA, notch = F)+
    geom_jitter(shape=21,size=3, aes(fill=Condition), width = 0.1)+
    stat_pvalue_manual(stat.test.mito, label = "p.format", y.position = 5.1, step.increase = 0.08) +
    scale_fill_manual(values = jco)+
    scale_color_manual(values = jco)+
    theme_bw()+
    theme(panel.grid = element_blank(),axis.title.x = element_blank(), legend.position="none")+
    ylab("Mitochondrial ATP levels")


# Nuc Volume Adhesion Susp Conf --------------------

nuc_vol_raw <- read_xlsx("supplementary_data_1.xlsx", sheet = "mitoloc")

nuc_vol_oligo_raw <- read_excel("supplementary_data_1.xlsx", sheet = "atp_fret")  %>%
    filter(Staining == "nuc", 
           Cell_line == "hela",
           Condition %in% c("susp", "conf"),
           Treatment %in% c("Control", "Oligomycin"))

nuc_vol_latA <-  nuc_vol_raw %>%  
    filter(Cell_line == "HeLa", 
           Treatment %in% c("Latrunculin A (500 nM)"),
           Condition %in% c("Susp.", "Conf.", "susp", "conf")) %>%
    group_by(Date, Exp, Condition, Treatment, CellID, Media) %>% 
    summarise(vol = sum(zstack*Area), .groups = "drop") %>% 
    group_by(Condition, Treatment) %>% 
    add_tally() %>% 
    ungroup()

nuc_vol_oligo <- nuc_vol_oligo_raw %>% 
    filter(Cell_line == "hela",
           Condition %in% c("Susp.", "Conf.", "susp", "conf")) %>% 
    group_by(Date, Exp, Condition, Treatment, CellID, Media) %>% 
    summarise(vol = sum(zstack*Area), .groups = "drop") %>% 
    group_by(Condition, Treatment, Media) %>% 
    add_tally() %>% 
    ungroup()

nuc_vol_ctrl <- nuc_vol_raw %>% 
    filter(Cell_line == "HeLa",
           Treatment == "Control",
           Condition %in% c("Susp.", "Conf.", "susp", "conf")) %>% 
    group_by(Date, Exp, Condition, Treatment, CellID, Media) %>% 
    summarise(vol = sum(zstack*Area), .groups = "drop") %>% 
    group_by(Condition, Treatment) %>% 
    add_tally() %>% 
    ungroup() 


combined_vol_data <- rbind(nuc_vol_latA, nuc_vol_oligo) %>% 
    mutate(Condition = case_when(Condition == "conf" ~ "Conf.",
                                 Condition == "susp" ~ "Susp.",
                                 TRUE ~ Condition),
           labelwithN = paste0(Treatment, "\n",Media, "\n", Condition, "\nn=", n)) %>% 
    mutate(across(labelwithN, factor, levels = c("Control\nDMEM\nSusp.\nn=142",
                                                 "Control\nDMEM\nConf.\nn=175",
                                                 "Oligomycin\nDMEM\nSusp.\nn=60",
                                                 "Oligomycin\nDMEM\nConf.\nn=83",
                                                 "Latrunculin A (500 nM)\nDMEM\nSusp.\nn=26",
                                                 "Latrunculin A (500 nM)\nDMEM\nConf.\nn=36",
                                                 "Control\nDMEM-Glu\nSusp.\nn=35",
                                                 "Control\nDMEM-Glu\nConf.\nn=46",
                                                 "Oligomycin\nDMEM-Glu\nSusp.\nn=46",
                                                 "Oligomycin\nDMEM-Glu\nConf.\nn=53")),
           across(Condition, factor, levels = c("Susp.", "Conf."))) %>% 
    filter(vol < 2500)

summary_data <- combined_vol_data %>% 
    group_by(Condition, Treatment, labelwithN) %>% 
    summarise(mean = mean(vol),
              sd = sd(vol))

stats.nuc.vol <- compare_means(data = combined_vol_data, vol ~ labelwithN, method = "t.test") %>% 
    slice(1,
          18,
          31,
          40,
          45)

stats.nuc.vol_2 <- compare_means(data = combined_vol_data, vol ~ labelwithN, method = "t.test") %>% 
    slice(2,4,
          11,13)

stats.nuc.vol_3 <- compare_means(data = combined_vol_data, vol ~ labelwithN, method = "t.test") %>% 
    slice(41, 44)

ggplot(combined_vol_data, aes(x = labelwithN, y = vol)) +
    geom_boxplot(aes(fill=Condition), outlier.shape = NA) +
    geom_jitter(width=0.1, aes(fill = Condition), shape = 21, size = 2)+
    scale_fill_manual(values = jco) +
    theme_bw() +
    stat_pvalue_manual(stats.nuc.vol, label = "p.format", y.position = 2500) +
    stat_pvalue_manual(stats.nuc.vol_2, label = "p.format", y.position = 2800, step.increase = 0.09) +
    stat_pvalue_manual(stats.nuc.vol_3, label = "p.format", y.position = 2800, step.increase = 0.09) +
    ylab("Nuclear Volume (?m3) x10^3")

ggplot(summary_data, aes(x = labelwithN, y = mean)) +
    geom_col(aes(fill = Condition))+
    geom_jitter(data = combined_vol_data, width=0.15, aes(x = labelwithN, y = vol), shape = 21, size = 2, fill = "black", alpha = 0.2)+
    geom_errorbar(aes(ymin = mean-sd, ymax = mean+sd), width = 0.2)+
    scale_fill_manual(values = jco) +
    theme_bw() +
    stat_pvalue_manual(stats.nuc.vol, label = "p.format", y.position = 2500) +
    stat_pvalue_manual(stats.nuc.vol_2, label = "p.format", y.position = 2800, step.increase = 0.09) +
    stat_pvalue_manual(stats.nuc.vol_3, label = "p.format", y.position = 2800, step.increase = 0.09) +
    ylab("Nuclear Volume (?m3) x10^3")




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


# NAM WT ------------------------------
mitoloc <- read_xlsx("supplementary_data_1.xlsx", sheet = "mitoloc") 

# Compute means
meanNAM <- mitoloc %>%
    dplyr::filter(Treatment=="Control",
                  Cell_line %in% c("HeLa","MiaPaca2","MDA-MB-231","U2OS"),
                  Condition %in% c("Susp.", "Conf."),
                  Exp %in% c("Exp001","Exp002","Exp003","Exp004",
                             "Exp005","Exp006","Exp007","Exp008",
                             "Exp009","Exp010","Exp011","Exp012",
                             "Exp013","Exp014","Exp015","Exp016"),
                  Concentration_stain %in% c("500 nM", "200 nM", "100 nM", "50 nM")) %>% 
    group_by(Treatment, Condition, Exp, Cell_line, CellID, Concentration_stain) %>% 
    summarise(meanRID = mean(RawIntDen),
              minSolidity = min(Solidity)) %>% 
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
                            level=c("Susp.\nn=81", "Conf.\nn=148",
                                    "Susp.\nn=29", "Conf.\nn=23",
                                    "Susp.\nn=19", "Conf.\nn=24",
                                    "Susp.\nn=37", "Conf.\nn=53")),
                   y=meanRID))+
    geom_boxplot(aes(fill=Condition),lwd=0.8, fatten=0.8, outlier.shape = NA)+
    geom_jitter(shape=21,size=3, aes(fill=Condition), width = 0.1, alpha=1)+
    scale_fill_manual(values=jco)+
    theme_bw()+
    theme(
        legend.position = "none",
        axis.title.x=element_blank(),
        axis.text.x=element_text(size=9))+
    ylab("NAM\n(Total mitotracker intensity)")+
    facet_wrap(~Cell_line, scales="free", nrow=1)+
    stat_pvalue_manual(stat.meanNAM,
                       label = "p.format",
                       y.position= "ypos")



# NAM Drugs ------------------------------
mitoloc <- read_xlsx("supplementary_data_1.xlsx", sheet = "mitoloc") 

# Only experiments with BAPTA, LAT (500) CK666 SMIFh2 (100) Jasplakinolide (500nM 30 min)
meanNAM_drugs <- mitoloc %>% 
    dplyr::filter(Cell_line == "HeLa", 
                  Exp %in% c("Exp008",
                             "Exp009",
                             "Exp018",
                             "Exp019",
                             "Exp020",
                             "Exp021",
                             "Exp022",
                             "Exp023",
                             "Exp033",
                             "Exp035",
                             "Exp036",
                             "Exp049"),
                  Treatment %in% c("Control", 
                                   "BAPTA-AM", 
                                   "Latrunculin A (500 nM)", 
                                   "CK666", 
                                   "SMIFH2",
                                   "Jasplakinolide",
                                   "Bam15"
                  ),
                  Concentration_stain %in% c("500 nM", "200 nM", "100 nM", "50 nM"))

meanNAM_drugs_separate <- meanNAM_drugs %>% 
    group_by(Treatment, Condition, Exp, CellID, Concentration_drug) %>% 
    summarise(meanRID = mean(RawIntDen)) %>% 
    ungroup() %>% 
    
    mutate(Concentration_drug = case_when(Concentration_drug == "500nM_20m" ~ "500nM_30m",
                                          Concentration_drug != "500nM_20m" ~ Concentration_drug)) %>% 
    
    group_by(Treatment, Condition, Concentration_drug) %>% 
    add_tally() %>% 
    ungroup() %>% 
    mutate(labelwithN = paste0(Treatment, "\n", Concentration_drug, "\n", Condition,"\nn=",n),
           across(Condition, factor, levels=c("Susp.", "Conf.")))


# Perform statistics for both the categories
stat.meanNAM.separate <- compare_means(data = meanNAM_drugs_separate,
                                       meanRID ~ labelwithN,
                                       method = "wilcox.test") 



# Plot the NAM for all various drug treatments, but independently, with their own respective controls

ggplot(meanNAM_drugs_separate %>% mutate(Filt = paste0(Treatment,meanRID)) %>% 
           # Not plotting extreme outlier for nicer plot. Note, they were used for stats.
           dplyr::filter(Filt != "BAPTA-AM2311570.69230769",
                         Filt != "Latrunculin A (500 nM)852615.75"),
       aes(
           x=factor(
               labelwithN,level=c(
                   "Control\n0?M\nSusp.\nn=74",
                   "Control\n0?M\nConf.\nn=144",
                   "BAPTA-AM\n10?M\nSusp.\nn=25",
                   "BAPTA-AM\n10?M\nConf.\nn=47",
                   "Latrunculin A (500 nM)\n500nM\nSusp.\nn=26",
                   "Latrunculin A (500 nM)\n500nM\nConf.\nn=36",
                   "CK666\n100?M\nSusp.\nn=32",
                   "CK666\n100?M\nConf.\nn=64",
                   "SMIFH2\n100?M\nSusp.\nn=34",
                   "SMIFH2\n100?M\nConf.\nn=65",
                   "Jasplakinolide\n500nM_30m\nSusp.\nn=25",
                   "Jasplakinolide\n500nM_30m\nConf.\nn=34",
                   "Bam15\n10?M\nSusp.\nn=18",
                   "Bam15\n10?M\nConf.\nn=24",
                   "Bam15\n50?M\nSusp.\nn=23",
                   "Bam15\n50?M\nConf.\nn=30"
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
                                             "SMIFH2",
                                             "Jasplakinolide",
                                             "Bam15")), nrow=2, scales="free")+
    ylab("NAM\n(Total mitotracker intensity)")
# stat_pvalue_manual(stat.meanNAM.separate,
#                    label = "p.signif")



#Normalizing values to present data with single control suspension
#Zero mean normalisation
meanNAM_drugs_control_susp = dplyr::filter(meanNAM_drugs_separate, Treatment=="Control", Condition=="Susp.") %>% 
    mutate(meanRID_mean = mean(meanRID),sdRID = sd(meanRID),
           newRID = ((meanRID - meanRID_mean)/sdRID)
    )
meanNAM_drugs_control_conf = dplyr::filter(meanNAM_drugs_separate, Treatment=="Control", Condition=="Conf.") %>% 
    mutate(
        newRID = 
            ((meanRID-(meanNAM_drugs_control_susp$meanRID_mean %>% unique()))/(meanNAM_drugs_control_susp$sdRID %>% unique()))
    )

meanNAM_drugs_bapta_susp = dplyr::filter(meanNAM_drugs_separate, Treatment=="BAPTA-AM", Condition=="Susp.") %>% 
    mutate(meanRID_mean = mean(meanRID),sdRID = sd(meanRID),
           newRID = ((meanRID - meanRID_mean)/sdRID)
    )
meanNAM_drugs_bapta_conf = dplyr::filter(meanNAM_drugs_separate, Treatment=="BAPTA-AM", Condition=="Conf.") %>% 
    mutate(
        newRID = 
            ((meanRID-(meanNAM_drugs_bapta_susp$meanRID_mean %>% unique()))/(meanNAM_drugs_bapta_susp$sdRID %>% unique()))
    )

meanNAM_drugs_lat_susp = dplyr::filter(meanNAM_drugs_separate, Treatment=="Latrunculin A (500 nM)", Condition=="Susp.") %>% 
    mutate(meanRID_mean = mean(meanRID),sdRID = sd(meanRID),
           newRID = ((meanRID - meanRID_mean)/sdRID)
    )
meanNAM_drugs_lat_conf = dplyr::filter(meanNAM_drugs_separate, Treatment=="Latrunculin A (500 nM)", Condition=="Conf.") %>% 
    mutate(
        newRID = 
            ((meanRID-(meanNAM_drugs_lat_susp$meanRID_mean %>% unique()))/(meanNAM_drugs_lat_susp$sdRID %>% unique()))
    )

meanNAM_drugs_ck_susp = dplyr::filter(meanNAM_drugs_separate, Treatment=="CK666", Condition=="Susp.") %>% 
    mutate(meanRID_mean = mean(meanRID),sdRID = sd(meanRID),
           newRID = ((meanRID - meanRID_mean)/sdRID)
    )
meanNAM_drugs_ck_conf = dplyr::filter(meanNAM_drugs_separate, Treatment=="CK666", Condition=="Conf.") %>% 
    mutate(
        newRID = 
            ((meanRID-(meanNAM_drugs_ck_susp$meanRID_mean %>% unique()))/(meanNAM_drugs_ck_susp$sdRID %>% unique()))
    )

meanNAM_drugs_smifh2_susp_100 = dplyr::filter(meanNAM_drugs_separate, Treatment=="SMIFH2", Condition=="Susp.", labelwithN == "SMIFH2\n100?M\nSusp.\nn=34") %>% 
    mutate(meanRID_mean = mean(meanRID), sdRID = sd(meanRID),
           newRID = ((meanRID - meanRID_mean)/sdRID)
    )
meanNAM_drugs_smifh2_conf_100 = dplyr::filter(meanNAM_drugs_separate, Treatment=="SMIFH2", Condition=="Conf.", labelwithN == "SMIFH2\n100?M\nConf.\nn=65") %>%
    mutate(
        newRID = 
            ((meanRID-(meanNAM_drugs_smifh2_susp_100$meanRID_mean %>% unique()))/(meanNAM_drugs_smifh2_susp_100$sdRID %>% unique()))
    )


meanNAM_drugs_Jasplak_susp_30min = dplyr::filter(meanNAM_drugs_separate, Treatment=="Jasplakinolide", Condition=="Susp.", labelwithN == "Jasplakinolide\n500nM_30m\nSusp.\nn=25") %>% 
    mutate(meanRID_mean = mean(meanRID), sdRID = sd(meanRID),
           newRID = ((meanRID - meanRID_mean)/sdRID)
    )
meanNAM_drugs_Jasplak_conf_30min = dplyr::filter(meanNAM_drugs_separate, Treatment=="Jasplakinolide", Condition=="Conf.", labelwithN == "Jasplakinolide\n500nM_30m\nConf.\nn=34") %>% 
    mutate(
        newRID = 
            ((meanRID-(meanNAM_drugs_Jasplak_susp_30min$meanRID_mean %>% unique()))/(meanNAM_drugs_Jasplak_susp_30min$sdRID %>% unique()))
    )

meanNAM_drugs_bam_susp_10 = dplyr::filter(meanNAM_drugs_separate, Treatment=="Bam15", Condition=="Susp.", labelwithN == "Bam15\n10?M\nSusp.\nn=18") %>% 
    mutate(meanRID_mean = mean(meanRID), sdRID = sd(meanRID),
           newRID = ((meanRID - meanRID_mean)/sdRID)
    )
meanNAM_drugs_bam_conf_10 = dplyr::filter(meanNAM_drugs_separate, Treatment=="Bam15", Condition=="Conf.", labelwithN == "Bam15\n10?M\nConf.\nn=24") %>% 
    mutate(
        newRID = 
            ((meanRID-(meanNAM_drugs_bam_susp_10$meanRID_mean %>% unique()))/(meanNAM_drugs_bam_susp_10$sdRID %>% unique()))
    )


meanNAM_drugs_bam_susp_50 = dplyr::filter(meanNAM_drugs_separate, Treatment=="Bam15", Condition=="Susp.", labelwithN == "Bam15\n50?M\nSusp.\nn=23") %>% 
    mutate(meanRID_mean = mean(meanRID), sdRID = sd(meanRID),
           newRID = ((meanRID - meanRID_mean)/sdRID)
    )
meanNAM_drugs_bam_conf_50 = dplyr::filter(meanNAM_drugs_separate, Treatment=="Bam15", Condition=="Conf.", labelwithN == "Bam15\n50?M\nConf.\nn=30") %>% 
    mutate(
        newRID = 
            ((meanRID-(meanNAM_drugs_bam_susp_50$meanRID_mean %>% unique()))/(meanNAM_drugs_bam_susp_50$sdRID %>% unique()))
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
                                meanNAM_drugs_smifh2_conf_100,
                                meanNAM_drugs_Jasplak_susp_30min %>% dplyr::select(-meanRID_mean, -sdRID),
                                meanNAM_drugs_Jasplak_conf_30min,
                                meanNAM_drugs_bam_susp_10 %>% dplyr::select(-meanRID_mean, -sdRID),
                                meanNAM_drugs_bam_conf_10,
                                meanNAM_drugs_bam_susp_50 %>% dplyr::select(-meanRID_mean, -sdRID),
                                meanNAM_drugs_bam_conf_50) 

meanNAM_drugs_corrected = meanNAM_drugs_corrected %>% 
    mutate(labelwithN_comb = case_when(grepl("Susp.",labelwithN) == TRUE ~ "Susp.\nn=257", #Sum of all suspensions in this plot
                                       grepl("Conf.",labelwithN) == TRUE ~ labelwithN),
           across(Condition, factor, levels=c("Susp.", "Conf.")))



# Plot single Susp ctrl - corrected
stat.meanNAM.corrected <- compare_means(data = meanNAM_drugs_corrected,
                                        newRID ~ labelwithN_comb,
                                        method = "wilcox.test") %>% 
    slice(1:15)


# Plot all drugs standardised with controls, and hence plot single control.
ggplot(meanNAM_drugs_corrected, aes(x=factor(labelwithN_comb, 
                                             level=c("Susp.\nn=257",
                                                     "Control\n0?M\nConf.\nn=144",
                                                     "BAPTA-AM\n10?M\nConf.\nn=47",
                                                     "Latrunculin A (500 nM)\n500nM\nConf.\nn=36",
                                                     "CK666\n100?M\nConf.\nn=64",
                                                     "SMIFH2\n100?M\nConf.\nn=65",
                                                     "Jasplakinolide\n500nM_30m\nConf.\nn=34",
                                                     "Bam15\n10?M\nConf.\nn=24",
                                                     "Bam15\n50?M\nConf.\nn=30")),
                                    y=newRID))+
    geom_boxplot(lwd=0.8, fatten=0.8, outlier.shape = NA, aes(fill=Condition))+
    geom_jitter(shape=21,size=2, aes(fill=Condition), width = 0.1)+
    scale_fill_manual(values=jco)+
    theme_bw()+
    theme(
        panel.grid = element_blank(),
        legend.position = "none",
        axis.title.x=element_blank(),
        axis.text.x=element_text(size=12))+
    ylab("NAM\n(Total mitotracker intensity)")+
    stat_pvalue_manual(stat.meanNAM.corrected,  
                       label = "p.format", 
                       y.position=15,
                       step.increase = 0.05)
# geom_text(data = image_files_NAM_drugs, aes(label = paste0(Exp, "_", CellID)), check_overlap = F)


# Image files for FP. These are decided based on being within 10% of the median to be representative.
# image_files_NAM_drugs = meanNAM_drugs_corrected %>%
#     group_by(labelwithN_comb) %>%
#     mutate(median_NAM = median(meanRID),
#            median_NAM_normal = median(newRID)) %>%
#     ungroup() %>%
#     mutate(image_file_meanRID = case_when(meanRID < 1.1*median_NAM & meanRID > 0.9*median_NAM ~ 1,
#                                           meanRID >= 1.1*median_NAM | meanRID <= 0.9*median_NAM ~ 0),
#            image_file_newRID = case_when(newRID < 1.1*median_NAM_normal & newRID > 0.9*median_NAM_normal ~ 1,
#                                          newRID >= 1.1*median_NAM_normal | newRID <= 0.9*median_NAM_normal ~ 0)) %>%
#     filter(image_file_meanRID + image_file_newRID > 0)
# 
# write_csv(image_files_NAM_drugs, "image_files_NAM_drugs.csv", quote_escape = "none")


# Only BAM15
bam15_xclusive <- meanNAM_drugs_corrected %>% 
    filter(Exp %in% c("Exp019", "Exp022", # Controls
                      "Exp049"),
           Treatment %in% c("Control", "Bam15")) %>% # Bam15
    dplyr::select(-labelwithN, - labelwithN_comb, -n) %>%
    group_by(Treatment, Condition, Concentration_drug) %>% 
    add_tally() %>% 
    ungroup() %>% 
    mutate(labelwithN_comb = case_when(Condition == "Susp." ~ paste0("Control", "\n", Condition, "\nn = 69"),
                                       Condition == "Conf." ~ paste0(Treatment, "\n", Concentration_drug, "\n", Condition, "\nn = ", n)),
           labelwithN_comb = factor(labelwithN_comb, levels = c("Control\nSusp.\nn = 69",
                                                                "Control\n0?M\nConf.\nn = 38",
                                                                "Bam15\n10?M\nConf.\nn = 24",
                                                                "Bam15\n50?M\nConf.\nn = 30")))

# Plot single Susp ctrl - corrected
stat.bam15 <- compare_means(data = bam15_xclusive,
                            newRID ~ labelwithN_comb,
                            method = "wilcox.test") 


ggplot(bam15_xclusive, aes(x = labelwithN_comb, y = newRID))+
    geom_boxplot(lwd=0.8, fatten=0.8, outlier.shape = NA, aes(fill=Condition))+
    geom_jitter(shape=21,size=2, aes(fill=Condition), width = 0.1)+
    scale_fill_manual(values=jco)+
    theme_bw() +
    theme(panel.grid = element_blank()) +
    ylab("NAM\n(Total mitotracker intensity)")+
    stat_pvalue_manual(stat.bam15,
                       label = "p.format",
                       y.position=6,
                       step.increase = 0.09)

# NAM Oligo ------------------------
meanNAM_oligo <- mitoloc %>%
    dplyr::filter(Treatment %in% c("Control","Oligomycin"),
                  Cell_line %in% c("HeLa"),
                  Condition %in% c("Susp.", "Conf."),
                  Exp %in% c("Exp019","Exp022", "Exp037", "Exp051"),
                  Concentration_stain %in% c("500 nM", "200 nM", "100 nM", "50 nM")) %>% 
    group_by(Treatment, Condition, Exp, Cell_line, CellID, Concentration_stain) %>% 
    summarise(meanRID = mean(RawIntDen)) %>% 
    ungroup() %>% 
    group_by(Treatment, Condition, Cell_line) %>% 
    add_tally() %>% 
    ungroup() %>% 
    mutate(labelwithN=paste0(Treatment, "\n",Condition,"\nn=",n),
           across(Condition, factor, levels=c("Susp.", "Conf.")))


# Perform statistics for both the categories
stat.meanNAM.oligo <- compare_means(data = meanNAM_oligo,
                                    meanRID ~ labelwithN,
                                    method = "wilcox.test") 



# Plot the NAM for all various drug treatments, but independently, with their own respective controls
ggplot(meanNAM_oligo,
       aes(
           x=factor(
               labelwithN, level=c(
                   "Control\nSusp.\nn=51",
                   "Control\nConf.\nn=68",
                   "Oligomycin\nSusp.\nn=34", 
                   "Oligomycin\nConf.\nn=34" 
               )),
           y=meanRID
           
       ))+
    geom_boxplot(lwd=0.8, fatten=0.8, outlier.shape = NA, aes(fill=Condition))+
    scale_fill_manual(values=jco)+
    theme_bw()+
    theme(
        legend.position = "none",
        axis.title.x=element_blank(),
        axis.text.x=element_text(size=12))+
    geom_jitter(shape=21,size=3, aes(fill=Condition), width = 0.1)+
    ylab("NAM\n(Total mitotracker intensity)")+
    stat_pvalue_manual(stat.meanNAM.oligo, label = "p.format", y.position = 6e+05,
                       step.increase = 0.1)

# NAM KO ------------------------------
mitoloc <- read_xlsx("supplementary_data_1.xlsx", sheet = "mitoloc") 

meanNAM_KO <- mitoloc %>% 
    dplyr::filter(Cell_line %in% c("HeLa","HeLa_drp1KO", "HeLa_mfn1KO", "HeLa_fis1KO"), 
                  Exp %in% c("Exp016","Exp017","Exp018","Exp021",
                             "Exp038","Exp039","Exp040","Exp041",
                             "Exp042","Exp043","Exp044","Exp045",
                             "Exp046"),
                  Concentration_stain == "50 nM",
                  Treatment == "Control") %>% 
    group_by(Treatment, Condition, Exp, CellID, Cell_line) %>% 
    summarise(meanRID = mean(RawIntDen)) %>% 
    ungroup() %>% 
    group_by(Treatment, Condition, Cell_line) %>%
    
    # Filter data within 3SD becasue very extreme outliers are present in the data.
    mutate(median_meanRID = median(meanRID),
           sd_meanRID = sd(meanRID)) %>%
    filter(meanRID >= median_meanRID - (3 * sd_meanRID) &
               meanRID <= median_meanRID + (3 * sd_meanRID)) %>%
    ungroup() %>% 
    
    # Add n number
    group_by(Treatment, Condition, Cell_line) %>% 
    add_tally() %>% 
    ungroup() %>% 
    mutate(
        # Cell_line = case_when(Cell_line == "HeLa" ~ "HeLa",
        #                          Cell_line == "HeLa_drp1KO" ~ "HeLa",
        #                          Cell_line == "HeLa_mfn1KO" ~ "HeLa",
        #                          Cell_line == "HeLa_fis1KO" ~ "HeLa"),
        
        labelwithN = paste0(Cell_line, "\n", Condition,"\nn=",n),
        across(Condition, factor, levels=c("Susp.", "Conf.")))



# Perform statistics for both the categories
stat.meanNAM.KO <- compare_means(data = meanNAM_KO,
                                 meanRID ~ labelwithN,
                                 method = "wilcox.test") %>%
    slice(1,3,2, 23,25,24)  %>% 
    cbind(y.position = 1350000)

stat.meanNAM.KO_susp_conf <- compare_means(data = meanNAM_KO,
                                           meanRID ~ labelwithN,
                                           method = "wilcox.test") %>%
    slice(4, 11, 17, 22)  %>% 
    cbind(y.position = 1250000)


ggplot(meanNAM_KO, aes(
    x=factor(
        labelwithN,level=c(
            "HeLa\nSusp.\nn=42",
            "HeLa\nConf.\nn=67",
            "HeLa_drp1KO\nSusp.\nn=76",
            "HeLa_drp1KO\nConf.\nn=101",
            "HeLa_fis1KO\nSusp.\nn=67",
            "HeLa_fis1KO\nConf.\nn=92",
            "HeLa_mfn1KO\nSusp.\nn=58",
            "HeLa_mfn1KO\nConf.\nn=59"
        )),
    y=meanRID))+
    geom_boxplot(lwd=0.8, fatten=0.8, outlier.shape = NA, aes(fill=Condition))+
    scale_fill_manual(values=jco)+
    theme_bw()+
    theme(legend.position = "none",
          axis.title.x=element_blank(),
          axis.text.x=element_text(size=12))+
    geom_jitter(shape=21,size=3, aes(fill=Condition), width = 0.1)+
    ylab("NAM\n(Total mitotracker intensity)")+
    stat_pvalue_manual(stat.meanNAM.KO,
                       label = "p.format",
                       step.increase = 0.05,
                       tip.length = 0.01)+
    stat_pvalue_manual(stat.meanNAM.KO_susp_conf,
                       label = "p.format",
                       tip.length = 0.01)
# geom_text(data = image_files_NAM_KO, aes(label = paste0(Exp, "_", CellID)), check_overlap = F)


#For figures to be representative, we used the following cells:

# 
# image_files_NAM_KO = meanNAM_KO %>% 
#         group_by(labelwithN) %>%
#         # mutate(median_NAM = median(meanRID)) %>%
#         mutate(first=quantile(meanRID,probs=0.25),
#               second=quantile(meanRID,probs=0.5),
#               third=quantile(meanRID,probs=0.75)) %>% 
#         ungroup() %>%
#         mutate(image_file_meanRID = case_when(labelwithN == "HeLa\nSusp.\nn=42" &
#                                                   meanRID < 1.1*first & 
#                                                   meanRID > 0.9*first ~ 1,
#                                               
#                                               labelwithN == "HeLa\nConf.\nn=67" &
#                                                   meanRID < 1.3*third & 
#                                                   meanRID > 1.1*third ~ 1,
#                                               
#                                               labelwithN == "HeLa_drp1KO\nSusp.\nn=76" &
#                                                   meanRID < 1.1*second & 
#                                                   meanRID > 0.9*second ~ 1,
#                                               
#                                               labelwithN == "HeLa_drp1KO\nConf.\nn=101" &
#                                                   meanRID < 0.8*third & 
#                                                   meanRID > 0.6*third ~ 1,
#                                               
#                                               labelwithN == "HeLa_fis1KO\nSusp.\nn=67" &
#                                                   meanRID < 1.1*first & 
#                                                   meanRID > 0.9*first ~ 1,
#                                               
#                                               labelwithN == "HeLa_fis1KO\nConf.\nn=92" &
#                                                   meanRID < 1.1*third & 
#                                                   meanRID > 0.9*third ~ 1,
#                                               
#                                               labelwithN == "HeLa_mfn1KO\nSusp.\nn=58" &
#                                                   meanRID < 1.1*third & 
#                                                   meanRID > 0.9*third ~ 1,
#                                               
#                                               labelwithN == "HeLa_mfn1KO\nConf.\nn=59" &
#                                                   meanRID < 1.1*second & 
#                                                   meanRID > 0.9*second ~ 1)) %>% 
#     filter(image_file_meanRID == 1) %>% 
#     arrange(Cell_line, Condition)


# write_csv(image_files_NAM_KO, "image_files_NAM_KO.csv", quote_escape = "none")

# NAM KO z-score (classical) -----------------------------
nam_ko_norm <- list()
for(j in meanNAM_KO$Cell_line %>% unique){
    # for(k in meanNAM_KO$Exp %>% unique){
    for(l in meanNAM_KO$Treatment %>% unique){
        
        # Calculate mean of susp per cell line per experiment
        norm_mean = meanNAM_KO %>% filter(
            Condition == "Susp.",
            Treatment == l,
            Cell_line == j,
            # Exp == k
        ) %>%
            .$meanRID %>% mean(na.rm = TRUE)
        
        # Calculate sd of susp per cell line per experiment
        norm_sd = meanNAM_KO %>% filter(
            Condition == "Susp.",
            Treatment == l,
            Cell_line == j,
            # Exp == k
        )  %>%
            .$meanRID %>% sd(na.rm = TRUE)
        
        # Normalise
        # nam_ko_norm[[paste(j, k, l, sep="_")]] <- meanNAM_KO %>% 
        nam_ko_norm[[paste(j, l, sep="_")]] <- meanNAM_KO %>% 
            filter(
                Treatment == l,
                Cell_line == j,
                # Exp == k
            ) %>%
            mutate(meanRID_normalised = (meanRID-norm_mean)/norm_sd)
    }
}
# }
nam_ko_norm_df <- bind_rows(nam_ko_norm) 

# # Adding the minimum value as a constant to help visualisation 
# mutate(meanRID_normalised = meanRID_normalised + abs(min(meanRID_normalised, na.rm=T)))

# Perform statistics for both the categories
stat.meanNAM_normal.KO <- compare_means(data = nam_ko_norm_df,
                                        meanRID_normalised ~ labelwithN,
                                        method = "wilcox.test") %>%
    slice(9, 13, 11, 2, 6, 4)  %>% 
    cbind(y.position = 14.5)

stat.meanNAM_normal.KO_susp_conf <- compare_means(data = nam_ko_norm_df,
                                                  meanRID_normalised ~ labelwithN,
                                                  method = "wilcox.test") %>%
    slice(1, 14, 23, 28)  %>% 
    cbind(y.position = 13.5)


ggplot(nam_ko_norm_df, aes(
    x=factor(
        labelwithN,level=c(
            "HeLa\nSusp.\nn=42",
            "HeLa\nConf.\nn=67",
            "HeLa_drp1KO\nSusp.\nn=76",
            "HeLa_drp1KO\nConf.\nn=101",
            "HeLa_fis1KO\nSusp.\nn=67",
            "HeLa_fis1KO\nConf.\nn=92",
            "HeLa_mfn1KO\nSusp.\nn=58",
            "HeLa_mfn1KO\nConf.\nn=59"
        )),
    y=meanRID_normalised))+
    geom_boxplot(lwd=0.8, fatten=0.8, outlier.shape = NA, aes(fill=Condition))+
    scale_fill_manual(values=jco)+
    theme_bw()+
    theme(legend.position = "none",
          axis.title.x=element_blank(),
          axis.text.x=element_text(size=12))+
    geom_jitter(shape=21,size=3, aes(fill=Condition), width = 0.1)+
    ylab("NAM\n(Total mitotracker intensity)") +
    stat_pvalue_manual(stat.meanNAM_normal.KO,
                       label = "p.format",
                       step.increase = 0.05,
                       tip.length = 0.01) +
    stat_pvalue_manual(stat.meanNAM_normal.KO_susp_conf,
                       label = "p.format",
                       tip.length = 0.01)

# Indentation vs NAM ------------------------------
data_indent_nam <- read_excel("supplementary_data_1.xlsx", sheet = "indentation_vs_nam")

# Convert data to counts per category
data_summary_indent_nam <- data_indent_nam %>%
    group_by(Cell_line,Nuclear_Indentation, NAM) %>%
    summarise(Count = n()) %>%
    ungroup()

# # Create a stacked bar plot
# ggplot(data_summary, aes(x = interaction(Nuclear_Indentation, NAM), y = Count, fill = Nuclear_Indentation)) +
#     geom_bar(stat = "identity") +
#     labs(x = "Nuclear Indentation & NAM Combination", y = "Count", fill = "Nuclear Indentation") +
#     theme_minimal() +
#     scale_fill_brewer(palette = "Set2") +
#     ggtitle("Distribution of Nuclear Indentation and NAM Combinations")+
#     facet_wrap(~ Cell_line, scales="free")

ind_nam_percent <- data_summary_indent_nam %>% filter(NAM == "No")

ind_nam_percent <- ind_nam_percent %>%
    group_by(Cell_line) %>%
    mutate(Percentage = Count / sum(Count) * 100)

# Create stacked bar chart
ggplot(ind_nam_percent, aes(x = Cell_line, y = Percentage, fill = Nuclear_Indentation)) +
    geom_bar(stat = "identity") +
    labs(x = "Cell Line", y = "Percentage", title = "Nuclear Indentation\nPercentage per Cell Line") +
    scale_y_continuous(labels = scales::percent_format(scale = 1)) +
    theme_bw()+
    scale_fill_manual(values = c("gray", jco[2]))



# Time and NAM levels ------------------------------
time_NAM_df <- meanNAM_drugs_corrected %>% 
    dplyr::filter(Treatment == "Control", Exp!= "Exp019")

time_stamp <- read_xlsx("supplementary_data_1.xlsx", sheet = "time_nam")  %>% 
    dplyr::select(Treatment, Condition, labs, Exp, CellID, Time_rel)

time_NAM_df_stamp <- left_join(time_NAM_df, time_stamp) %>% 
    mutate(time = as.numeric(format(as.POSIXct(Time_rel), format = "%M"))) %>% 
    dplyr::filter(Condition == "Conf.")

ggplot(time_NAM_df_stamp %>% na.omit(), 
       aes(x=time, y=meanRID))+
    geom_point()+
    geom_smooth(method = lm)+
    theme_bw()+
    facet_wrap(~Condition, scales="free")+
    theme(panel.grid = element_blank())+
    stat_regline_equation(label.y = 12,label.x = 5, aes(label = ..rr.label..))+
    xlab("Time (h)")+
    ylab("NAM\n(Total mitotracker intensity)")


# Nuclear Shape ---------------------------------------------
nucShape_raw <- read_excel("supplementary_data_1.xlsx", sheet = "nucShape_data")

nucShape_n <- nucShape_raw %>% 
    filter(Shape_classification_confidence == 1) %>% 
    group_by(Condition, Shape) %>% 
    summarise(n=n()) %>%
    mutate(total=sum(n)) %>% 
    ungroup() %>% 
    mutate(percentage = n/total*100)

ggplot(nucShape_n, aes(x=Condition, y=n, fill=Shape))+
    geom_bar(position="fill", stat="identity")+
    scale_fill_manual(values = jco)+
    theme_bw()+
    theme(axis.title.x = element_blank())+
    ylab("Percentage")



# Shape and NAM -----------------------------------------------------------

nucShape_NAM_raw <- read_excel("supplementary_data_1.xlsx", sheet = "nucShape_NAM")

nucShape_NAM <- nucShape_NAM_raw %>% 
    filter(Satisfaction != 0) %>% 
    group_by(Condition, Shape) %>% 
    add_tally()

ggshape_mitotracker <- ggplot(nucShape_NAM,
                              aes(x=Shape, y=meanRID, fill=Shape,))+
    geom_boxplot(outlier.shape=21,
                 lwd=0.8, 
                 fatten=0.8, 
                 color="black")+
    geom_point(shape=21,size=3, aes(fill=Shape))+
    facet_wrap(~Condition)+
    theme_bw()+
    theme(axis.title.x=element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank())+
    # theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    color_palette(jco)+
    fill_palette(jco)+
    ylab("NAM\n(Mean raw-int-den)")


stat.nucshape.nam <- compare_means(data = nucShape_NAM, meanRID~Shape, group = "Condition")


# Perinuclear Actin -------------------------------------------------------

actin_data_raw <- read_xlsx("supplementary_data_1.xlsx", sheet = "actin") 

actin_data <- actin_data_raw %>% 
    group_by(Condition, Treatment, Filename) %>% 
    summarise(meanActin_inner = mean(Intensity_MeanIntensity_actin.1_inner),
              meanActin_outer = mean(Intensity_MeanIntensity_actin.1_outer),
              meanActin_inner_outer = mean(Intensity_MeanIntensity_actin.1_inner+Intensity_MeanIntensity_actin.1_outer),
              meanNuclearActin = mean(Intensity_MeanIntensity_actin_mask),
              .groups = "drop") %>% 
    group_by(Condition, Treatment) %>% 
    add_tally() %>% 
    ungroup() %>% 
    mutate(labelwithN = paste0(Condition,"\n", Treatment, "\nn=", n)) %>% 
    mutate(labelwithN = factor(labelwithN, levels = c("Adhesion\nCTRL\nn=20",
                                                      "Susp\nCTRL\nn=19",
                                                      "Conf\nCTRL\nn=28",
                                                      
                                                      "Adhesion\nLATA\nn=22",
                                                      "Susp\nLATA\nn=23",
                                                      "Conf\nLATA\nn=23",
                                                      
                                                      "Adhesion\nFIS1\nn=25",
                                                      "Susp\nFIS1\nn=16",
                                                      "Conf\nFIS1\nn=16",
                                                      
                                                      "Adhesion\nMFN1\nn=25",
                                                      "Susp\nMFN1\nn=16",
                                                      "Conf\nMFN1\nn=16")),
           
           Condition = factor(Condition, levels = c("Adhesion", "Susp", "Conf"))) %>% 
    pivot_longer(cols = 4:7, names_to = "analysis", values_to = "values") %>% 
    filter(analysis == "meanActin_inner_outer")

# stat.actin <- compare_means(data = actin_data, values ~ labelwithN, test = "wilcox")
# 
# ggplot(actin_data , aes(x = labelwithN, y = values))+
#     geom_boxplot(outlier.shape = NA, aes(fill = Condition)) +
#     geom_jitter(width = 0.2, shape = 21, aes(fill = Condition)) +
#     theme_bw() +
#     scale_fill_manual(values = c(jco[3], jco[1], jco[2])) +
#     stat_pvalue_manual(stat.actin, label = "p.format", y.position = 0.35, step.increase = 0.04)

#Since different cell lines were used, we will normalise the actin levels to the adhesion conditions.

# Zero mean normalisation to correct for differences in value range in different experiments. Normalisation is done within each cell line per experiment, and then combined.
periAct_norm <- list()

for(j in actin_data$Treatment %>% unique){
    
    # Calculate mean of susp per cell line per experiment
    norm_mean = actin_data %>% filter(
        Condition == "Adhesion",
        Treatment == j,
    ) %>%
        .$values %>% mean(na.rm = TRUE)
    
    # Calculate sd of susp per cell line per experiment
    norm_sd = actin_data %>% filter(
        Condition == "Adhesion",
        Treatment == j,
    ) %>%
        .$values %>% sd(na.rm = TRUE)
    
    # Normalise
    periAct_norm[[j]] <- actin_data %>%
        filter(
            Treatment == j
        ) %>%
        mutate(periAct_normalised = (values - norm_mean)/norm_sd)
}

periAct_norm_df <- bind_rows(periAct_norm) %>% 
    filter(periAct_normalised < 12) 

stat.actin.norm <- compare_means(data = periAct_norm_df, periAct_normalised ~ labelwithN, test = "wilcox")

ggplot(periAct_norm_df , aes(x = labelwithN, y = periAct_normalised))+
    geom_boxplot(outlier.shape = NA, aes(fill = Condition)) +
    geom_jitter(width = 0.2, shape = 21, aes(fill = Condition)) +
    theme_bw() +
    scale_fill_manual(values = c(jco[3], jco[1], jco[2])) +
    stat_pvalue_manual(stat.actin.norm, label = "p.format", y.position = 0.09, step.increase = 0.05, size = 2)

# Actin_Dapi_lineProfile --------------------------------------------------

totalsignal <- read_xlsx("supplementary_data_1.xlsx", sheet = "actin_dapi_latA") %>% 
    group_by(Treatment, Condition, Filename, Signal) %>% 
    mutate(rolling_avg = rollmean(yVal, k=3, fill=NA, align='right'),
           nVal = (yVal-min(yVal))/(max(yVal)-min(yVal)),
           nDistance = (Distance-min(Distance))/(max(Distance)-min(Distance))) %>% 
    ungroup()  %>% 
    mutate(across(Condition, factor, levels = c("Adhesion", "Susp", "Conf")),
           across(Treatment, factor, levels = c("Hela", "LatA", "Fis1KO", "Mfn1KO")))

# In the case where the nucleus is one side, by having nucleus on either side, the signal
# gets averaged and the trend is not visible. Insted, here, if the nucleus is on one side, we place it
# always on the right, and the signal becomes more evident. It is simply flipping the x coordinates.

flipCells_latSusp <- totalsignal %>% filter(Treatment == "LatA", Condition == "Susp") %>% 
    .$cellID %>% 
    as_tibble() %>% 
    unique() %>% 
    slice(5, 8, 11, 16) 

flipCells_latConf = totalsignal %>% filter(Treatment == "LatA", Condition == "Conf") %>% 
    .$cellID %>% 
    as_tibble() %>% 
    unique() %>% 
    slice(4,5,7,8,18,19,23)  

flipCells_fisSusp <- totalsignal %>% filter(Treatment == "Fis1KO", Condition == "Susp") %>%
    .$cellID %>%
    as_tibble() %>%
    unique() %>%
    slice(3, 12, 13, 14, 15)

flipCells_fisConf <- totalsignal %>% filter(Treatment == "Fis1KO", Condition == "Conf") %>%
    .$cellID %>%
    as_tibble() %>%
    unique() %>%
    slice(1,14)


flipCells_mfnAdh <- totalsignal %>% filter(Treatment == "Mfn1KO", Condition == "Adhesion") %>%
    .$cellID %>%
    as_tibble() %>%
    unique() %>%
    slice(3, 6, 11, 17)

flipCells_mfnSusp <- totalsignal %>% filter(Treatment == "Mfn1KO", Condition == "Susp") %>%
    .$cellID %>%
    as_tibble() %>%
    unique() %>%
    slice(10, 13) 

flipCells_mfnConf <- totalsignal %>% filter(Treatment == "Mfn1KO", Condition == "Conf") %>%
    .$cellID %>%
    as_tibble() %>%
    unique() %>%
    slice(2, 5, 16) 

flipCells = rbind(flipCells_latConf, flipCells_latSusp, 
                  flipCells_mfnAdh, flipCells_mfnSusp, flipCells_mfnConf
) %>% .$value

totalsignal_flip <- totalsignal %>% 
    mutate(flip_nDist = case_when(cellID %in% flipCells ~ 1-nDistance,
                                  !cellID %in% flipCells ~ nDistance))


# ggplot(totalsignal_flip  %>% filter(Treatment == "Fis1KO", Condition == "Conf"),
#        aes(x = nDistance, y = nVal, color = Signal, fill = Signal))+
#   # geom_line()+
#   # facet_grid(Treatment ~ Condition)+
#   facet_wrap(~ cellID)+
#   geom_smooth(stat = "smooth", position = "identity")+
#   scale_color_manual(values = c("Green", "Blue"))+
#   scale_fill_manual(values = c("Green", "Blue"))+
#   theme_bw()

ggplot(totalsignal_flip , 
       aes(x = flip_nDist, y = nVal, color = Signal, fill = Signal))+
    # geom_line()+
    facet_grid(Treatment ~ Condition)+
    # facet_wrap(~ cellID)+
    geom_smooth(stat = "smooth", position = "identity")+
    scale_color_manual(values = c("Green", "Blue"))+
    scale_fill_manual(values = c("Green", "Blue"))+
    theme_bw() 

# Kymograph Standard Deviation --------------------------------------------

kymo_data_raw <- read_xlsx("supplementary_data_1.xlsx", sheet = "er_kymograph_std_dev")

kymo_data <- kymo_data_raw %>% 
    group_by(condition, cell_id) %>% 
    mutate(plotting_group = paste0(condition, cell_id),
           StdIntensity = (Intensity-min(Intensity))/(max(Intensity)-min(Intensity)),
           StdPP = (Pixel_position-min(Pixel_position))/(max(Pixel_position)-min(Pixel_position)),
           cutPP = cut(StdPP, breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1))) %>%
    group_by(condition, cell_id, cutPP) %>%
    mutate(stdev = sd(StdIntensity),
           mInt = mean(StdIntensity)) %>% 
    ungroup() %>% 
    na.omit() %>% 
    dplyr::select(condition, cell_id, plotting_group, Time_point, cutPP, stdev, mInt, cell_line) %>% 
    unique() 

ggplot(kymo_data, aes(x=cutPP %>% as.numeric(), y=stdev))+
    # geom_point(aes(colour=condition))+
    geom_smooth(aes(linetype=condition))+
    # geom_line(aes(linetype = condition, colour = cell_id))+
    theme_bw() +
    ylab("Average standard deviation")+
    xlab("Standardized pixel position")


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
        compare_means(data =.,  coeffvar ~ Treat_Cond)  %>% 
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

# Stats
coeffvar_test_df <- cv_process %>% 
    group_by(Treat_Cond, Bin) %>% 
    summarise(mean_coeffvar = mean(coeffvar),.groups = "drop") %>% 
    pivot_wider(names_from = Treat_Cond, values_from = mean_coeffvar)

wilcox.test(coeffvar_test_df$`ControlConf.`,
            coeffvar_test_df$`ControlSusp.`)

wilcox.test(coeffvar_test_df$`OligomycinConf.`,
            coeffvar_test_df$`OligomycinSusp.`)

wilcox.test(coeffvar_test_df$`OligomycinSusp.`,
            coeffvar_test_df$`ControlSusp.`)

wilcox.test(coeffvar_test_df$`ControlConf.`,
            coeffvar_test_df$`ControlSusp.`)


# 53BP1 Eto/IR control experiments ----------------------
dna_damage_control <- read_xlsx("supplementary_data_1.xlsx", sheet = "eto_ir_dna_damage")

stat.eto.ir <- compare_means(data=dna_damage_control, 
                             Spot_Count~Treatment_Condition, 
                             group.by = "Treatment_group") %>% 
    slice(1:2,
          5:8,
          11:12) %>% 
    mutate(y.position = c(30,34,38,30,
                          44,48,52,44))

ggplot(dna_damage_control, aes(x = Treatment_Condition, y = Spot_Count))+
    geom_jitter(width = 0.2, shape = 19, alpha = 0.1)+
    geom_boxplot(outlier.shape = NA, aes(fill = Treatment), alpha = 0.9)+
    facet_wrap(~Treatment_group, scales = "free")+
    scale_fill_manual(values = jco)+
    theme_bw() +
    stat_pvalue_manual(stat.eto.ir, label.size = 2, label = "p.adj")+
    ylab("53BP1 foci count")

# 53BP1 under confinement ------------------------------------------------

# Oligomycin A
damage_data_raw <- read_xlsx("supplementary_data_1.xlsx", sheet = "53bp1_conf_oligo")

damage_data_ <- damage_data_raw %>% 
    dplyr::select(Repetition, Condition, Treatment, Area, CellID) %>% 
    dplyr::filter(Area < 3) %>%
    group_by(Repetition, Condition, Treatment, CellID) %>% 
    add_tally() %>% 
    group_by(Repetition, Condition, Treatment, CellID) %>% 
    mutate(mean_Area = mean(Area)) %>% 
    ungroup() %>% 
    dplyr::select(-Area,-CellID) %>% 
    unique() 

# Adding the zeros; 
# Susp Control = 13
# Conf Control = 6
# Susp Oligo = 9
# Conf Oligo = 16

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


damage_data <- damage_data_ %>% 
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
                       label = "p.format", 
                       y.position=c(42,44,46,42))+
    theme(panel.grid = element_blank(),
          axis.title.x = element_blank(),
          legend.position="none")


# Latrunculin A
lat_53bp1_data_raw <- read_xlsx("supplementary_data_1.xlsx", sheet = "53bp1_conf_latA")


lat_53bp1_data <- lat_53bp1_data_raw %>% 
    filter(Area < 1000) %>% 
    as_tibble() %>% 
    group_by(Treatment, Condition, Filename, Replicate) %>% 
    tally(name = "foci_count") %>%   
    ungroup() 

# Adding the zeros; 
# NOTE: In the Oligomycin analysis, where we manually counted the cells with no spots; with the Latrunculin A experiment performed during revisions, we updated the Fiji analysis macro to automatically create very large Area of >1000 if no spots are detected. 
lat_53bp1_data_zeros <- lat_53bp1_data_raw %>% 
    filter(Area >= 1000) %>% 
    as_tibble() %>% 
    group_by(Treatment, Condition, Filename, Replicate) %>% 
    tally(name = "foci_count") %>%   
    ungroup() %>% 
    mutate(foci_count = 0)


lat_53bp1 <- rbind(lat_53bp1_data, lat_53bp1_data_zeros) %>% 
    group_by(Treatment, Condition) %>% 
    add_tally() %>% 
    ungroup() %>% 
    mutate(labelwithN = paste0(Treatment, "\n", Condition, "\nn=", n),
           labelwithN = factor(labelwithN, levels = c("Ctrl\nsusp\nn=62",
                                                      "Ctrl\nconf\nn=89",
                                                      "LatA\nsusp\nn=61",
                                                      "LatA\nconf\nn=79"))) %>% 
    group_by(labelwithN) %>% 
    mutate(median_foci = median(foci_count),
           sd_foci = sd(foci_count)) %>% 
    ungroup() %>% 
    filter(foci_count >= (median_foci - 3 * sd_foci) & foci_count <= (median_foci + 3 * sd_foci))


stat.test.damage <- compare_means(data = lat_53bp1, 
                                  foci_count ~ labelwithN,
                                  method = 'wilcox.test') %>% 
    dplyr::filter(row_number() %in% c(1,2,5,6))

ggplot(lat_53bp1, aes(x = labelwithN, y = foci_count))+
    geom_boxplot(aes(fill=Condition),
                 colour='black',
                 width=0.8,
                 outlier.shape=NA)+
    geom_jitter(shape=21,size=3, aes(fill=Condition), width = 0.1)+
    theme_bw()+
    scale_fill_manual(values=jco)+
    scale_color_manual(values=jco)+
    stat_pvalue_manual(stat.test.damage,  
                       label = "p.format",
                       y.position=44,
                       step.increase = 0.1)+
    theme(panel.grid = element_blank(),
          axis.title.x = element_blank(),
          legend.position="none")

#BAM15, Pyr, Pyr-Oli
damage_data_bam15_raw <- read_xlsx("supplementary_data_1.xlsx", sheet = "53bp1_conf_bam15_pyr")

data_nonZero <- damage_data_bam15_raw %>% filter(Area < 500) %>% 
    group_by(Date, CellName, Treatment, Condition) %>% 
    tally() 

data_Zero <- damage_data_bam15_raw %>% filter(Area >= 500) %>% 
    group_by(Date, CellName, Treatment, Condition) %>% 
    tally() %>% 
    mutate(n = 0)

data_combine <- rbind(data_nonZero, data_Zero) %>%
    as_tibble() %>% 
    mutate()

data_combine_filtered <- data_combine %>% 
    group_by(Treatment, Condition) %>% 
    add_tally(name = "count") %>% 
    ungroup() %>% 
    mutate(across(Treatment, factor, levels = c("Control", "Pyr", "PyrOli", "BAM15")),
           across(Condition, factor, levels = c("Susp.", "Conf.")),
           labelwithN = paste(Treatment, Condition, count, sep="\n"),
           across(labelwithN, factor, levels = c("Control\nSusp.\n80",
                                                 "Control\nConf.\n82",
                                                 "BAM15\nSusp.\n75",
                                                 "BAM15\nConf.\n82",
                                                 "Pyr\nSusp.\n75",
                                                 "Pyr\nConf.\n79",
                                                 "PyrOli\nSusp.\n50",
                                                 "PyrOli\nConf.\n55"
           )))


stat.53bp1.df <- compare_means(data = data_combine_filtered, n ~ labelwithN)

ggplot(data_combine_filtered, aes(x = labelwithN, y = n))+
    geom_boxplot(outlier.shape = NA, aes(fill = Condition)) +
    geom_jitter(width=0.15, shape = 21, size = 2, aes(fill = Condition)) +
    theme_bw() +
    scale_fill_manual(values = jco) +
    scale_colour_manual(values = jco)+ 
    stat_pvalue_manual(stat.53bp1.df, label = "p.format", y.position = 40, step.increase = 0.08)

# 53BP1 post-confinement -------------------------------------
# Oligomycin 
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
    # geom_jitter(height=0.5, width=0.2, alpha=0.02, size=1)+
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

stat.ctrl.oligo <- wilcox.test(test_data$`Control\nUnconfined`,
                               test_data$`Oligomycin\nUnconfined`)

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
                       label = "p.adj",
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



# Proliferation Pyruvate --------------------------------------------------
pyru_prol_raw <- read_excel("supplementary_data_1.xlsx", sheet = "pyruvate_proliferation")

pyru_prol <- pyru_prol_raw %>% 
    group_by(Treatment, Condition, Rep, Date) %>% 
    mutate(Normalised_proliferation = Cell_count/Cell_count[Time == min(Time)])
# 
# ggplot(pyru_prol, aes(x = Time, y = Normalised_proliferation))+
#     # geom_smooth(aes(colour = Treatment, fill = Treatment, linetype = Condition))+
#     geom_path(aes(colour = Treatment, fill = Treatment, linetype = Condition))
#     theme_bw()+
#     scale_colour_manual(values = c(jco[1:3], "black"))+
#     scale_fill_manual(values = jco)

# One subset of data start at 3h. Performing scaling to adjust the hours as 1h to be minimum.
# Transform the Time column within each group
scaled_data <- pyru_prol %>%
    group_by(Treatment, Condition, Date, Rep) %>%
    mutate(Time = if (min(Time) == 3) {
        # Rescale Time when starting at 3
        round((Time - min(Time)) * (23 / (max(Time) - min(Time))) + 1)
    } else {
        # Leave Time unchanged if starting at 1
        Time
    }) %>%
    ungroup()

ggplot(scaled_data, aes(x = Time, y = Normalised_proliferation))+
    # geom_path(aes(colour = Treatment, fill = Treatment, linetype = Condition))
    geom_smooth(aes(colour = Treatment, fill = Treatment, linetype = Condition))+
    theme_bw()+
    scale_colour_manual(values = c(jco[1:3], "black"))+
    scale_fill_manual(values = jco)


# Wilcox Test on averaged data
test_data = pyru_prol %>% 
    mutate(Treatments = paste(Treatment, Condition, sep="\n")) %>% 
    group_by(Treatments, Time) %>% 
    mutate(SD=sd(Normalised_proliferation),
           Mean = mean(Normalised_proliferation)) %>% 
    ungroup() %>% 
    dplyr::select(Time, Treatments, Mean) %>% 
    unique() %>% 
    pivot_wider(names_from = Treatments, id_cols=Time, values_from = Mean)

# Initialize an empty data frame to store the results
prol_pyru_stat <- data.frame(Treatment1 = character(),
                             Treatment2 = character(),
                             p_value = numeric(),
                             stringsAsFactors = FALSE)

# Get the names of the treatment columns (columns 2 to 9)
treatment_columns <- colnames(test_data)[2:9]

# Loop over all combinations of the treatment columns
for (i in 1:(length(treatment_columns) - 1)) {
    for (j in (i + 1):length(treatment_columns)) {
        # Extract the column names for the two treatments being compared
        treatment1 <- treatment_columns[i]
        treatment2 <- treatment_columns[j]
        
        # Perform the Wilcoxon test between the two treatment columns
        stat.test <- wilcox.test(test_data[[treatment1]], test_data[[treatment2]])
        
        # Extract the p-value
        p_value <- stat.test$p.value
        
        # Append the result to prol_pyru_stat
        prol_pyru_stat <- rbind(prol_pyru_stat, 
                                data.frame(Treatment1 = treatment1,
                                           Treatment2 = treatment2,
                                           p_value = p_value,
                                           stringsAsFactors = FALSE))
    }
}

# Display the resulting data frame
print(prol_pyru_stat)


# ROS ---------------------------------------------------------------------
ros_raw <- read_excel("supplementary_data_1.xlsx", sheet = "ros") 

# data_ros_h2o2_control <- ros_raw %>% 
#     filter(Exp == "Exp1") %>% 
#     group_by(Treatment, Condition, Date, Exp, CellID) %>% 
#     summarise(meanIntDen = mean(IntDen),
#               meanRID = mean(RawIntDen), .groups = "drop") %>% 
#     group_by(Treatment, Condition) %>% 
#     add_tally() %>% 
#     mutate(labelwithN = paste(Treatment, Condition, n, sep = "\n")) %>% 
#     ungroup() %>% 
#     mutate(labelwithN = factor(labelwithN, levels = c("Control\nsusp\n17",
#                                                       "Control\nconf\n32",
#                                                       "H2O2\nsusp\n9",
#                                                       "H2O2\nconf\n15")),
#            Condition = factor(Condition, levels = c("susp", "conf")))
# 
# stat.h2o2 <- compare_means(meanIntDen ~ labelwithN, data = data_ros_h2o2_control, test = "wilcox")
# 
# ggplot(data_ros_h2o2_control, aes(x = labelwithN, y = meanIntDen)) +
#     geom_boxplot(outlier.shape = NA, aes(fill = Condition)) +
#     geom_jitter(width = 0.2, aes(fill = Condition), shape=21) +
#     theme_bw() +
#     facet_wrap(~Date, scales = "free") +
#     scale_fill_manual(values = jco) +
#     stat_pvalue_manual(stat.h2o2, label = "p.format", y.position = 4500, step.increase = 0.04)

data_ros <- ros_raw %>% 
    filter(Exp != "Exp1") %>% 
    group_by(Treatment, Condition, Date, Exp, CellID) %>% 
    summarise(meanIntDen = mean(IntDen),
              meanRID = mean(RawIntDen), .groups = "drop") %>% 
    group_by(Treatment, Condition) %>% 
    add_tally() %>% 
    mutate(labelwithN = paste(Treatment, Condition, n, sep = "\n")) %>% 
    ungroup() %>% 
    mutate(Condition = factor(Condition, levels = c("susp", "conf")),
           labelwithN = factor(labelwithN, levels = c("Control\nsusp\n54",
                                                      "Control\nconf\n65",
                                                      "Oligomycin\nsusp\n55",
                                                      "Oligomycin\nconf\n61"
           )))


data_ros_norm <- list()
for(j in data_ros$Treatment %>% unique){
    
    # Calculate mean of susp per cell line per experiment
    norm_mean = data_ros %>% filter(
        Condition == "susp",
        Treatment == j,
    ) %>%
        .$meanIntDen %>% mean(na.rm = TRUE)
    
    # Calculate sd of susp per cell line per experiment
    norm_sd = data_ros %>% filter(
        Condition == "susp",
        Treatment == j,
    ) %>%
        .$meanIntDen %>% sd(na.rm = TRUE)
    
    # Normalise
    data_ros_norm[[j]] <- data_ros %>%
        filter(
            Treatment == j
        ) %>%
        mutate(zscore = (meanIntDen - norm_mean)/norm_sd)
}

data_ros_z <- bind_rows(data_ros_norm) %>% 
    
    #Remove outside 3 sd
    group_by(Treatment, Condition) %>% 
    mutate(mean_z = mean(zscore, na.rm = TRUE),
           sd_z = sd(zscore, na.rm = TRUE)) %>% 
    ungroup() %>% 
    filter(zscore >= (mean_z - 2 * sd_z) & zscore <= (mean_z + 2 * sd_z))

stat.ros <- compare_means(zscore ~ labelwithN, data = data_ros_z, test = "wilcox")

ggplot(data_ros_z, aes(x = labelwithN, y = zscore)) +
    geom_boxplot(outlier.shape = NA, aes(fill = Condition)) +
    geom_jitter(width = 0.1, aes(fill = Condition), shape=21, size =2) +
    theme_bw() +
    facet_wrap(~Date, scales = "free")+
    scale_fill_manual(values = jco) +
    stat_pvalue_manual(stat.ros, label = "p.format", y.position = 10, step.increase = 0.06, size = 3) 


# Mitochondrial_morphology ------------------------------------------------
mito_morph_validation <- read_excel("supplementary_data_1.xlsx", sheet = "mito_morph_validation") 

data_morph <- mito_morph_validation %>%
    pivot_longer(cols = 4:16, names_to = "parameter", values_to = "measurement") %>% 
    # dplyr::filter(!grepl(".mito", parameter)) %>% 
    filter(parameter %in% c("Branch.Junctions", "Total.Branch.Length", "Sphericity..Weighted.")) 

compare_means(data = data_morph, measurement ~ Condition, group.by = "parameter") %>% view()

ggplot(data_morph, aes(x=Condition, y=measurement))+
    geom_boxplot(outlier.shape=NA)+
    geom_jitter(alpha=0.5, width=0.1)+
    facet_wrap(~parameter, scales="free")+
    theme_bw()

# CENPA_Nucleoli ----------------------------------------------------------

nucleoli_cenpa_data <- read_excel("supplementary_data_1.xlsx", sheet = "cenpa_nucleoli") %>% 
    group_by(Treatment, Condition) %>% 
    add_tally() %>% 
    mutate(Condition = factor(Condition, levels=c("Susp.", "Conf.")),
           labelwithN = paste(Treatment, Condition, n, sep="\n"),
           across(labelwithN, factor, levels = c("Control\nSusp.\n11",
                                                 "Control\nConf.\n11",
                                                 "Oligomycin\nSusp.\n12",
                                                 "Oligomycin\nConf.\n14"))
    ) %>% 
    ungroup() %>% 
    pivot_longer(cols = 8:9, names_to = "analysis", values_to = "quantification") 

stat.cenpa <- compare_means(data = nucleoli_cenpa_data, quantification ~ labelwithN, group.by = "analysis")

ggplot(data = nucleoli_cenpa_data, aes(x = labelwithN, y= quantification))+
    geom_boxplot(aes(fill =Condition), outlier.shape = NA)+
    geom_jitter(width= 0.1, size=3, shape =21, aes(fill= Condition))+
    facet_wrap(~analysis, scales="free")+
    scale_fill_manual(values = jco)+
    theme_bw() +
    stat_pvalue_manual(stat.ssd, label = "p.format", y.position = 9, step.increase = 0.06)

# Metabolomics ------------------------------------------------------------
de_novo_data_raw <-  read_excel("supplementary_data_1.xlsx", sheet = "de_novo_nucleotides_metabolomic") %>%
    pivot_longer(cols = 2:13, names_to = "sample", values_to = "rel_abundance") %>% 
    mutate(rep = str_sub(sample, -1),
           condition = case_when(grepl("WT_DMSO", sample) ~ "DMSO",
                                 grepl("WT_etoposide_0", sample) ~ "Eto_0h",
                                 grepl("WT_etoposide_8h", sample) ~ "Eto_8h",
                                 grepl("WT_etoposide_24h", sample) ~ "Eto_24h"),
           rel_abundance = as.numeric(rel_abundance, na.rm = T)) %>%
    group_by(metabolite, condition) %>% 
    summarise(mean_rel_abundance = mean(rel_abundance, na.rm = T), .groups = "drop") %>% 
    group_by(metabolite) %>% 
    mutate(mean = mean(mean_rel_abundance, na.rm = T),
           sd = sd(mean_rel_abundance, na.rm = T),
           normalised_abundance = (mean_rel_abundance - mean)/sd)
# mutate(min = min(mean_rel_abundance, na.rm = T),
#        max = max(mean_rel_abundance, na.rm = T),
#        normalised_abundance = (mean_rel_abundance - min)/(max - min))


# selected_metabolites <- c(
#     "adenine", "adenosine_cyclic_monophosphate", "adenosine_monophosphate",
#     "adenosine_triphosphate", "cytidine", "cytidine_triphosphate", 
#     "deoxy_methylthio_adenosine", "deoxyadenosine_triphosphate", 
#     "deoxycytidine_monophosphate", "deoxythymidine_triphosphate", 
#     "deoxyuridine", "flavin_adenine_dinucleotide", "guanosine", 
#     "guanosine_triphosphate", "hypoxanthine", "inosine_monophosphate", 
#     "inosine_triphosphate", "nicotinamide_adenine_dinucleotide", 
#     "orotic_acid", "uridine", "uridine_diphosphohexose", 
#     "uridine_monophosphate", "uridine_triphosphate", "xanthine", 
#     "xanthosine", "ribose_5_phosphate", "phospho_serine", 
#     "phosphoglyceric_acid", "pyridoxal_5_phosphate", "thymidine"
# )

selected_metabolites <- c(
    "inosine_monophosphate",
    "adenosine_monophosphate",    
    "uridine_monophosphate",
    "uridine_triphosphate",
    "inosine_triphosphate",
    "adenosine_triphosphate"
)


de_novo_data <- de_novo_data_raw %>% filter(metabolite %in% selected_metabolites)

de_novo_data_mat <- de_novo_data %>% 
    dplyr::select(1,2,6) %>% 
    pivot_wider(values_from = "normalised_abundance", names_from = "condition") %>%
    # filter(metabolite %in% c("adenosine_triphosphate", "adenosine_monophosphate", "cystine")) %>% 
    column_to_rownames(var = "metabolite") %>% 
    dplyr::select(DMSO, Eto_0h, Eto_8h, Eto_24h) 
# rowwise() %>% 
# mutate(DMSO = DMSO+abs(DMSO),
#        Eto_0h = Eto_0h+abs(DMSO),
#        Eto_8h = Eto_8h+abs(DMSO),
#        Eto_24h = Eto_24h+abs(DMSO))

max_val = max(max(de_novo_data_mat$DMSO, na.rm=T),
              max(de_novo_data_mat$Eto_0h, na.rm=T),
              max(de_novo_data_mat$Eto_8h, na.rm=T),
              max(de_novo_data_mat$Eto_24h, na.rm=T))

min_val = min(min(de_novo_data_mat$DMSO, na.rm=T), 
              min(de_novo_data_mat$Eto_0h, na.rm=T),
              min(de_novo_data_mat$Eto_8h, na.rm=T),
              min(de_novo_data_mat$Eto_24h, na.rm=T))

mean_val = mean(mean(de_novo_data_mat$DMSO, na.rm=T),
                mean(de_novo_data_mat$Eto_0h, na.rm=T),
                mean(de_novo_data_mat$Eto_8h, na.rm=T),
                mean(de_novo_data_mat$Eto_24h, na.rm=T))

div = (max_val-min_val)/4

col_fun3 = colorRamp2(c(max_val, 
                        # max_val-div, 
                        # min_val+div,
                        0,
                        min_val), 
                      # c("#cb4335","#fadbd8"))
                      # c("black","white","blue"))
                      c(jco[1],"white",jco[2]))

Heatmap(de_novo_data_mat %>% as.matrix(), 
        column_order = c("DMSO", "Eto_0h", "Eto_8h", "Eto_24h"),
        row_names_gp = grid::gpar(fontsize = 9),
        col = col_fun3,
        row_order = selected_metabolites)



# ATAC volcano plots ------------------------------------------------------

ctrl_atac <- read.delim("supplementary_data_4.txt") %>% 
    as_tibble() %>% 
    mutate(signif = case_when(PValue < 0.05 & logFC > 2 ~ "signif_conf",
                              PValue < 0.05 & logFC < -2 ~ "signif_susp",
                              TRUE ~ "non-signif"))


plot_atac_volcano_control <- ggplot(ctrl_atac, aes(x = logFC, y = -log10(PValue)))+
    theme_bw() +
    theme(panel.grid = element_blank())+
    geom_hline(yintercept=0, colour = "darkgray")+
    geom_vline(xintercept=0,colour = "darkgray")+
    # geom_point(alpha = 0.5, colour = "black") +
    geom_point(aes(colour = signif, size = signif, alpha= signif), shape = 16) +
    scale_colour_manual(values = c(jco[4], "black", "black"))+
    scale_size_manual(values = c(2,2,2))+
    scale_alpha_manual(values = c(0.1,0.7,0.7))

pdf("atac_volcano_ctrl_ConfvsSusp_signif_rg_black.pdf", height = 5, width = 7)
plot_atac_volcano_control
dev.off()


oligo_atac <- read.delim("supplementary_data_5.txt") %>% 
    as_tibble() %>% 
    mutate(signif = case_when(PValue < 0.05 & logFC > 2 ~ "signif_conf",
                              PValue < 0.05 & logFC < -2 ~ "signif_susp",
                              TRUE ~ "non-signif")) 
# %>% 
#     sample_n(10000)

plot_atac_volcano_oligo <- ggplot(oligo_atac, aes(x = logFC, y = -log10(PValue)))+
    theme_bw() +
    theme(panel.grid = element_blank())+
    geom_hline(yintercept=0, colour = "darkgray")+
    geom_vline(xintercept=0,colour = "darkgray")+
    # geom_point(alpha = 0.5, colour = "black") +
    geom_point(aes(colour = signif, size = signif, alpha= signif), shape = 16) +
    scale_colour_manual(values = c(jco[4], "black", "black"))+
    scale_size_manual(values = c(2,3,3))+
    scale_alpha_manual(values = c(0.1,0.7,0.7))

pdf("atac_volcano_oligo_ConfvsSusp_signif_rg_black.pdf", height = 5, width = 7)
plot_atac_volcano_oligo
dev.off()




# ATAC significant genes ORA ----------------------------------------------
ctrl_df <- read_excel("supplementary_data_1.xlsx", sheet = "CU_vs_UU_diff_exp_annotated_pea")
oligo_df <- read_excel("supplementary_data_1.xlsx", sheet =  "CT_vs_UT_diff_exp_annotated_pea")

background_genes = rbind(ctrl_df, oligo_df) %>% 
    select(ENSEMBL) %>% 
    unique() %>% 
    .$ENSEMBL

#Ctrl
ctrl_df_signif_pos <- ctrl_df %>% 
    filter(PValue < 0.05,
           logFC >= 2) %>% 
    dplyr::select(ENSEMBL) %>% 
    unique()

ctrl_df_signif_neg <- ctrl_df %>% 
    filter(PValue < 0.05,
           logFC <= -2) %>% 
    dplyr::select(ENSEMBL) %>% 
    unique()

ctrl_df_signif_both <- ctrl_df %>% 
    filter(PValue < 0.05,
           abs(logFC) >= 2) %>% 
    dplyr::select(ENSEMBL) %>% 
    unique()

ctrl_signif_ora <- enrichGO(
    gene          = ctrl_df_signif_both$ENSEMBL,  # Vector of genes
    OrgDb         = org.Hs.eg.db,        # The database for human genes
    keyType       = "ENSEMBL",           # The identifier type (UniProt IDs)
    ont           = "ALL",               # Ontology (ALL includes BP, CC, and MF)
    pAdjustMethod = "BH",                # p-value adjustment method
    pvalueCutoff  = 0.05,                # p-value threshold
    qvalueCutoff  = 0.05,                # q-value threshold for multiple testing correction
    minGSSize     = 10,                  # Minimum gene set size
    maxGSSize     = 500,                  # Maximum gene set size
    universe    = background_genes
)

ggplot(ctrl_signif_ora, aes(x = reorder(Description, Count), y = Count, fill = pvalue)) +
    geom_bar(stat = "identity") +
    coord_flip() +  # Flip the coordinates to make the bars horizontal
    scale_fill_gradient(low = jco[2], high = jco[1], name = "p-value") +
    theme_bw()


#Oligo
oligo_df_signif_pos <- oligo_df %>% 
    filter(PValue < 0.05,
           logFC >= 2) %>% 
    dplyr::select(ENSEMBL) %>% 
    unique()

oligo_df_signif_neg <- oligo_df %>% 
    filter(PValue < 0.05,
           logFC <= -2) %>% 
    dplyr::select(ENSEMBL) %>% 
    unique()

oligo_df_signif_both <- oligo_df %>% 
    filter(PValue < 0.05,
           abs(logFC) >= 2) %>% 
    dplyr::select(ENSEMBL) %>% 
    unique()

oligo_signif_ora <- enrichGO(
    gene          = oligo_df_signif_both$ENSEMBL,  # Vector of genes
    OrgDb         = org.Hs.eg.db,        # The database for human genes
    keyType       = "ENSEMBL",           # The identifier type (UniProt IDs)
    ont           = "ALL",               # Ontology (ALL includes BP, CC, and MF)
    pAdjustMethod = "BH",                # p-value adjustment method
    pvalueCutoff  = 0.05,                # p-value threshold
    qvalueCutoff  = 0.05,                # q-value threshold for multiple testing correction
    minGSSize     = 10,                  # Minimum gene set size
    maxGSSize     = 500,                  # Maximum gene set size
    universe    = background_genes
)



# Using venn diagrams to check gene overlap


# df <- oligo_signif_ora@result %>% 
#     separate_rows(geneID, sep = "/") 

# gene_sets <- list(
#     "postsynaptic membrane" = unlist(strsplit(df$geneID[df$Description == "postsynaptic membrane"], "/")),
#     "intrinsic component of synaptic membrane" = unlist(strsplit(df$geneID[df$Description == "intrinsic component of synaptic membrane"], "/")),
#     # "leading edge membrane" = unlist(strsplit(df$geneID[df$Description == "leading edge membrane"], "/")),
#     "intrinsic component of postsynaptic membrane" = unlist(strsplit(df$geneID[df$Description == "intrinsic component of postsynaptic membrane"], "/")),
#     "integral component of postsynaptic membrane" = unlist(strsplit(df$geneID[df$Description == "integral component of postsynaptic membrane"], "/"))
#     # "intrinsic component of presynaptic membrane" = unlist(strsplit(df$geneID[df$Description == "intrinsic component of presynaptic membrane"], "/")),
#     # "ion gated channel activity" = unlist(strsplit(df$geneID[df$Description == "ion gated channel activity"], "/"))
# )
# 
# venn.plot <- venn.diagram(
#     x = gene_sets,
#     category.names = names(gene_sets),
#     filename = NULL,  # Display the Venn plot in R directly
#     fill = c("red", "blue", "green", "yellow"),
#     alpha = 0.5,
#     cex = 1.5,
#     cat.cex = 1.2,
#     cat.pos = 0
# )
# 
# grid.draw(venn.plot)

oligo_ora_results <- oligo_signif_ora@result %>% 
    separate_rows(geneID, sep = "/") %>% 
    filter(!Description %in% c("intrinsic component of postsynaptic membrane",
                               "integral component of postsynaptic membrane")) %>% 
    select(-geneID) %>% 
    unique()


ggplot(oligo_ora_results, aes(x = factor(Description, levels = c(
    "ion gated channel activity",
    "intrinsic component of presynaptic membrane",
    "intrinsic component of synaptic membrane",
    "leading edge membrane",
    "postsynaptic membrane")),
    y = Count, fill = pvalue)) +
    geom_bar(stat = "identity") +
    coord_flip() +  # Flip the coordinates to make the bars horizontal
    scale_fill_gradient(low = jco[2], high = jco[1], name = "p-value") +
    theme_bw()


# Shortlisted GO heatmaps ----------------------------------------------------------
# read in the diff exp data frames
ctrl_diff <- read_excel("supplementary_data_1.xlsx", sheet = "cu_vs_uu_diff_exp_collapsed_ann")
oligo_diff <- read_excel("supplementary_data_1.xlsx", sheet =  "ct_vs_ut_diff_exp_collapsed_ann")

# combine the diff exp data frames
logFC_df <- left_join(
    ctrl_diff %>% select(Ensemble,Symbol, Identity,  logFC) %>% rename(logFC_ctrl = logFC),
    oligo_diff %>% select(Ensemble, Symbol, Identity, logFC) %>% rename(logFC_oligo = logFC),
    by = c("Ensemble", "Symbol", "Identity")
) %>% 
    unique() %>% 
    rename(SYMBOL = Symbol,
           Annotation = Identity)

# GOs of interest
# Updated list of GO terms
go_data <- oligo_ora_results %>% 
    select(ID, Description) %>% 
    set_names("GOALL", "go_name")

reference_go_data <- data.frame(GOALL = c("GO:0007049","GO:0030036"),
                                go_name = c("cell cycle", "actin cytoskeleton organisation")) %>% 
    as_tibble()

# Query Ensembl for genes associated with these GO terms
genes_by_go <- AnnotationDbi::select(org.Hs.eg.db, keytype="GOALL", keys=go_data$GOALL, columns=c("SYMBOL"))

reference_genes_by_go <- AnnotationDbi::select(org.Hs.eg.db, keytype="GOALL", keys=reference_go_data$GOALL, columns=c("SYMBOL")) %>% 
    left_join(reference_go_data)

genes_by_go_complete <- left_join(genes_by_go, go_data) %>%
    as_tibble() %>%
    dplyr::select(1,4,5) %>%
    unique()

# dir = "C:/Users/rghose/OneDrive - CRG - Centre de Regulacio Genomica/Rito_Fabio/Submission/NCB_Revision/Revision_Data/ATAC/Ilir_NoBatchEffect_Final/genes_collapsed_regions/significant_go_genes/"
# Define a function to get the Entrez ID and summary for a gene symbol
get_gene_summary <- function(gene_symbol) {
    # Search Entrez for the gene symbol in the 'gene' database
    search_result <- entrez_search(db = "gene", term = paste0(gene_symbol, "[symbol] AND Homo sapiens[orgn]"))
    
    # If the gene was found, fetch its summary
    if (length(search_result$ids) > 0) {
        entrez_id <- search_result$ids[1]  # Take the first result (most relevant)
        gene_summary <- entrez_summary(db = "gene", id = entrez_id)
        return(gene_summary$summary)
    } else {
        return(NA)  # Return NA if the gene symbol was not found
    }
}

for(i in go_data$GOALL){
    
    go_description = go_data %>% filter(GOALL == i) %>% .$go_name
    
    go_genes = genes_by_go_complete %>% filter(GOALL == i) %>% 
        set_names("GOALL", "SYMBOL", "go_id")
    
    # logFC df select GO genes
    df_heatmap_logfc = left_join(go_genes, logFC_df) %>%
        filter(Annotation == "Promoter") %>% 
        na.omit() %>%
        select(SYMBOL, logFC_ctrl, logFC_oligo) %>%
        unique() 
    
    actin_ref_genes <- reference_genes_by_go %>% filter(go_name == "actin cytoskeleton organisation") %>% .$SYMBOL %>% unique()
    cellcycle_ref_genes <- reference_genes_by_go %>% filter(go_name == "cell cycle") %>% .$SYMBOL %>% unique()
    
    
    gene_annotation <- df_heatmap_logfc %>%
        mutate(cell_cycle = case_when(SYMBOL %in% cellcycle_ref_genes ~ "yes", TRUE ~ "no"),
               actin_regulation = case_when(SYMBOL %in% actin_ref_genes ~ "yes", TRUE ~ "no"))
    
    # Define colors for the yes/no values
    yes_no_colors <- c("yes" = "black", 
                       "no"  = "gray")  
    
    # Create the row annotations for dna_damage, cell_cycle, and actin_regulation
    row_ha <- rowAnnotation(
        `Cell Cycle` = gene_annotation$cell_cycle,
        `Actin Regulation` = gene_annotation$actin_regulation,
        col = list(
            `DNA Damage` = yes_no_colors,
            `Cell Cycle` = yes_no_colors,
            `Actin Regulation` = yes_no_colors
        )
    )
    
    
    # col_fun <- colorRamp2(
    #     c(min(range(df_heatmap_logfc, na.rm = TRUE)), 
    #       0, 
    #       max(range(df_heatmap_logfc, na.rm = TRUE))), 
    #     c(jco[2], "white", jco[1]))
    col_fun <- colorRamp2(
        c(-1.5, -0.75, 0, 0.75, 1.5), 
        c("#0c4a7b", jco[2], "white", jco[1], "#faa61b"))
    
    
    
    
    # Create heatmap for this annotation
    ht_logFC_promoter <- Heatmap(df_heatmap_logfc %>% column_to_rownames(var="SYMBOL"),
                                 cluster_columns = F,
                                 # cluster_rows =F,
                                 name = "logFC",
                                 column_title = go_description,
                                 col = col_fun,
                                 row_names_gp = gpar(fontsize = 5),
                                 right_annotation = row_ha)
    
    # pdf(paste0(dir,go_description,"_logFC_promoter_heatmap.pdf"))
    draw(ht_logFC_promoter)
    # dev.off()
}


# ATAC_Seq Chromosome Heatmap -----------------------------------------------------------------

# Release 37.1 of the human genome sequence includes multiple assemblies of the chromosomes from different source data. The tables on chromosome size and base composition are based solely on the 24 reference human chromosome sequences. Source: http://www.cshlp.org/ghg5_all/section/dna.shtml

chr_gene_density <- data.frame(
    seqnames = c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", 
                 "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", 
                 "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", 
                 "chrX"),
    Size_Mb = c(249.3, 243.2, 198.0, 191.2, 180.9, 171.1, 159.1, 146.4, 141.2, 
                135.5, 135.0, 133.9, 115.2, 107.3, 102.5, 90.4, 81.2, 78.1, 59.1, 
                63.0, 48.1, 51.3, 155.3),
    Sequenced_Mb = c(225.3, 238.2, 194.8, 187.7, 177.7, 167.4, 155.4, 142.9, 120.1, 
                     131.3, 131.1, 130.5, 95.6, 88.3, 81.7, 78.9, 77.8, 74.7, 55.8, 
                     59.5, 35.1, 34.9, 151.1),
    Genes = c(1959, 1184, 1029, 721, 835, 1002, 855, 638, 748, 714, 1236, 987, 
              305, 577, 547, 783, 1111, 257, 1332, 518, 213, 418, 806),
    Genes_per_Mb = c(7.86, 4.87, 5.20, 3.77, 4.62, 5.86, 5.37, 4.36, 5.30, 5.27, 
                     9.16, 7.37, 2.65, 5.37, 5.33, 8.67, 13.68, 3.29, 22.53, 8.22, 
                     4.43, 8.15, 5.19),
    Genes_per_Sequenced_Mb = c(8.70, 4.97, 5.28, 3.84, 4.70, 5.99, 5.50, 4.47, 6.23, 
                               5.44, 9.43, 7.56, 3.19, 6.54, 6.70, 9.93, 14.28, 
                               3.44, 23.87, 8.71, 6.07, 11.98, 5.33)
)




ctrl_atac <- read.delim("supplementary_data_4.txt") %>% 
    as_tibble() %>% 
    mutate(signif = case_when(PValue < 0.05 & abs(logFC) > 2 ~ "signif",
                              TRUE ~ "non-signif")) %>% 
    filter(seqnames != "chrY")

# Merge the chromosome sizes with the atac-seq data to ensure proper limits
ctrl_atac <- merge(ctrl_atac, chr_gene_density %>% select(1, 2), by = "seqnames")

ctrl_mat <- ctrl_atac %>% 
    mutate(peak_width = end-start) %>% 
    group_by(seqnames, direction, Size_Mb) %>%
    summarise(total_peak_width = sum(peak_width), .groups = "drop") %>% 
    mutate(normalised_peaks = total_peak_width/(Size_Mb*1000000),
           enriched_in = case_when(direction == "down" ~ "Control Susp.",
                                   direction == "up" ~"Control Conf.")) %>% 
    select(normalised_peaks, seqnames, enriched_in) %>% 
    pivot_wider(names_from = "seqnames", values_from = "normalised_peaks") %>% 
    column_to_rownames(var = "enriched_in")


# Oligo
oligo_atac <- read.delim("supplementary_data_5.txt") %>% 
    as_tibble() %>% 
    mutate(signif = case_when(PValue < 0.05 & abs(logFC) > 2 ~ "signif",
                              TRUE ~ "non-signif"))%>% 
    filter(seqnames != "chrY")

# Merge the chromosome sizes with the atac-seq data to ensure proper limits
oligo_atac <- merge(oligo_atac, chr_gene_density %>% select(1, 2), by = "seqnames")


oligo_mat <- oligo_atac %>%
    mutate(peak_width = end-start) %>% 
    group_by(seqnames, direction, Size_Mb) %>%
    summarise(total_peak_width = sum(peak_width), .groups = "drop") %>% 
    mutate(normalised_peaks = total_peak_width/(Size_Mb*1000000),
           enriched_in = case_when(direction == "down" ~ "Oligo Susp.",
                                   direction == "up" ~"Oligo Conf.")) %>% 
    select(normalised_peaks, seqnames, enriched_in) %>% 
    pivot_wider(names_from = "seqnames", values_from = "normalised_peaks") %>% 
    column_to_rownames(var = "enriched_in")

# Heatmap non ratio

comb_mat <- t(rbind(ctrl_mat, oligo_mat))
comb_mat <- comb_mat[chr_gene_density$seqnames,]


col_fun <- colorRamp2(
    c(min(range(comb_mat, na.rm = TRUE)), mean(range(comb_mat, na.rm = TRUE)), max(range(comb_mat, na.rm = TRUE))), 
    c("#ffeeee", "#a65a82", "#00008b")
)

set.seed(100)
Chr_Density_anno <- rowAnnotation(`Chr_Density` = chr_gene_density$Genes_per_Mb)

Chr_SizeMb_anno <- rowAnnotation(`Chr_Size\n(Mb)` = chr_gene_density$Size_Mb)


Heatmap(comb_mat, 
        cluster_columns = F, 
        cluster_rows = F, 
        name = "Normalised\npeak count", col = col_fun, 
        column_title = "ATAC-Seq Peak count normalised by Chromosome Length",
        right_annotation = c(Chr_Density_anno, Chr_SizeMb_anno))


# Heatmap ratio

# perform the ratio
comb_mat_norm <- log10(cbind(comb_mat[,2]/comb_mat[,1],comb_mat[,4]/comb_mat[,3]))
#rename ratio columns
colnames(comb_mat_norm) <- c("Control\nratio", "Oligo\nratio")

#colour scale for the ratio
col_fun = colorRamp2(c(max(range(comb_mat_norm)), 
                       mean(range(comb_mat_norm)), 
                       min(range(comb_mat_norm))), 
                     c("white", "#0a7acc", "darkblue"))

#colour scale for gene density per chromosome
col_fun2 = colorRamp2(c(max(chr_gene_density$Genes_per_Mb),
                        mean(chr_gene_density$Genes_per_Mb),
                        min(chr_gene_density$Genes_per_Mb)), 
                      c( "darkgreen","lightgreen", "white"))

#colour scale for chromosome size
col_fun3 = colorRamp2(c(max(chr_gene_density$Size_Mb), 
                        mean(chr_gene_density$Size_Mb), 
                        min(chr_gene_density$Size_Mb)), 
                      c("red3", "#EFCB68", "#0a7acc"))


Chr_anno <- HeatmapAnnotation(`Chr_Density` = chr_gene_density$Genes_per_Mb,
                              `Chr_Size` = chr_gene_density$Size_Mb,
                              col = list(`Chr_Density` = col_fun2, 
                                         `Chr_Size` = col_fun3))


# set.seed(100)
Heatmap(t(comb_mat_norm), 
        cluster_columns = F, 
        cluster_rows = F, 
        name = "Normalised\npeak count (log10)", col = col_fun, 
        column_title = "ATAC-Seq Peak count normalised by Chromosome Length",
        bottom_annotation = Chr_anno)


# ATAC-Seq Chromomsome Circos plots ATAC-Seq ---------------------------------------------------

# Define the color gradient function based on the matrix values
col_fun <- colorRamp2(
    c(min(range(comb_mat, na.rm = TRUE)), mean(range(comb_mat, na.rm = TRUE)), max(range(comb_mat, na.rm = TRUE))), 
    c("#ffeeee", "#a65a82", "#00008b")
)

# Circos Control
# Calculate the midpoint between start and end
ctrl_atac <- read.delim("CU_vs_UU_csaw_DA_windows_all.txt") %>% 
    as.data.frame() %>% 
    as_tibble() %>% 
    mutate(signif = case_when(PValue < 0.05 & abs(logFC) > 2 ~ "signif",
                              TRUE ~ "non-signif")) %>% 
    filter(seqnames != "chrY") %>%
    # rowwise() %>%
    mutate(location = (start + end) / 2)  # Calculate midpoint

# Bin size: 100,000
# bin_size <- 10000000


# Plot the tracks
pdf("circos_control.pdf", width = 5, height = 5)  # Open a PDF device to save the plot

# Initialize the circos plot with ideograms
circos.clear()  # Clear any previous plots
circos.par("start.degree" = 90, 
           "gap.degree" = rep(c(2, 4), 12), 
           "track.margin" = c(0.01, 0))  # Customize gaps between chromosomes
circos.initializeWithIdeogram(species = "hg19")  # Initialize with human ideogram (hg19)

# Plot lines for Positive logFC values using circos.segments() and coloring from comb_mat[,1:2]
circos.par("track.height" = 0.3)

# Plot the positive logFC values
circos.track(ylim = c(0, max(ctrl_atac$logFC, na.rm = TRUE)), panel.fun = function(x, y) {
    # Get the current chromosome/sector being processed
    sector.index = CELL_META$sector.index
    
    # Filter the data for the current chromosome/sector
    chr_data <- ctrl_atac[ctrl_atac$seqnames == sector.index, ]
    
    # Only keep positive logFC values
    chr_data$logFC_pos <- ifelse(chr_data$logFC > 0, chr_data$logFC, 0)
    
    # Get the 'Control Conf.' values from comb_mat[,1:2]
    if (sector.index %in% rownames(comb_mat[,1:2])) {
        chr_color <- col_fun(comb_mat[,1:2][sector.index, "Control Conf."])
    } else {
        chr_color <- "grey"
    }
    
    if (nrow(chr_data) > 0) {
        # Plot lines for positive logFC values
        for (i in 1:nrow(chr_data)) {
            circos.segments(chr_data$location[i], 0, chr_data$location[i], chr_data$logFC_pos[i], 
                            col = chr_color, lwd = 2)  # Color from comb_mat[,1:2]
        }
    }
}, bg.border = NA)

# Plot lines for Negative logFC values using circos.segments() and coloring from comb_mat[,1:2]
circos.par("track.height" = 0.3)

# Plot the negative logFC values
circos.track(ylim = c(min(ctrl_atac$logFC, na.rm = TRUE), 0), panel.fun = function(x, y) {
    # Get the current chromosome/sector being processed
    sector.index = CELL_META$sector.index
    
    # Filter the data for the current chromosome/sector
    chr_data <- ctrl_atac[ctrl_atac$seqnames == sector.index, ]
    
    # Only keep negative logFC values
    chr_data$logFC_neg <- ifelse(chr_data$logFC < 0, chr_data$logFC, 0)
    
    # Get the 'Control Susp.' values from comb_mat[,1:2]
    if (sector.index %in% rownames(comb_mat[,1:2])) {
        chr_color <- col_fun(comb_mat[,1:2][sector.index, "Control Susp."])
    } else {
        chr_color <- "grey"
    }
    
    if (nrow(chr_data) > 0) {
        # Plot lines for negative logFC values
        for (i in 1:nrow(chr_data)) {
            circos.segments(chr_data$location[i], chr_data$logFC_neg[i], chr_data$location[i], 0, 
                            col = chr_color, lwd = 2)  # Color from comb_mat[,1:2]
        }
    }
}, bg.border = NA)


# Close the PDF device to save the plot
dev.off()




# Circos Oligomycin Plot
oligo_atac <- read.delim("CT_vs_UT_csaw_DA_windows_all.txt") %>% 
    as.data.frame() %>% 
    as_tibble() %>% 
    mutate(signif = case_when(PValue < 0.05 & abs(logFC) > 2 ~ "signif", TRUE ~ "non-signif")) %>% 
    filter(seqnames != "chrY") %>%
    mutate(location = (start + end) / 2) # Calculate midpoint

# Plot the tracks
pdf("circos_oligo.pdf", width = 5, height = 5)  # Open a PDF device to save the plot

# Initialize the circos plot with ideograms
circos.clear()  # Clear any previous plots
circos.par("start.degree" = 90, 
           "gap.degree" = rep(c(2, 4), 12), 
           "track.margin" = c(0.01, 0))  # Customize gaps between chromosomes
circos.initializeWithIdeogram(species = "hg19")  # Initialize with human ideogram (hg19)



# Plot lines for Positive logFC values using circos.segments() and coloring from comb_mat[,3:4]
circos.par("track.height" = 0.3)

# Plot the positive logFC values
circos.track(ylim = c(0, max(oligo_atac$logFC, na.rm = TRUE)), panel.fun = function(x, y) {
    sector.index = CELL_META$sector.index
    
    # Filter the data for the current chromosome/sector
    chr_data <- oligo_atac[oligo_atac$seqnames == sector.index, ]
    
    # Only keep positive logFC values
    chr_data$logFC_pos <- ifelse(chr_data$logFC > 0, chr_data$logFC, 0)
    
    # Get the 'Oligo Conf.' values from comb_mat[,3:4]
    if (sector.index %in% rownames(comb_mat[,3:4])) {
        chr_color <- col_fun(comb_mat[,3:4][sector.index, "Oligo Conf."])
    } else {
        chr_color <- "grey"
    }
    
    if (nrow(chr_data) > 0) {
        # Plot lines for positive logFC values
        for (i in 1:nrow(chr_data)) {
            circos.segments(chr_data$location[i], 0, chr_data$location[i], chr_data$logFC_pos[i], 
                            col = chr_color, lwd = 2)  # Color from comb_mat[,3:4]
        }
    }
}, bg.border = NA)

# Plot lines for Negative logFC values using circos.segments() and coloring from comb_mat[,3:4]
circos.par("track.height" = 0.3)

# Plot the negative logFC values
circos.track(ylim = c(min(oligo_atac$logFC, na.rm = TRUE), 0), panel.fun = function(x, y) {
    sector.index = CELL_META$sector.index
    
    # Filter the data for the current chromosome/sector
    chr_data <- oligo_atac[oligo_atac$seqnames == sector.index, ]
    
    # Only keep negative logFC values
    chr_data$logFC_neg <- ifelse(chr_data$logFC < 0, chr_data$logFC, 0)
    
    # Get the 'Oligo Susp.' values from comb_mat[,3:4]
    if (sector.index %in% rownames(comb_mat[,3:4])) {
        chr_color <- col_fun(comb_mat[,3:4][sector.index, "Oligo Susp."])
    } else {
        chr_color <- "grey"
    }
    
    if (nrow(chr_data) > 0) {
        # Plot lines for negative logFC values
        for (i in 1:nrow(chr_data)) {
            circos.segments(chr_data$location[i], chr_data$logFC_neg[i], chr_data$location[i], 0, 
                            col = chr_color, lwd = 2)  # Color from comb_mat[,3:4]
        }
    }
}, bg.border = NA)

# Close the PDF device to save the plot
dev.off()




