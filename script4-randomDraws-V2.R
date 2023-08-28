library("dplyr")
library("ggplot2")
library("epitools")
library("descr")
library("ggstance")
library("ggpubr")
adpa <- read.csv("data_phase 3/ADPA.csv")
adsl <- read.csv("data_phase 3/ADSL.csv")
#####################################################
############### MCAR ################################
#####################################################
############### calculate OR via glm ################
MCAR_get_OR <- function(SEED,percentage){
  set.seed(SEED)
  subjid_sample<-sample(1:1831,1831*percentage)
  ##### completers ######
  completersData <- adpa %>%
    filter(!(SUBJID %in% subjid_sample) & PARAMN == 10 & !(AVISITN == 2)) %>%
    left_join(select(adsl,SUBJID,SEX),by="SUBJID") %>%
    select(c(1:20,29)) %>%
    filter(AVISITN == 6)
  # threeWayTableCompleters<-ftable(table(completersData$TRTPN,completersData$SEX,completersData$PCHGCA1N))
  model1 <- glm(PCHGCA1N~SEX+TRTP, family = binomial(),data = completersData)
  coef.model1 <- exp(coefficients(model1))
  ORcompletersSexM2F <- as.numeric(round(coef.model1[2],4)) # OR male to female
  ORcompletersAvsP <- as.numeric(round(1/coef.model1[3],4)) # OR active control vs placebo
  ORcompleters140vsP <- as.numeric(round(coef.model1[4]/coef.model1[3],4)) # OR 140mg vs placebo
  ORcompleters210vsP <- as.numeric(round(coef.model1[5]/coef.model1[3],4)) # OR 140mg vs placebo
  ORcompleters <- list(SexM2F = ORcompletersSexM2F,AvsP = ORcompletersAvsP,T140vsP = ORcompleters140vsP, T210vsP =ORcompleters210vsP)
  ##### imputed ######
  imputedData <- adpa %>%
    filter(PARAMN == 10 & !(AVISITN == 2)) %>%
    mutate(PCHGCA1N = ifelse(SUBJID %in% subjid_sample,0,PCHGCA1N))%>%
    left_join(select(adsl,SUBJID,SEX),by="SUBJID")  %>%
    select(c(1:20,29)) %>%
    filter(AVISITN == 6)
  # threeWayTableImputed<-ftable(table(imputedData$TRTPN,imputedData$SEX,imputedData$PCHGCA1N))
  model2 <- glm(PCHGCA1N~SEX+TRTP, family = binomial(),data = imputedData)
  coef.model2 <- exp(coefficients(model2))
  ORimputedSexM2F <- as.numeric(round(coef.model2[2],4)) # OR male to female
  ORimputedAvsP <- as.numeric(round(1/coef.model2[3],4)) # OR active control vs placebo
  ORimputed140vsP <- as.numeric(round(coef.model2[4]/coef.model2[3],4)) # OR 140mg vs placebo
  ORimputed210vsP <- as.numeric(round(coef.model2[5]/coef.model2[3],4)) # OR 140mg vs placebo
  ORimputed <- list(SexM2F = ORimputedSexM2F,AvsP = ORimputedAvsP,T140vsP = ORimputed140vsP, T210vsP =ORimputed210vsP)
  ##### LOCF #####
  imputedLOCFData <- adpa %>%
    filter(PARAMN == 10 & !(AVISITN == 2)) %>%
    group_by(SUBJID) %>%
    mutate(LAGPCHGCA1N = lag(PCHGCA1N,n=1,order_by = SUBJID)) %>%
    left_join(select(adsl,SUBJID,SEX),by="SUBJID") %>%
    mutate(PCHGCA1N = ifelse(SUBJID %in% subjid_sample,LAGPCHGCA1N,PCHGCA1N)) %>%
    select(c(1:20,29,30)) %>%
    filter(AVISITN == 6)
  model3 <- glm(PCHGCA1N~SEX+TRTP, family = binomial(),data = imputedLOCFData)
  coef.model3 <- exp(coefficients(model3))
  ORimputedSexM2F <- as.numeric(round(coef.model3[2],4)) # OR male to female
  ORimputedAvsP <- as.numeric(round(1/coef.model3[3],4)) # OR active control vs placebo
  ORimputed140vsP <- as.numeric(round(coef.model3[4]/coef.model3[3],4)) # OR 140mg vs placebo
  ORimputed210vsP <- as.numeric(round(coef.model3[5]/coef.model3[3],4)) # OR 140mg vs placebo
  ORlocf <- list(SexM2F = ORimputedSexM2F,AvsP = ORimputedAvsP,T140vsP = ORimputed140vsP, T210vsP =ORimputed210vsP)
  ##########
  # return #
  c(list(ORcompleters = ORcompleters,ORimputed=ORimputed,ORlocf = ORlocf))
}
n = 1000 # take seeds 1 to 1000
############# 10% missing ##################
MCAR_10P_result<-lapply(1:1000,
                        FUN = MCAR_get_OR,
                        0.1)
tmp.1<-unlist(MCAR_10P_result)
MCAR_10P_meanOR<- c(by(tmp.1, names(tmp.1), mean))
MCAR_10P_varOR<- c(by(tmp.1, names(tmp.1), var))
MCAR_10P_result <- matrix(,nrow=12,ncol = 4)
MCAR_10P_result[,1]<-MCAR_10P_meanOR
MCAR_10P_result[,2]<-MCAR_10P_varOR
# lower CI = mean - phi(0.975) * sd/sqrt(n)
MCAR_10P_result[,3]<-MCAR_10P_meanOR-qnorm(0.975) * sqrt(MCAR_10P_varOR)/sqrt(n)
# upper CI = mean + phi(0.975) * sd/sqrt(n)
MCAR_10P_result[,4]<-MCAR_10P_meanOR+qnorm(0.975) * sqrt(MCAR_10P_varOR)/sqrt(n)
MCAR_10P_result = data.frame(MCAR_10P_result)
names(MCAR_10P_result)<- c("mean","var","low","high")
row.names(MCAR_10P_result)<-names(MCAR_10P_meanOR)
# get rid of sex
MCAR_10P_result1 <-MCAR_10P_result[c(-2,-6,-10),]
MCAR_10P_result1$var <- round(MCAR_10P_result1$var,2)
MCAR_10P_result1$OR = factor(c("Active Control vs Placebo","140 mg Test Drug vs Placebo","210 mg Test Drug vs Placebo","Active Control vs Placebo","140 mg Test Drug vs Placebo","210 mg Test Drug vs Placebo","Active Control vs Placebo","140 mg Test Drug vs Placebo","210 mg Test Drug vs Placebo"))
MCAR_10P_result1$dataHandling = factor(c("completers","completers","completers","composite","composite","composite","LOCF","LOCF","LOCF"))
MCAR10PvarPlot <- ggplot(MCAR_10P_result1, aes(fill=OR, y=var, x=dataHandling)) + 
  geom_bar(position="dodge", stat="identity")+theme_bw()+
  theme(panel.grid.minor = element_blank(),panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"))+
  theme(plot.title = element_text(hjust = 0.5,size=13))+
  ylab("Var") +
  xlab("Strategies") +
  ggtitle("10% Missing, MCAR")+ 
  scale_color_discrete(
    guide = guide_legend(override.aes = list(alpha=0))
  ) + theme(
    legend.text = element_text(color = "white"),
    legend.title = element_text(color = "white"),
    legend.key = element_rect(fill = "white")
  )+ylim(0,162)+
  geom_text(aes(x = dataHandling, y = var,label=var),size=3,vjust=-0.5,position = position_dodge(width=0.9))
MCAR10PvarPlot 
########### 20% missing ###################
MCAR_20P_result<-lapply(1:1000,
                        FUN = MCAR_get_OR,
                        0.2)
tmp.2<-unlist(MCAR_20P_result)
MCAR_20P_meanOR<- c(by(tmp.2, names(tmp.2), mean))
MCAR_20P_varOR<- c(by(tmp.2, names(tmp.2), var))
MCAR_20P_result <- matrix(,nrow=12,ncol = 4)
MCAR_20P_result[,1]<-MCAR_20P_meanOR
MCAR_20P_result[,2]<-MCAR_20P_varOR
MCAR_20P_result[,3]<-MCAR_20P_meanOR-qnorm(0.975) * sqrt(MCAR_20P_varOR)/sqrt(n)
MCAR_20P_result[,4]<-MCAR_20P_meanOR+qnorm(0.975) * sqrt(MCAR_20P_varOR)/sqrt(n)
MCAR_20P_result = data.frame(MCAR_20P_result)
names(MCAR_20P_result)<- c("mean","var","low","high")
row.names(MCAR_20P_result)<-names(MCAR_20P_meanOR)
# Var plot #
# get rid of sex
MCAR_20P_result1 <-MCAR_20P_result[c(-2,-6,-10),]
MCAR_20P_result1$var <- round(MCAR_20P_result1$var,2)
MCAR_20P_result1$OR = factor(c("Active Control vs Placebo","140 mg Test Drug vs Placebo","210 mg Test Drug vs Placebo","Active Control vs Placebo","140 mg Test Drug vs Placebo","210 mg Test Drug vs Placebo","Active Control vs Placebo","140 mg Test Drug vs Placebo","210 mg Test Drug vs Placebo"))
MCAR_20P_result1$dataHandling = factor(c("completers","completers","completers","composite","composite","composite","LOCF","LOCF","LOCF"))
MCAR20PvarPlot <- ggplot(MCAR_20P_result1, aes(fill=OR, y=var, x=dataHandling)) + 
  geom_bar(position="dodge", stat="identity")+theme_bw()+
  theme(panel.grid.minor = element_blank(),panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"))+
  theme(plot.title = element_text(hjust = 0.5,size=13))+
  ylab("Var") +
  xlab("Strategies") +
  ggtitle("20% Missing, MCAR")+ 
  scale_color_discrete(
    guide = guide_legend(override.aes = list(color = "white"))
  ) + theme(
    legend.text = element_text(color = "white"),
    legend.title = element_text(color = "white"),
    legend.key = element_rect(fill = "white")
  )+ ylim(0,162)+
  geom_text(aes(x = dataHandling, y = var,label=var),size=3,vjust=-0.5,position = position_dodge(width=0.9))
MCAR20PvarPlot 
########### 30% missing ###################
MCAR_30P_result<-lapply(1:1000,
                        FUN = MCAR_get_OR,
                        0.3)
tmp.3<-unlist(MCAR_30P_result)
MCAR_30P_meanOR<- c(by(tmp.3, names(tmp.3), mean))
MCAR_30P_varOR<- c(by(tmp.3, names(tmp.3), var))
MCAR_30P_result <- matrix(,nrow=12,ncol = 4)
MCAR_30P_result[,1]<-MCAR_30P_meanOR
MCAR_30P_result[,2]<-MCAR_30P_varOR
MCAR_30P_result[,3]<-MCAR_30P_meanOR-qnorm(0.975) * sqrt(MCAR_30P_varOR)/sqrt(n)
MCAR_30P_result[,4]<-MCAR_30P_meanOR+qnorm(0.975) * sqrt(MCAR_30P_varOR)/sqrt(n)
MCAR_30P_result = data.frame(MCAR_30P_result)
names(MCAR_30P_result)<- c("mean","var","low","high")
row.names(MCAR_30P_result)<-names(MCAR_30P_meanOR)
# Var plot #
# get rid of sex
MCAR_30P_result1 <-MCAR_30P_result[c(-2,-6,-10),]
MCAR_30P_result1$var <- round(MCAR_30P_result1$var,2)
MCAR_30P_result1$OR = factor(c("Active Control vs Placebo","140 mg Test Drug vs Placebo","210 mg Test Drug vs Placebo","Active Control vs Placebo","140 mg Test Drug vs Placebo","210 mg Test Drug vs Placebo","Active Control vs Placebo","140 mg Test Drug vs Placebo","210 mg Test Drug vs Placebo"))
MCAR_30P_result1$dataHandling = factor(c("completers","completers","completers","composite","composite","composite","LOCF","LOCF","LOCF"))
MCAR30PvarPlot <- ggplot(MCAR_30P_result1, aes(fill=OR, y=var, x=dataHandling)) + 
  geom_bar(position="dodge", stat="identity")+theme_bw()+
  theme(panel.grid.minor = element_blank(),panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"))+
  theme(plot.title = element_text(hjust = 0.5,size=13))+
  ylab("Var") +
  xlab("Strategies") +
  ggtitle("30% Missing, MCAR")+ylim(0,162)+
  geom_text(aes(x = dataHandling, y = var,label=var),size=3,vjust=-0.5,position = position_dodge(width=0.9))
MCAR30PvarPlot
#####################################################
############### MAR #################################
#####################################################
MAR_get_OR <- function(SEED){
  set.seed(SEED)
  # subjects with 30% missing visit4
  subjidMARvisit4high<- adpa %>%
    filter(AVISITN == 3 & PARAMN == 10 & PCHG < 10) %>%
    pull(SUBJID)
  subjidMARvisit4high_sample <- sample(subjidMARvisit4high,0.3*length(subjidMARvisit4high))
  allSubjid <- 1:1831
  subjidMARvisit4low <- setdiff(allSubjid,subjidMARvisit4high) 
  subjidMARvisit4low_sample <- sample(subjidMARvisit4low,0.05*length(subjidMARvisit4low))
  adpa_mar <- adpa %>% 
    filter(PARAMN == 10) %>%
    left_join(select(adsl,SUBJID,SEX),by="SUBJID")%>%
    group_by(SUBJID) %>%
    mutate(LAGPCHGCA1N = lag(PCHGCA1N,n=1,order_by = SUBJID)) %>%
    mutate(AVAL = ifelse(SUBJID %in% subjidMARvisit4high_sample & PARAMN == 10 & AVISITN %in% c(4,5,6),NA,AVAL)) %>%
    mutate(PCHG = ifelse(SUBJID %in% subjidMARvisit4high_sample & PARAMN == 10 & AVISITN %in% c(4,5,6),NA,PCHG)) %>%
    mutate(AVAL = ifelse(SUBJID %in% subjidMARvisit4low_sample & PARAMN == 10 & AVISITN %in% c(4,5,6),NA,AVAL)) %>%
    mutate(PCHG = ifelse(SUBJID %in% subjidMARvisit4low_sample & PARAMN == 10 & AVISITN %in% c(4,5,6),NA,PCHG)) %>%
    mutate(PCHGCA1N = ifelse(SUBJID %in% subjidMARvisit4low_sample & PARAMN == 10 & AVISITN %in% c(4,5,6),NA,PCHGCA1N)) %>%
    mutate(PCHGCA1N = ifelse(SUBJID %in% subjidMARvisit4high_sample & PARAMN == 10 & AVISITN %in% c(4,5,6),NA,PCHGCA1N))
  subjidMARvisit5high<- adpa_mar %>%
    filter(AVISITN == 4 & PARAMN == 10 & PCHG < 10) %>%
    pull(SUBJID)
  subjidMARvisit5high_sample <- sample(subjidMARvisit5high,0.3*length(subjidMARvisit5high))
  # subjects with 5% missing visit5
  subjidMARvisit5low <- setdiff(allSubjid,subjidMARvisit5high) 
  subjidMARvisit5low_sample <- sample(subjidMARvisit5low,0.05*length(subjidMARvisit5low))
  adpa_mar <- adpa_mar %>% 
    mutate(AVAL = ifelse(SUBJID %in% subjidMARvisit5high_sample & PARAMN == 10 & AVISITN %in% c(5,6),NA,AVAL)) %>%
    mutate(PCHG = ifelse(SUBJID %in% subjidMARvisit5high_sample & PARAMN == 10 & AVISITN %in% c(5,6),NA,PCHG)) %>%
    mutate(AVAL = ifelse(SUBJID %in% subjidMARvisit5low_sample & PARAMN == 10 & AVISITN %in% c(5,6),NA,AVAL)) %>%
    mutate(PCHG = ifelse(SUBJID %in% subjidMARvisit5low_sample & PARAMN == 10 & AVISITN %in% c(5,6),NA,PCHG)) %>%
    mutate(PCHGCA1N = ifelse(SUBJID %in% subjidMARvisit5low_sample & PARAMN == 10 & AVISITN %in% c(5,6),NA,PCHGCA1N)) %>%
    mutate(PCHGCA1N = ifelse(SUBJID %in% subjidMARvisit5high_sample & PARAMN == 10 & AVISITN %in% c(5,6),NA,PCHGCA1N)) 
  # subjects with 30% missing visit6
  subjidMARvisit6high<- adpa_mar %>%
    filter(AVISITN == 5 & PARAMN == 10 & PCHG < 10) %>%
    pull(SUBJID)
  subjidMARvisit6high_sample <- sample(subjidMARvisit6high,0.3*length(subjidMARvisit6high))
  subjidMARvisit6low <- setdiff(allSubjid,subjidMARvisit6high) 
  subjidMARvisit6low_sample <- sample(subjidMARvisit6low,0.05*length(subjidMARvisit6low))
  adpa_mar <- adpa_mar %>% 
    mutate(AVAL = ifelse(SUBJID %in% subjidMARvisit6high_sample & PARAMN == 10 & AVISITN == 6,NA,AVAL)) %>%
    mutate(PCHG = ifelse(SUBJID %in% subjidMARvisit6high_sample & PARAMN == 10 & AVISITN ==6,NA,PCHG)) %>%
    mutate(AVAL = ifelse(SUBJID %in% subjidMARvisit6low_sample & PARAMN == 10 & AVISITN == 6,NA,AVAL)) %>%
    mutate(PCHG = ifelse(SUBJID %in% subjidMARvisit6low_sample & PARAMN == 10 & AVISITN == 6,NA,PCHG)) %>%
    mutate(PCHGCA1N = ifelse(SUBJID %in% subjidMARvisit6low_sample & PARAMN == 10 & AVISITN == 6,NA,PCHGCA1N)) %>%
    mutate(PCHGCA1N = ifelse(SUBJID %in% subjidMARvisit6high_sample & PARAMN == 10 & AVISITN == 6,NA,PCHGCA1N)) %>%
    select(c(1:20,29,30))%>%
    filter(PARAMN == 10 & !(AVISITN == 2)) 
  adpa_mar_completer <- adpa_mar %>%
    filter(!is.na(PCHGCA1N) & AVISITN == 6)
  adpa_mar_imputed <- adpa_mar%>% 
    filter(AVISITN == 6) %>%
    mutate(PCHGCA1N=ifelse(is.na(PCHGCA1N),0,PCHGCA1N))
  adpa_mar_imputed_LOCF <- adpa_mar %>%
    mutate(PCHGCA1N=ifelse(is.na(PCHGCA1N) & AVISIT == 3,0,PCHGCA1N)) %>%
    mutate(PCHGCA1N = ifelse(is.na(PCHGCA1N) & !(AVISIT == 3),LAGPCHGCA1N,PCHGCA1N))%>%
    filter(AVISITN == 6)
  # completers
  model1 <- glm(PCHGCA1N~SEX+TRTP, family = binomial(),data = adpa_mar_completer)
  coef.model1 <- exp(coefficients(model1))
  ORcompletersSexM2F <- as.numeric(round(coef.model1[2],4)) # OR male to female
  ORcompletersAvsP <- as.numeric(round(1/coef.model1[3],4)) # OR active control vs placebo
  ORcompleters140vsP <- as.numeric(round(coef.model1[4]/coef.model1[3],4)) # OR 140mg vs placebo
  ORcompleters210vsP <- as.numeric(round(coef.model1[5]/coef.model1[3],4)) # OR 140mg vs placebo
  ORcompleters <- list(SexM2F = ORcompletersSexM2F,AvsP = ORcompletersAvsP,T140vsP = ORcompleters140vsP, T210vsP =ORcompleters210vsP)
  # imputed
  model2 <- glm(PCHGCA1N~SEX+TRTP, family = binomial(),data = adpa_mar_imputed)
  coef.model2 <- exp(coefficients(model2))
  ORimputedSexM2F <- as.numeric(round(coef.model2[2],4)) # OR male to female
  ORimputedAvsP <- as.numeric(round(1/coef.model2[3],4)) # OR active control vs placebo
  ORimputed140vsP <- as.numeric(round(coef.model2[4]/coef.model2[3],4)) # OR 140mg vs placebo
  ORimputed210vsP <- as.numeric(round(coef.model2[5]/coef.model2[3],4)) # OR 140mg vs placebo
  ORimputed <- list(SexM2F = ORimputedSexM2F,AvsP = ORimputedAvsP,T140vsP = ORimputed140vsP, T210vsP =ORimputed210vsP)
  # LOCF
  model3 <- glm(PCHGCA1N~SEX+TRTP, family = binomial(),data = adpa_mar_imputed_LOCF)
  coef.model3 <- exp(coefficients(model3))
  ORimputedSexM2F <- as.numeric(round(coef.model3[2],4)) # OR male to female
  ORimputedAvsP <- as.numeric(round(1/coef.model3[3],4)) # OR active control vs placebo
  ORimputed140vsP <- as.numeric(round(coef.model3[4]/coef.model3[3],4)) # OR 140mg vs placebo
  ORimputed210vsP <- as.numeric(round(coef.model3[5]/coef.model3[3],4)) # OR 140mg vs placebo
  ORlocf <- list(SexM2F = ORimputedSexM2F,AvsP = ORimputedAvsP,T140vsP = ORimputed140vsP, T210vsP =ORimputed210vsP)
  ##########
  # return  
  c(list(ORcompleters = ORcompleters,ORimputed=ORimputed,ORlocf = ORlocf))
}
######### Run MAR_get_OR with 1000 seeds ######
MAR_result<-lapply(1:1000,FUN = MAR_get_OR)
##### get mean, var, CI for MAR_result ######
tmp.mar<-unlist(MAR_result)
MAR_meanOR<- c(by(tmp.mar, names(tmp.mar), mean))
MAR_varOR<- c(by(tmp.mar, names(tmp.mar), var))
MAR_result <- matrix(,nrow=12,ncol = 4)
MAR_result[,1]<-MAR_meanOR
MAR_result[,2]<-MAR_varOR
MAR_result[,3]<-MAR_meanOR-qnorm(0.975) * sqrt(MAR_result[,2]<-MAR_varOR)/sqrt(n)
MAR_result[,4]<-MAR_meanOR+qnorm(0.975) * sqrt(MAR_result[,2]<-MAR_varOR)/sqrt(n)
MAR_result = data.frame(MAR_result)
names(MAR_result)<- c("mean","var","low","high")
row.names(MAR_result)<-names(MAR_meanOR)
###### complete ######
completeData<-adpa %>%
  filter(PARAMN == 10 & !(AVISITN == 2)) %>%
  left_join(select(adsl,SUBJID,SEX),by="SUBJID") %>%
  select(c(1:20,29)) %>%
  filter(AVISITN == 6)
model3 <- glm(PCHGCA1N~SEX+TRTP, family = binomial(),data = completeData)
coef.model3 <- exp(coefficients(model3))
ORcompleteSexM2F <- as.numeric(round(coef.model3[2],4)) # OR male to female
ORcompleteAvsP <- as.numeric(round(1/coef.model3[3],4)) # OR active control vs placebo
ORcomplete140vsP <- as.numeric(round(coef.model3[4]/coef.model3[3],4)) # OR 140mg vs placebo
ORcomplete210vsP <- as.numeric(round(coef.model3[5]/coef.model3[3],4)) # OR 140mg vs placebo
ORcomplete <- list(SexM2F = ORcompleteSexM2F,AvsP = ORcompleteAvsP,T140vsP = ORcomplete140vsP, T210vsP =ORcomplete210vsP)
threeWayTableComplete<-ftable(table(completeData$TRTPN,completeData$SEX,completeData$PCHGCA1N))
############## MCAR plot 10p ###########################
boxLabels=factor(c(1,2,3,4,1,2,3,4,1,2,3,4,1,2,3,4))
boxLabels=c('TRTP Active control \n vs Placebo','TRTP Test drug 140mg \n vs Placebo ','TRTP Test drug 210mg \n vs Placebo ','Sex(M vs F)')
# plot 1: MCAR
df.10p = data.frame(
  yAxis = boxLabels,
  boxOdds = c(MCAR_10P_result$mean[1],MCAR_10P_result$mean[3],MCAR_10P_result$mean[4],MCAR_10P_result$mean[2],
              MCAR_10P_result$mean[5],MCAR_10P_result$mean[7],MCAR_10P_result$mean[8],MCAR_10P_result$mean[6],
              MCAR_10P_result$mean[9],MCAR_10P_result$mean[11],MCAR_10P_result$mean[12],MCAR_10P_result$mean[10],
              ORcomplete$AvsP,ORcomplete$T140vsP,ORcomplete$T210vsP,ORcomplete$SexM2F),
  boxCILow = c(MCAR_10P_result$low[1],MCAR_10P_result$low[3],MCAR_10P_result$low[4],MCAR_10P_result$low[2],
               MCAR_10P_result$low[5],MCAR_10P_result$low[7],MCAR_10P_result$low[8],MCAR_10P_result$low[6],
               MCAR_10P_result$low[9],MCAR_10P_result$low[11],MCAR_10P_result$low[12],MCAR_10P_result$low[10],
               ORcomplete$AvsP,ORcomplete$T140vsP,ORcomplete$T210vsP,ORcomplete$SexM2F),
  boxCIHigh = c(MCAR_10P_result$high[1],MCAR_10P_result$high[3],MCAR_10P_result$high[4],MCAR_10P_result$high[2],
                MCAR_10P_result$high[5],MCAR_10P_result$high[7],MCAR_10P_result$high[8],MCAR_10P_result$high[6],
                MCAR_10P_result$high[9],MCAR_10P_result$high[11],MCAR_10P_result$high[12],MCAR_10P_result$high[10],
                ORcomplete$AvsP,ORcomplete$T140vsP,ORcomplete$T210vsP,ORcomplete$SexM2F),
  group=factor(c('Completers Data','Completers Data','Completers Data','Completers Data',
                 'Imputed Data (Composite)','Imputed Data (Composite)','Imputed Data (Composite)','Imputed Data (Composite)',
                 'Imputed Data (LOCF)','Imputed Data (LOCF)','Imputed Data (LOCF)','Imputed Data (LOCF)',
                 'Complete Data','Complete Data','Complete Data','Complete Data')),
  VAL=c(round(MCAR_10P_result$mean[1],3),round(MCAR_10P_result$mean[3],3),round(MCAR_10P_result$mean[4],3),round(MCAR_10P_result$mean[2],3),
        round(MCAR_10P_result$mean[5],3),round(MCAR_10P_result$mean[7],3),round(MCAR_10P_result$mean[8],3),round(MCAR_10P_result$mean[6],3),
        round(MCAR_10P_result$mean[9],3),round(MCAR_10P_result$mean[11],3),round(MCAR_10P_result$mean[12],3),round(MCAR_10P_result$mean[10],3),
        round(ORcomplete$AvsP,3),round(ORcomplete$T140vsP,3),round(ORcomplete$T210vsP,3),round(ORcomplete$SexM2F,3))
)

p = ggplot(df.10p, aes(x = boxOdds, y = yAxis, colour = group))
p.10p <-p + geom_errorbarh(aes(y = yAxis, xmax = boxCIHigh, xmin = boxCILow), size = .5, height = .2, position = position_dodgev(height=0.9)) +
  geom_point(aes(x = boxOdds, y = yAxis),size = 2.5,  position = position_dodgev(height=0.9))+
  geom_text(aes(x = boxOdds, y = yAxis,label=VAL),size=3,vjust=-0.5,position = position_dodgev(height=0.9))+
  theme_bw()+
  theme(panel.grid.minor = element_blank(),panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"))+
  theme(plot.title = element_text(hjust = 0.5,size=13))+
  scale_x_continuous(breaks = c(1,seq(0,100,20)))+
  coord_trans(x = 'log10')+
  ylab("") +
  xlab("Odds ratio (log scale)") +
  ggtitle("10% missing, MCAR")+ 
  scale_color_discrete(
    guide = guide_legend(override.aes = list(color = "white"))
  ) + theme(
    legend.text = element_text(color = "white"),
    legend.title = element_text(color = "white"),
    legend.key = element_rect(fill = "white")
  )
p.10p
##################### MCAR 20p ######################
boxLabels=factor(c(1,2,3,4,1,2,3,4,1,2,3,4,1,2,3,4))
boxLabels=c('TRTP Active control \n vs Placebo','TRTP Test drug 140mg \n vs Placebo ','TRTP Test drug 210mg \n vs Placebo ','Sex(M vs F)')
# plot 1: MCAR
df.20p = data.frame(
  yAxis = boxLabels,
  boxOdds = c(MCAR_20P_result$mean[1],MCAR_20P_result$mean[3],MCAR_20P_result$mean[4],MCAR_20P_result$mean[2],
              MCAR_20P_result$mean[5],MCAR_20P_result$mean[7],MCAR_20P_result$mean[8],MCAR_20P_result$mean[6],
              MCAR_20P_result$mean[9],MCAR_20P_result$mean[11],MCAR_20P_result$mean[12],MCAR_20P_result$mean[10],
              ORcomplete$AvsP,ORcomplete$T140vsP,ORcomplete$T210vsP,ORcomplete$SexM2F),
  boxCILow = c(MCAR_20P_result$low[1],MCAR_20P_result$low[3],MCAR_20P_result$low[4],MCAR_20P_result$low[2],
               MCAR_20P_result$low[5],MCAR_20P_result$low[7],MCAR_20P_result$low[8],MCAR_20P_result$low[6],
               MCAR_20P_result$low[9],MCAR_20P_result$low[11],MCAR_20P_result$low[12],MCAR_20P_result$low[10],
               ORcomplete$AvsP,ORcomplete$T140vsP,ORcomplete$T210vsP,ORcomplete$SexM2F),
  boxCIHigh = c(MCAR_20P_result$high[1],MCAR_20P_result$high[3],MCAR_20P_result$high[4],MCAR_20P_result$high[2],
                MCAR_20P_result$high[5],MCAR_20P_result$high[7],MCAR_20P_result$high[8],MCAR_20P_result$high[6],
                MCAR_20P_result$high[9],MCAR_20P_result$high[11],MCAR_20P_result$high[12],MCAR_20P_result$high[10],
                ORcomplete$AvsP,ORcomplete$T140vsP,ORcomplete$T210vsP,ORcomplete$SexM2F),
  group=factor(c('Completers Data','Completers Data','Completers Data','Completers Data',
                 'Imputed Data (Composite)','Imputed Data (Composite)','Imputed Data (Composite)','Imputed Data (Composite)',
                 'Imputed Data (LOCF)','Imputed Data (LOCF)','Imputed Data (LOCF)','Imputed Data (LOCF)',
                 'Complete Data','Complete Data','Complete Data','Complete Data')),
  VAL=c(round(MCAR_20P_result$mean[1],3),round(MCAR_20P_result$mean[3],3),round(MCAR_20P_result$mean[4],3),round(MCAR_20P_result$mean[2],3),
        round(MCAR_20P_result$mean[5],3),round(MCAR_20P_result$mean[7],3),round(MCAR_20P_result$mean[8],3),round(MCAR_20P_result$mean[6],3),
        round(MCAR_20P_result$mean[9],3),round(MCAR_20P_result$mean[11],3),round(MCAR_20P_result$mean[12],3),round(MCAR_20P_result$mean[10],3),
        round(ORcomplete$AvsP,3),round(ORcomplete$T140vsP,3),round(ORcomplete$T210vsP,3),round(ORcomplete$SexM2F,3))
)

p = ggplot(df.20p, aes(x = boxOdds, y = yAxis, colour = group))
p.20p <-p + geom_errorbarh(aes(y = yAxis, xmax = boxCIHigh, xmin = boxCILow), size = .5, height = .2, position = position_dodgev(height=0.9)) +
  geom_point(aes(x = boxOdds, y = yAxis),size = 2.5,  position = position_dodgev(height=0.9))+
  geom_text(aes(x = boxOdds, y = yAxis,label=VAL),size=3,vjust=-0.5,position = position_dodgev(height=0.9))+
  theme_bw()+
  theme(panel.grid.minor = element_blank(),panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"))+
  theme(plot.title = element_text(hjust = 0.5,size=13))+
  scale_x_continuous(breaks = c(1,seq(0,100,20)))+
  coord_trans(x = 'log10')+
  ylab("") +
  xlab("Odds ratio (log scale)") +
  ggtitle("20% missing, MCAR") + 
  scale_color_discrete(
    guide = guide_legend(override.aes = list(color = "white"))
  ) + theme(
    legend.text = element_text(color = "white"),
    legend.title = element_text(color = "white"),
    legend.key = element_rect(fill = "white")
  )
p.20p
################## MCAR 30p ###################
#################################################
boxLabels=factor(c(1,2,3,4,1,2,3,4,1,2,3,4,1,2,3,4))
boxLabels=c('TRTP Active control \n vs Placebo','TRTP Test drug 140mg \n vs Placebo ','TRTP Test drug 210mg \n vs Placebo ','Sex(M vs F)')
# plot 1: MCAR
df.30p = data.frame(
  yAxis = boxLabels,
  boxOdds = c(MCAR_30P_result$mean[1],MCAR_30P_result$mean[3],MCAR_30P_result$mean[4],MCAR_30P_result$mean[2],
              MCAR_30P_result$mean[5],MCAR_30P_result$mean[7],MCAR_30P_result$mean[8],MCAR_30P_result$mean[6],
              MCAR_30P_result$mean[9],MCAR_30P_result$mean[11],MCAR_30P_result$mean[12],MCAR_30P_result$mean[10],
              ORcomplete$AvsP,ORcomplete$T140vsP,ORcomplete$T210vsP,ORcomplete$SexM2F),
  boxCILow = c(MCAR_30P_result$low[1],MCAR_30P_result$low[3],MCAR_30P_result$low[4],MCAR_30P_result$low[2],
               MCAR_30P_result$low[5],MCAR_30P_result$low[7],MCAR_30P_result$low[8],MCAR_30P_result$low[6],
               MCAR_30P_result$low[9],MCAR_30P_result$low[11],MCAR_30P_result$low[12],MCAR_30P_result$low[10],
               ORcomplete$AvsP,ORcomplete$T140vsP,ORcomplete$T210vsP,ORcomplete$SexM2F),
  boxCIHigh = c(MCAR_30P_result$high[1],MCAR_30P_result$high[3],MCAR_30P_result$high[4],MCAR_30P_result$high[2],
                MCAR_30P_result$high[5],MCAR_30P_result$high[7],MCAR_30P_result$high[8],MCAR_30P_result$high[6],
                MCAR_30P_result$high[9],MCAR_30P_result$high[11],MCAR_30P_result$high[12],MCAR_30P_result$high[10],
                ORcomplete$AvsP,ORcomplete$T140vsP,ORcomplete$T210vsP,ORcomplete$SexM2F),
  group=factor(c('Completers Data','Completers Data','Completers Data','Completers Data',
                 'Imputed Data (Composite)','Imputed Data (Composite)','Imputed Data (Composite)','Imputed Data (Composite)',
                 'Imputed Data (LOCF)','Imputed Data (LOCF)','Imputed Data (LOCF)','Imputed Data (LOCF)',
                 'Complete Data','Complete Data','Complete Data','Complete Data')),
  VAL=c(round(MCAR_30P_result$mean[1],3),round(MCAR_30P_result$mean[3],3),round(MCAR_30P_result$mean[4],3),round(MCAR_30P_result$mean[2],3),
        round(MCAR_30P_result$mean[5],3),round(MCAR_30P_result$mean[7],3),round(MCAR_30P_result$mean[8],3),round(MCAR_30P_result$mean[6],3),
        round(MCAR_30P_result$mean[9],3),round(MCAR_30P_result$mean[11],3),round(MCAR_30P_result$mean[12],3),round(MCAR_30P_result$mean[10],3),
        round(ORcomplete$AvsP,3),round(ORcomplete$T140vsP,3),round(ORcomplete$T210vsP,3),round(ORcomplete$SexM2F,3))
)

p = ggplot(df.30p, aes(x = boxOdds, y = yAxis, colour = group))
p.30p <-p + geom_errorbarh(aes(y = yAxis, xmax = boxCIHigh, xmin = boxCILow), size = .5, height = .2, position = position_dodgev(height=0.9)) +
  geom_point(aes(x = boxOdds, y = yAxis),size = 2.5,  position = position_dodgev(height=0.9))+
  geom_text(aes(x = boxOdds, y = yAxis,label=VAL),size=3,vjust=-0.5,position = position_dodgev(height=0.9))+
  theme_bw()+
  theme(panel.grid.minor = element_blank(),panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"))+
  theme(plot.title = element_text(hjust = 0.5,size=13))+
  scale_x_continuous(breaks = c(1,seq(0,100,20)))+
  coord_trans(x = 'log10')+
  ylab("") +
  xlab("Odds ratio (log scale)") +
  ggtitle("30% missing, MCAR")
p.30p
################# MAR plot ##############
boxLabels=factor(c(1,2,3,4,1,2,3,4,1,2,3,4,1,2,3,4))
boxLabels=c('TRTP Active control \n vs Placebo','TRTP Test drug 140mg \n vs Placebo ','TRTP Test drug 210mg \n vs Placebo ','Sex(M vs F)')
# plot 1: MCAR
df.mar = data.frame(
  yAxis = boxLabels,
  boxOdds = c(MAR_result$mean[1],MAR_result$mean[3],MAR_result$mean[4],MAR_result$mean[2],
              MAR_result$mean[5],MAR_result$mean[7],MAR_result$mean[8],MAR_result$mean[6],
              MAR_result$mean[9],MAR_result$mean[11],MAR_result$mean[12],MAR_result$mean[10],
              ORcomplete$AvsP,ORcomplete$T140vsP,ORcomplete$T210vsP,ORcomplete$SexM2F),
  boxCILow = c(MAR_result$low[1],MAR_result$low[3],MAR_result$low[4],MAR_result$low[2],
               MAR_result$low[5],MAR_result$low[7],MAR_result$low[8],MAR_result$low[6],
               MAR_result$low[9],MAR_result$low[11],MAR_result$low[12],MAR_result$low[10],
               ORcomplete$AvsP,ORcomplete$T140vsP,ORcomplete$T210vsP,ORcomplete$SexM2F),
  boxCIHigh = c(MAR_result$high[1],MAR_result$high[3],MAR_result$high[4],MAR_result$high[2],
                MAR_result$high[5],MAR_result$high[7],MAR_result$high[8],MAR_result$high[6],
                MAR_result$high[9],MAR_result$high[11],MAR_result$high[12],MAR_result$high[10],
                ORcomplete$AvsP,ORcomplete$T140vsP,ORcomplete$T210vsP,ORcomplete$SexM2F),
  group=factor(c('Completers Data','Completers Data','Completers Data','Completers Data',
                 'Imputed Data (Composite)','Imputed Data (Composite)','Imputed Data (Composite)','Imputed Data (Composite)',
                 'Imputed Data (LOCF)','Imputed Data (LOCF)','Imputed Data (LOCF)','Imputed Data (LOCF)',
                 'Complete Data','Complete Data','Complete Data','Complete Data')),
  VAL=c(round(MAR_result$mean[1],3),round(MAR_result$mean[3],3),round(MAR_result$mean[4],3),round(MAR_result$mean[2],3),
        round(MAR_result$mean[5],3),round(MAR_result$mean[7],3),round(MAR_result$mean[8],3),round(MAR_result$mean[6],3),
        round(MAR_result$mean[9],3),round(MAR_result$mean[11],3),round(MAR_result$mean[12],3),round(MAR_result$mean[10],3),
        round(ORcomplete$AvsP,3),round(ORcomplete$T140vsP,3),round(ORcomplete$T210vsP,3),round(ORcomplete$SexM2F,3))
)

p = ggplot(df.mar, aes(x = boxOdds, y = yAxis, colour = group))
p.mar <-p + geom_errorbarh(aes(y = yAxis, xmax = boxCIHigh, xmin = boxCILow), size = .5, height = .2, position = position_dodgev(height=0.9)) +
  geom_point(aes(x = boxOdds, y = yAxis),size = 2.5,  position = position_dodgev(height=0.9))+
  geom_text(aes(x = boxOdds, y = yAxis,label=VAL),size=3,vjust=-0.5,position = position_dodgev(height=0.9))+
  theme_bw()+
  theme(panel.grid.minor = element_blank(),panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"))+
  theme(plot.title = element_text(hjust = 0.5,size=13))+
  scale_x_continuous(breaks = c(1,seq(0,100,20)))+
  coord_trans(x = 'log10')+
  ylab("") +
  xlab("Odds ratio (log scale)") +
  ggtitle("MAR, Mean of OR")
p.mar
# Var plot #
MAR_result1 <- MAR_result[c(-2,-6,-10),]
MAR_result1$var <- round(MAR_result1$var,2)
MAR_result1$OR = factor(c("Active Control vs Placebo","140 mg Test Drug vs Placebo","210 mg Test Drug vs Placebo","Active Control vs Placebo","140 mg Test Drug vs Placebo","210 mg Test Drug vs Placebo","Active Control vs Placebo","140 mg Test Drug vs Placebo","210 mg Test Drug vs Placebo"))
MAR_result1$dataHandling = factor(c("completers","completers","completers","composite","composite","composite","LOCF","LOCF","LOCF"))
MARvarPlot <- ggplot(MAR_result1, aes(fill=OR, y=var, x=dataHandling)) + 
  geom_bar(position="dodge", stat="identity")+theme_bw()+
  theme(panel.grid.minor = element_blank(),panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"))+
  theme(plot.title = element_text(hjust = 0.5,size=13))+
  ylab("Var") +
  xlab("Strategies") +
  ggtitle("MAR, Variance of OR")+
  geom_text(aes(x = dataHandling, y = var,label=var),size=3,vjust=-0.5,position = position_dodge(width=0.9))
MARvarPlot 
# Combine Plots
theme_set(theme_pubr())
figure <- ggarrange(
  ggarrange(p.10p, p.20p),
  p.30p,
  nrow=2
)
figure

figure <- ggarrange(
  p.10p,p.20p,p.30p,ncol=3
)
