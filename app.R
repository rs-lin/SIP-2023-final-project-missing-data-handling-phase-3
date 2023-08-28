library(shiny)
library(dplyr)
library(ggplot2)
library(epitools)
library(descr)
library(ggstance)

adpa <- read.csv("data_phase 3/ADPA.csv")
adsl <- read.csv("data_phase 3/ADSL.csv")
############### Complete Data ##################
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
completeData.Result.df <- data.frame(mean = c(ORcompleteAvsP,ORcompleteSexM2F,ORcomplete140vsP,ORcomplete210vsP))
row.names(completeData.Result.df) <- c('ORcomplete.AvsP','ORcomplete.SexM2F','ORcomplete.T140vsP','ORcomplete.T210vsP')
### plot complete ###
boxLabels=factor(c(1,2,3,4))
boxLabels=c('TRTP Active control \n vs Placebo','TRTP Test drug 140mg \n vs Placebo ','TRTP Test drug 210mg \n vs Placebo ','Sex(M vs F)')
# plot 1: MCAR
df.complete = data.frame(
  yAxis = boxLabels,
  boxOdds = c(completeData.Result.df$mean[1],completeData.Result.df$mean[3],completeData.Result.df$mean[4],completeData.Result.df$mean[2]),
  VAL=c(round(completeData.Result.df$mean[1],3),round(completeData.Result.df$mean[3],3),round(completeData.Result.df$mean[4],3),round(completeData.Result.df$mean[2],3))
)
p.complete = ggplot(df.complete, aes(x = boxOdds, y = yAxis))
p.complete <-p.complete  +
  geom_point(aes(x = boxOdds, y = yAxis),size = 2.5, color = "salmon", position = position_dodgev(height=0.9))+
  geom_text(aes(x = boxOdds, y = yAxis,label=VAL),size=3,vjust=-0.5,position = position_dodgev(height=0.9), color = "salmon")+
  theme_bw()+
  theme(panel.grid.minor = element_blank(),panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"))+
  theme(plot.title = element_text(hjust = 0.5,size=13))+
  scale_x_continuous(breaks = c(1,seq(0,100,20)))+
  coord_trans(x = 'log10')+
  ylab("") +
  xlab("Odds ratio (log scale)") +
  ggtitle("Odds Ratio: Complete Data")

ui<- shinyUI(pageWithSidebar(
  headerPanel("MCAR and MAR Example"),
  sidebarPanel(
    p("Please Select One Missing Mechanism"),
    wellPanel(
      checkboxInput("MCAR", "Proceed to MCAR Analysis", FALSE),
      conditionalPanel(
      condition="input.MCAR==true && input.MAR==false",
      selectInput(inputId = "variable1",label = "Select the Level of Missing for MCAR:", 
                  choices = c("30%" = "p1.MCAR",
                              "20%" = "p2.MCAR", 
                              "10%" = "p3.MCAR"),
                  selected = "30%"
      )
    )
  ),
    wellPanel(
      checkboxInput("MAR", "Proceed to MAR Analysis", FALSE),
      conditionalPanel(
        condition="input.MAR==true && input.MCAR==false",
        selectInput(inputId = "variable2", 
                    label = "Starting from visit 4, if the PASI score for visit 3 shows less than 10% improvement from baseline, then the probability of missing at visit 4 is",
                    choices = c("30%" = "p1.MAR.High",
                                "25%" = "p2.MAR.High", 
                                "20%" = "p3.MAR.High"),
                    selected = "30%"
        ),
        selectInput(inputId = "variable3", 
                    label = "Otherwise, the probability of missing at visit 4 is",
                    choices = c("10%" = "p1.MAR.Low",
                                "7.5%" = "p2.MAR.Low", 
                                "5%" = "p3.MAR.Low"),
                    selected = "10%"
        )
      )
    )
  ),
  
  mainPanel(
    h5("Output"),
    tabsetPanel(id ="analysisTabs",
                tabPanel(title = "True Result", value="panel_complete_data",
                         h4("True Result"),                       
                         plotOutput(outputId = "target")),
                tabPanel(title = "Completers Only", value="panel_completers",
                         h4("Completers Only"),                       
                         verbatimTextOutput("Completers Only"),
                         # p("Your result is", textOutput(outputId = "result", inline=T), "."),
                         plotOutput(outputId = "completers")),
                tabPanel(title = "Imputed with Composite Strategy", value="panel_composite",
                         h4("Imputed with Composite Strategy"),                       
                         verbatimTextOutput("Imputed with Composite Strategy"),
                         plotOutput(outputId = "imputed")),
                tabPanel(title = "Imputed with LOCF", value="panel_LOCF",
                         h4("Imputed with LOCF"),                       
                         verbatimTextOutput("Imputed with LOCF"),
                         plotOutput(outputId = "LOCF"))
    )
  )  
))
# end UI
server <- shinyServer(function(input, output, session) {
  set.seed(rnorm(1))
  ##############True Result Tab###############
  output$target <- renderPlot({p.complete})
  ################Extract Percentages###############
  MCARpercentage <- reactive({
    require(input$MCAR)
    require(input$MAR)
    if((input$MCAR == TRUE) & (input$MAR == FALSE)){
      require(input$variable1)
      if (input$variable1 == "p1.MCAR"){
        percent = 0.3
      }
      else if (input$variable1 == "p2.MCAR"){
        percent = 0.2
      }
      else{
        percent = 0.1
      }
    }
    else{
      percent <- NULL
    }
    percent
  })
  MARpercentage1 <- reactive({
    require(input$MCAR)
    require(input$MAR)
    if((input$MCAR == FALSE) & (input$MAR == TRUE)){
      require(input$variable2)
      require(input$variable3)
      if (input$variable2 == "p1.MAR.High"){
        percent1 = 0.3
      }
      else if (input$variable2 == "p2.MAR.High"){
        percent1 = 0.25
      }
      else{
        percent1 =0.2
      }
    }
    else{
      percent1<- NULL
    }
    percent1
  })
  MARpercent2 <- reactive({
    require(input$MCAR)
    require(input$MAR)
    if((input$MCAR == FALSE) & (input$MAR == TRUE)){
      require(input$variable2)
      require(input$variable3)
      if (input$variable3 == "p1.MAR.Low"){
        percent2 = 0.1
      }
      else if (input$variable2 == "p2.MAR.Low"){
        percent2 = 0.075
      }
      else{
        percent2 =0.05
      }
    }
    else{
      percent2 <- NULL
    }
    percent2
  })
  
  #output$result <- MARpercent2
  output$completers <- renderPlot({
    if((input$MCAR == TRUE) & (input$MAR == FALSE)){
      subjid_sample<-sample(1:1831,1831*MCARpercentage())
      completersData <- adpa %>%
        filter(!(SUBJID %in% subjid_sample) & PARAMN == 10 & !(AVISITN == 2)) %>%
        left_join(select(adsl,SUBJID,SEX),by="SUBJID") %>%
        select(c(1:20,29)) %>%
        filter(AVISITN == 6)
      model3<-glm(PCHGCA1N~SEX+TRTP, family = binomial(), data = completersData)
      coef.model3 <- exp(coefficients(model3))
      ORcompletersSexM2F <- as.numeric(round(coef.model3[2],4)) # OR male to female
      ORcompletersAvsP <- as.numeric(round(1/coef.model3[3],4)) # OR active control vs placebo
      ORcompleters140vsP <- as.numeric(round(coef.model3[4]/coef.model3[3],4)) # OR 140mg vs placebo
      ORcompleters210vsP <- as.numeric(round(coef.model3[5]/coef.model3[3],4)) # OR 140mg vs placebo
      ORcompleters <- list(SexM2F = ORcompletersSexM2F,AvsP = ORcompletersAvsP,T140vsP = ORcompleters140vsP, T210vsP =ORcompleters210vsP)
      # threeWayTable10PCompleters<-ftable(table(completers10PData$TRTPN,completers10PData$SEX,completers10PData$PCHGCA1N))
      completersData.Result.df <- data.frame(mean = c(ORcompletersAvsP,ORcompletersSexM2F,ORcompleters140vsP,ORcompleters210vsP))
      boxLabels=factor(c(1,2,3,4,1,2,3,4))
      boxLabels=c('TRTP Active control \n vs Placebo','TRTP Test drug 140mg \n vs Placebo ','TRTP Test drug 210mg \n vs Placebo ','Sex(M vs F)')
      # plot 1: MCAR completers
      df.completers = data.frame(
        yAxis = boxLabels,
        boxOdds = c(completersData.Result.df$mean[1],completersData.Result.df$mean[3],completersData.Result.df$mean[4],completersData.Result.df$mean[2],
                    completeData.Result.df$mean[1],completeData.Result.df$mean[3],completeData.Result.df$mean[4],completeData.Result.df$mean[2]),
        group=factor(c('Completers Data','Completers Data','Completers Data','Completers Data',
                       'Complete Data','Complete Data','Complete Data','Complete Data')),
        VAL=c(round(completersData.Result.df$mean[1],3),round(completersData.Result.df$mean[3],3),round(completersData.Result.df$mean[4],3),round(completersData.Result.df$mean[2],3),
              round(completeData.Result.df$mean[1],3),round(completeData.Result.df$mean[3],3),round(completeData.Result.df$mean[4],3),round(completeData.Result.df$mean[2],3))
      )
      p.completers = ggplot(df.completers, aes(x = boxOdds, y = yAxis, colour = group))
      p.completers.mcar <-p.completers  +
        geom_point(aes(x = boxOdds, y = yAxis),size = 2.5,  position = position_dodgev(height=0.9))+
        geom_text(aes(x = boxOdds, y = yAxis,label=VAL),size=3,vjust=-0.5,position = position_dodgev(height=0.9))+
        theme_bw()+
        theme(panel.grid.minor = element_blank(),panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"))+
        theme(plot.title = element_text(hjust = 0.5,size=13))+
        scale_x_continuous(breaks = c(1,seq(0,100,20)))+
        coord_trans(x = 'log10')+
        ylab("") +
        xlab("Odds ratio (log scale)") +
        ggtitle("Odds Ratio: Completers Data, MCAR")
      p.completers.mcar
    }
    else if((input$MCAR == FALSE) & (input$MAR == TRUE)){
      subjidMARvisit4high<- adpa %>%
        filter(AVISITN == 3 & PARAMN == 10 & PCHG < 10) %>%
        pull(SUBJID)
      subjidMARvisit4high_sample <- sample(subjidMARvisit4high,MARpercentage1()*length(subjidMARvisit4high))
      # subjects with 5% missing visit4
      allSubjid <- 1:1831
      subjidMARvisit4low <- setdiff(allSubjid,subjidMARvisit4high) 
      subjidMARvisit4low_sample <- sample(subjidMARvisit4low,MARpercent2()*length(subjidMARvisit4low))
      # set the selected subjects' PASI score to NA on visit 4 and on subsequent visits
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
      # subjects with 30% missing visit5
      subjidMARvisit5high<- adpa_mar %>%
        filter(AVISITN == 4 & PARAMN == 10 & PCHG < 10) %>%
        pull(SUBJID)
      subjidMARvisit5high_sample <- sample(subjidMARvisit5high,MARpercentage1()*length(subjidMARvisit5high))
      # subjects with 5% missing visit5
      subjidMARvisit5low <- setdiff(allSubjid,subjidMARvisit5high) 
      subjidMARvisit5low_sample <- sample(subjidMARvisit5low,MARpercent2()*length(subjidMARvisit5low))
      # set the selected subjects' PASI score to NA on visit 5 and 6
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
      subjidMARvisit6high_sample <- sample(subjidMARvisit6high,MARpercentage1()*length(subjidMARvisit6high))
      # subjects with 5% missing visit5
      subjidMARvisit6low <- setdiff(allSubjid,subjidMARvisit6high) 
      subjidMARvisit6low_sample <- sample(subjidMARvisit6low,MARpercent2()*length(subjidMARvisit6low))
      adpa_mar <- adpa_mar %>% 
        mutate(AVAL = ifelse(SUBJID %in% subjidMARvisit6high_sample & PARAMN == 10 & AVISITN == 6,NA,AVAL)) %>%
        mutate(PCHG = ifelse(SUBJID %in% subjidMARvisit6high_sample & PARAMN == 10 & AVISITN ==6,NA,PCHG)) %>%
        mutate(AVAL = ifelse(SUBJID %in% subjidMARvisit6low_sample & PARAMN == 10 & AVISITN == 6,NA,AVAL)) %>%
        mutate(PCHG = ifelse(SUBJID %in% subjidMARvisit6low_sample & PARAMN == 10 & AVISITN == 6,NA,PCHG)) %>%
        mutate(PCHGCA1N = ifelse(SUBJID %in% subjidMARvisit6low_sample & PARAMN == 10 & AVISITN == 6,NA,PCHGCA1N)) %>%
        mutate(PCHGCA1N = ifelse(SUBJID %in% subjidMARvisit6high_sample & PARAMN == 10 & AVISITN == 6,NA,PCHGCA1N)) %>%
        select(c(1:20,29,30))%>%
        filter(PARAMN == 10 & AVISITN == 6) 
      # MAR completers
      adpa_mar_completer <-  adpa_mar %>% filter(!is.na(PCHGCA1N))
      model3<-glm(PCHGCA1N~SEX+TRTP, family = binomial(), data = adpa_mar_completer)
      coef.model3 <- exp(coefficients(model3))
      ORcompletersSexM2F <- as.numeric(round(coef.model3[2],4)) # OR male to female
      ORcompletersAvsP <- as.numeric(round(1/coef.model3[3],4)) # OR active control vs placebo
      ORcompleters140vsP <- as.numeric(round(coef.model3[4]/coef.model3[3],4)) # OR 140mg vs placebo
      ORcompleters210vsP <- as.numeric(round(coef.model3[5]/coef.model3[3],4)) # OR 140mg vs placebo
      #ORmarC <- list(SexM2F = ORcompletersSexM2F,AvsP = ORcompletersAvsP,T140vsP = ORcompleters140vsP, T210vsP =ORcompleters210vsP)
      MARcompleter.Result.df <- data.frame(mean = c(ORcompletersAvsP,ORcompletersSexM2F,ORcompleters140vsP,ORcompleters210vsP))
      boxLabels=factor(c(1,2,3,4,1,2,3,4))
      boxLabels=c('TRTP Active control \n vs Placebo','TRTP Test drug 140mg \n vs Placebo ','TRTP Test drug 210mg \n vs Placebo ','Sex(M vs F)')
      # plot 3: LOCF MCAR
      df.completer.MAR = data.frame(
        yAxis = boxLabels,
        boxOdds = c(MARcompleter.Result.df$mean[1],MARcompleter.Result.df$mean[3],MARcompleter.Result.df$mean[4],MARcompleter.Result.df$mean[2],
                    completeData.Result.df$mean[1],completeData.Result.df$mean[3],completeData.Result.df$mean[4],completeData.Result.df$mean[2]),
        group=factor(c('Completers Data','Completers Data','Completers Data','Completers Data',
                       'Complete Data','Complete Data','Complete Data','Complete Data')),
        VAL=c(round(MARcompleter.Result.df$mean[1],3),round(MARcompleter.Result.df$mean[3],3),round(MARcompleter.Result.df$mean[4],3),round(MARcompleter.Result.df$mean[2],3),
              round(completeData.Result.df$mean[1],3),round(completeData.Result.df$mean[3],3),round(completeData.Result.df$mean[4],3),round(completeData.Result.df$mean[2],3))
      )
      p.completer = ggplot(df.completer.MAR, aes(x = boxOdds, y = yAxis, colour = group))
      p.completer.mar <-p.completer  +
        geom_point(aes(x = boxOdds, y = yAxis),size = 2.5,  position = position_dodgev(height=0.9))+
        geom_text(aes(x = boxOdds, y = yAxis,label=VAL),size=3,vjust=-0.5,position = position_dodgev(height=0.9))+
        theme_bw()+
        theme(panel.grid.minor = element_blank(),panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"))+
        theme(plot.title = element_text(hjust = 0.5,size=13))+
        scale_x_continuous(breaks = c(1,seq(0,100,20)))+
        coord_trans(x = 'log10')+
        ylab("") +
        xlab("Odds ratio (log scale)") +
        ggtitle("Odds Ratio: Completers")
      p.completer.mar
    }
  })
  ############### output imputed #########
  output$imputed <- renderPlot({
    if((input$MCAR == TRUE) & (input$MAR == FALSE)){
      subjid_sample<-sample(1:1831,1831*MCARpercentage())
      imputedData <- adpa %>%
        filter(PARAMN == 10 & !(AVISITN == 2)) %>%
        mutate(PCHGCA1N = ifelse(SUBJID %in% subjid_sample,0,PCHGCA1N))%>%
        left_join(select(adsl,SUBJID,SEX),by="SUBJID")  %>%
        select(c(1:20,29)) %>%
        filter(AVISITN == 6)
      model3<-glm(PCHGCA1N~SEX+TRTP, family = binomial(), data = imputedData)
      coef.model3 <- exp(coefficients(model3))
      ORcompletersSexM2F <- as.numeric(round(coef.model3[2],4)) # OR male to female
      ORcompletersAvsP <- as.numeric(round(1/coef.model3[3],4)) # OR active control vs placebo
      ORcompleters140vsP <- as.numeric(round(coef.model3[4]/coef.model3[3],4)) # OR 140mg vs placebo
      ORcompleters210vsP <- as.numeric(round(coef.model3[5]/coef.model3[3],4)) # OR 140mg vs placebo
      OR10Pimputed <- list(SexM2F = ORcompletersSexM2F,AvsP = ORcompletersAvsP,T140vsP = ORcompleters140vsP, T210vsP =ORcompleters210vsP)
      # threeWayTable10PCompleters<-ftable(table(completers10PData$TRTPN,completers10PData$SEX,completers10PData$PCHGCA1N))
      imputedData.Result.df <- data.frame(mean = c(ORcompletersAvsP,ORcompletersSexM2F,ORcompleters140vsP,ORcompleters210vsP))
      boxLabels=factor(c(1,2,3,4,1,2,3,4))
      boxLabels=c('TRTP Active control \n vs Placebo','TRTP Test drug 140mg \n vs Placebo ','TRTP Test drug 210mg \n vs Placebo ','Sex(M vs F)')
      # plot 1: MCAR
      df.imputed = data.frame(
        yAxis = boxLabels,
        boxOdds = c(imputedData.Result.df$mean[1],imputedData.Result.df$mean[3],imputedData.Result.df$mean[4],imputedData.Result.df$mean[2],
                    completeData.Result.df$mean[1],completeData.Result.df$mean[3],completeData.Result.df$mean[4],completeData.Result.df$mean[2]),
        group=factor(c('Imputed Data (Composite)','Imputed Data (Composite)','Imputed Data (Composite)','Imputed Data (Composite)',
                       'Complete Data','Complete Data','Complete Data','Complete Data')),
        VAL=c(round(imputedData.Result.df$mean[1],3),round(imputedData.Result.df$mean[3],3),round(imputedData.Result.df$mean[4],3),round(imputedData.Result.df$mean[2],3),
              round(completeData.Result.df$mean[1],3),round(completeData.Result.df$mean[3],3),round(completeData.Result.df$mean[4],3),round(completeData.Result.df$mean[2],3))
      )
      p.imputed = ggplot(df.imputed, aes(x = boxOdds, y = yAxis, colour = group))
      p.imputed.mcar <-p.imputed  +
        geom_point(aes(x = boxOdds, y = yAxis),size = 2.5,  position = position_dodgev(height=0.9))+
        geom_text(aes(x = boxOdds, y = yAxis,label=VAL),size=3,vjust=-0.5,position = position_dodgev(height=0.9))+
        theme_bw()+
        theme(panel.grid.minor = element_blank(),panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"))+
        theme(plot.title = element_text(hjust = 0.5,size=13))+
        scale_x_continuous(breaks = c(1,seq(0,100,20)))+
        coord_trans(x = 'log10')+
        ylab("") +
        xlab("Odds ratio (log scale)") +
        ggtitle("Odds Ratio: Imputed Data (Composite Strategy), MCAR")
      p.imputed.mcar
    }
    else if((input$MCAR == FALSE) & (input$MAR == TRUE)){
      subjidMARvisit4high<- adpa %>%
        filter(AVISITN == 3 & PARAMN == 10 & PCHG < 10) %>%
        pull(SUBJID)
      subjidMARvisit4high_sample <- sample(subjidMARvisit4high,MARpercentage1()*length(subjidMARvisit4high))
      # subjects with 5% missing visit4
      allSubjid <- 1:1831
      subjidMARvisit4low <- setdiff(allSubjid,subjidMARvisit4high) 
      subjidMARvisit4low_sample <- sample(subjidMARvisit4low,MARpercent2()*length(subjidMARvisit4low))
      # set the selected subjects' PASI score to NA on visit 4 and on subsequent visits
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
      # subjects with 30% missing visit5
      subjidMARvisit5high<- adpa_mar %>%
        filter(AVISITN == 4 & PARAMN == 10 & PCHG < 10) %>%
        pull(SUBJID)
      subjidMARvisit5high_sample <- sample(subjidMARvisit5high,MARpercentage1()*length(subjidMARvisit5high))
      # subjects with 5% missing visit5
      subjidMARvisit5low <- setdiff(allSubjid,subjidMARvisit5high) 
      subjidMARvisit5low_sample <- sample(subjidMARvisit5low,MARpercent2()*length(subjidMARvisit5low))
      # set the selected subjects' PASI score to NA on visit 5 and 6
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
      subjidMARvisit6high_sample <- sample(subjidMARvisit6high,MARpercentage1()*length(subjidMARvisit6high))
      # subjects with 5% missing visit5
      subjidMARvisit6low <- setdiff(allSubjid,subjidMARvisit6high) 
      subjidMARvisit6low_sample <- sample(subjidMARvisit6low,MARpercent2()*length(subjidMARvisit6low))
      adpa_mar <- adpa_mar %>% 
        mutate(AVAL = ifelse(SUBJID %in% subjidMARvisit6high_sample & PARAMN == 10 & AVISITN == 6,NA,AVAL)) %>%
        mutate(PCHG = ifelse(SUBJID %in% subjidMARvisit6high_sample & PARAMN == 10 & AVISITN ==6,NA,PCHG)) %>%
        mutate(AVAL = ifelse(SUBJID %in% subjidMARvisit6low_sample & PARAMN == 10 & AVISITN == 6,NA,AVAL)) %>%
        mutate(PCHG = ifelse(SUBJID %in% subjidMARvisit6low_sample & PARAMN == 10 & AVISITN == 6,NA,PCHG)) %>%
        mutate(PCHGCA1N = ifelse(SUBJID %in% subjidMARvisit6low_sample & PARAMN == 10 & AVISITN == 6,NA,PCHGCA1N)) %>%
        mutate(PCHGCA1N = ifelse(SUBJID %in% subjidMARvisit6high_sample & PARAMN == 10 & AVISITN == 6,NA,PCHGCA1N)) %>%
        select(c(1:20,29,30))%>%
        filter(PARAMN == 10 & AVISITN == 6) 
      # Composite
      adpa_mar_imputed <-  adpa_mar %>% mutate(PCHGCA1N=ifelse(is.na(PCHGCA1N),0,PCHGCA1N))
      model3<-glm(PCHGCA1N~SEX+TRTP, family = binomial(), data = adpa_mar_imputed)
      coef.model3 <- exp(coefficients(model3))
      ORcompletersSexM2F <- as.numeric(round(coef.model3[2],4)) # OR male to female
      ORcompletersAvsP <- as.numeric(round(1/coef.model3[3],4)) # OR active control vs placebo
      ORcompleters140vsP <- as.numeric(round(coef.model3[4]/coef.model3[3],4)) # OR 140mg vs placebo
      ORcompleters210vsP <- as.numeric(round(coef.model3[5]/coef.model3[3],4)) # OR 140mg vs placebo
      #ORmarC <- list(SexM2F = ORcompletersSexM2F,AvsP = ORcompletersAvsP,T140vsP = ORcompleters140vsP, T210vsP =ORcompleters210vsP)
      MARimputed.Result.df <- data.frame(mean = c(ORcompletersAvsP,ORcompletersSexM2F,ORcompleters140vsP,ORcompleters210vsP))
      boxLabels=factor(c(1,2,3,4,1,2,3,4))
      boxLabels=c('TRTP Active control \n vs Placebo','TRTP Test drug 140mg \n vs Placebo ','TRTP Test drug 210mg \n vs Placebo ','Sex(M vs F)')
      # plot 2: imputed MCAR
      df.imputed.MAR = data.frame(
        yAxis = boxLabels,
        boxOdds = c(MARimputed.Result.df$mean[1],MARimputed.Result.df$mean[3],MARimputed.Result.df$mean[4],MARimputed.Result.df$mean[2],
                    completeData.Result.df$mean[1],completeData.Result.df$mean[3],completeData.Result.df$mean[4],completeData.Result.df$mean[2]),
        group=factor(c('Imputed Data (Composite)','Imputed Data (Composite)','Imputed Data (Composite)','Imputed Data (Composite)',
                       'Complete Data','Complete Data','Complete Data','Complete Data')),
        VAL=c(round(MARimputed.Result.df$mean[1],3),round(MARimputed.Result.df$mean[3],3),round(MARimputed.Result.df$mean[4],3),round(MARimputed.Result.df$mean[2],3),
              round(completeData.Result.df$mean[1],3),round(completeData.Result.df$mean[3],3),round(completeData.Result.df$mean[4],3),round(completeData.Result.df$mean[2],3))
      )
      p.imputed = ggplot(df.imputed.MAR, aes(x = boxOdds, y = yAxis, colour = group))
      p.imputed.mar <-p.imputed +
        geom_point(aes(x = boxOdds, y = yAxis),size = 2.5,  position = position_dodgev(height=0.9))+
        geom_text(aes(x = boxOdds, y = yAxis,label=VAL),size=3,vjust=-0.5,position = position_dodgev(height=0.9))+
        theme_bw()+
        theme(panel.grid.minor = element_blank(),panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"))+
        theme(plot.title = element_text(hjust = 0.5,size=13))+
        scale_x_continuous(breaks = c(1,seq(0,100,20)))+
        coord_trans(x = 'log10')+
        ylab("") +
        xlab("Odds ratio (log scale)") +
        ggtitle("Odds Ratio: Imputed Data (Composite Strategy), MAR")
      p.imputed.mar
    }
  })
  ############### output LOCF #########
  output$LOCF <- renderPlot({
    if((input$MCAR == TRUE) & (input$MAR == FALSE)){
      subjid_sample<-sample(1:1831,1831*MCARpercentage())
      ########### MCAR: Imputed (LOCF) ####
      imputed10PLOCFData <- adpa %>%
        filter(PARAMN == 10 & !(AVISITN == 2)) %>%
        group_by(SUBJID) %>%
        mutate(LAGPCHGCA1N = lag(PCHGCA1N,n=1,order_by = SUBJID)) %>%
        left_join(select(adsl,SUBJID,SEX),by="SUBJID") %>%
        mutate(PCHGCA1N = ifelse(SUBJID %in% subjid_sample,LAGPCHGCA1N,PCHGCA1N)) %>%
        select(c(1:20,29,30)) %>%
        filter(AVISITN == 6)
      # threeWayTable10PLOCF<-ftable(table(imputed10PLOCFData$TRTPN,imputed10PLOCFData$SEX,imputed10PLOCFData$PCHGCA1N))
      model3<-glm(PCHGCA1N~SEX+TRTP, family = binomial(), data = imputed10PLOCFData)
      coef.model3 <- exp(coefficients(model3))
      ORcompletersSexM2F <- as.numeric(round(coef.model3[2],4)) # OR male to female
      ORcompletersAvsP <- as.numeric(round(1/coef.model3[3],4)) # OR active control vs placebo
      ORcompleters140vsP <- as.numeric(round(coef.model3[4]/coef.model3[3],4)) # OR 140mg vs placebo
      ORcompleters210vsP <- as.numeric(round(coef.model3[5]/coef.model3[3],4)) # OR 140mg vs placebo
      OR10PLOCF <- list(SexM2F = ORcompletersSexM2F,AvsP = ORcompletersAvsP,T140vsP = ORcompleters140vsP, T210vsP =ORcompleters210vsP)
      LOCFData.Result.df <- data.frame(mean = c(ORcompletersAvsP,ORcompletersSexM2F,ORcompleters140vsP,ORcompleters210vsP))
      boxLabels=factor(c(1,2,3,4,1,2,3,4))
      boxLabels=c('TRTP Active control \n vs Placebo','TRTP Test drug 140mg \n vs Placebo ','TRTP Test drug 210mg \n vs Placebo ','Sex(M vs F)')
      # plot 3: LOCF MCAR
      df.LOCF = data.frame(
        yAxis = boxLabels,
        boxOdds = c(LOCFData.Result.df$mean[1],LOCFData.Result.df$mean[3],LOCFData.Result.df$mean[4],LOCFData.Result.df$mean[2],
                    completeData.Result.df$mean[1],completeData.Result.df$mean[3],completeData.Result.df$mean[4],completeData.Result.df$mean[2]),
        group=factor(c('Imputed Data (LOCF)','Imputed Data (LOCF)','Imputed Data (LOCF)','Imputed Data (LOCF)',
                       'Complete Data','Complete Data','Complete Data','Complete Data')),
        VAL=c(round(LOCFData.Result.df$mean[1],3),round(LOCFData.Result.df$mean[3],3),round(LOCFData.Result.df$mean[4],3),round(LOCFData.Result.df$mean[2],3),
              round(completeData.Result.df$mean[1],3),round(completeData.Result.df$mean[3],3),round(completeData.Result.df$mean[4],3),round(completeData.Result.df$mean[2],3))
      )
      p.LOCF = ggplot(df.LOCF, aes(x = boxOdds, y = yAxis, colour = group))
      p.LOCF.mcar <-p.LOCF  +
        geom_point(aes(x = boxOdds, y = yAxis),size = 2.5,  position = position_dodgev(height=0.9))+
        geom_text(aes(x = boxOdds, y = yAxis,label=VAL),size=3,vjust=-0.5,position = position_dodgev(height=0.9))+
        theme_bw()+
        theme(panel.grid.minor = element_blank(),panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"))+
        theme(plot.title = element_text(hjust = 0.5,size=13))+
        scale_x_continuous(breaks = c(1,seq(0,100,20)))+
        coord_trans(x = 'log10')+
        ylab("") +
        xlab("Odds ratio (log scale)") +
        ggtitle("Odds Ratio: LOCF Data (Composite Strategy), MCAR")
      p.LOCF.mcar
    }
    else if((input$MCAR == FALSE) & (input$MAR == TRUE)){
      subjidMARvisit4high<- adpa %>%
        filter(AVISITN == 3 & PARAMN == 10 & PCHG < 10) %>%
        pull(SUBJID)
      subjidMARvisit4high_sample <- sample(subjidMARvisit4high,MARpercentage1()*length(subjidMARvisit4high))
      # subjects with 5% missing visit4
      allSubjid <- 1:1831
      subjidMARvisit4low <- setdiff(allSubjid,subjidMARvisit4high) 
      subjidMARvisit4low_sample <- sample(subjidMARvisit4low,MARpercent2()*length(subjidMARvisit4low))
      # set the selected subjects' PASI score to NA on visit 4 and on subsequent visits
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
      # subjects with 30% missing visit5
      subjidMARvisit5high<- adpa_mar %>%
        filter(AVISITN == 4 & PARAMN == 10 & PCHG < 10) %>%
        pull(SUBJID)
      subjidMARvisit5high_sample <- sample(subjidMARvisit5high,MARpercentage1()*length(subjidMARvisit5high))
      # subjects with 5% missing visit5
      subjidMARvisit5low <- setdiff(allSubjid,subjidMARvisit5high) 
      subjidMARvisit5low_sample <- sample(subjidMARvisit5low,MARpercent2()*length(subjidMARvisit5low))
      # set the selected subjects' PASI score to NA on visit 5 and 6
      adpa_mar <- adpa_mar %>% 
        mutate(AVAL = ifelse(SUBJID %in% subjidMARvisit5high_sample & PARAMN == 10 & AVISITN %in% c(5,6),NA,AVAL)) %>%
        mutate(PCHG = ifelse(SUBJID %in% subjidMARvisit5high_sample & PARAMN == 10 & AVISITN %in% c(5,6),NA,PCHG)) %>%
        mutate(AVAL = ifelse(SUBJID %in% subjidMARvisit5low_sample & PARAMN == 10 & AVISITN %in% c(5,6),NA,AVAL)) %>%
        mutate(PCHG = ifelse(SUBJID %in% subjidMARvisit5low_sample & PARAMN == 10 & AVISITN %in% c(5,6),NA,PCHG)) %>%
        mutate(PCHGCA1N = ifelse(SUBJID %in% subjidMARvisit5low_sample & PARAMN == 10 & AVISITN %in% c(5,6),NA,PCHGCA1N)) %>%
        mutate(PCHGCA1N = ifelse(SUBJID %in% subjidMARvisit5high_sample & PARAMN == 10 & AVISITN %in% c(5,6),NA,PCHGCA1N)) 
      #LOCF
      adpa_mar_imputed_LOCF <- adpa_mar %>%
        mutate(PCHGCA1N=ifelse(is.na(PCHGCA1N) & AVISIT == 3,0,PCHGCA1N)) %>%
        mutate(PCHGCA1N = ifelse(is.na(PCHGCA1N) & !(AVISIT == 3),LAGPCHGCA1N,PCHGCA1N))
      model3<-glm(PCHGCA1N~SEX+TRTP, family = binomial(), data = adpa_mar_imputed_LOCF)
      coef.model3 <- exp(coefficients(model3))
      ORcompletersSexM2F <- as.numeric(round(coef.model3[2],4)) # OR male to female
      ORcompletersAvsP <- as.numeric(round(1/coef.model3[3],4)) # OR active control vs placebo
      ORcompleters140vsP <- as.numeric(round(coef.model3[4]/coef.model3[3],4)) # OR 140mg vs placebo
      ORcompleters210vsP <- as.numeric(round(coef.model3[5]/coef.model3[3],4)) # OR 140mg vs placebo
      MARLOCF.Result.df <- data.frame(mean = c(ORcompletersAvsP,ORcompletersSexM2F,ORcompleters140vsP,ORcompleters210vsP))
      boxLabels=factor(c(1,2,3,4,1,2,3,4))
      boxLabels=c('TRTP Active control \n vs Placebo','TRTP Test drug 140mg \n vs Placebo ','TRTP Test drug 210mg \n vs Placebo ','Sex(M vs F)')
      # plot 3: LOCF MCAR
      df.LOCF.MAR = data.frame(
        yAxis = boxLabels,
        boxOdds = c(MARLOCF.Result.df$mean[1],MARLOCF.Result.df$mean[3],MARLOCF.Result.df$mean[4],MARLOCF.Result.df$mean[2],
                    completeData.Result.df$mean[1],completeData.Result.df$mean[3],completeData.Result.df$mean[4],completeData.Result.df$mean[2]),
        group=factor(c('Imputed Data (LOCF)','Imputed Data (LOCF)','Imputed Data (LOCF)','Imputed Data (LOCF)',
                       'Complete Data','Complete Data','Complete Data','Complete Data')),
        VAL=c(round(MARLOCF.Result.df$mean[1],3),round(MARLOCF.Result.df$mean[3],3),round(MARLOCF.Result.df$mean[4],3),round(MARLOCF.Result.df$mean[2],3),
              round(completeData.Result.df$mean[1],3),round(completeData.Result.df$mean[3],3),round(completeData.Result.df$mean[4],3),round(completeData.Result.df$mean[2],3))
      )
      p.LOCF = ggplot(df.LOCF.MAR, aes(x = boxOdds, y = yAxis, colour = group))
      p.LOCF.mar <-p.LOCF +
        geom_point(aes(x = boxOdds, y = yAxis),size = 2.5,  position = position_dodgev(height=0.9))+
        geom_text(aes(x = boxOdds, y = yAxis,label=VAL),size=3,vjust=-0.5,position = position_dodgev(height=0.9))+
        theme_bw()+
        theme(panel.grid.minor = element_blank(),panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"))+
        theme(plot.title = element_text(hjust = 0.5,size=13))+
        scale_x_continuous(breaks = c(1,seq(0,100,20)))+
        coord_trans(x = 'log10')+
        ylab("") +
        xlab("Odds ratio (log scale)") +
        ggtitle("Odds Ratio: Imputed Data (LOCF Strategy), MAR")
      p.LOCF.mar
    }
   })
  observe({  
    if ( (input$MCAR == TRUE) & (input$MAR == FALSE)) {
      print("MCAR only")
      updateTabsetPanel(session, inputId="analysisTabs", selected="panel_complete_data")
      updateTabsetPanel(session, inputId="analysisTabs", selected="panel_completers")
      updateTabsetPanel(session, inputId="analysisTabs", selected="panel_composite")
      updateTabsetPanel(session, inputId="analysisTabs", selected="panel_LOCF")
    } 
    else if((input$MCAR == FALSE) & (input$MAR == TRUE)) {
      print("MAR only")
      updateTabsetPanel(session, inputId="analysisTabs", selected="panel_complete_data")
      updateTabsetPanel(session, inputId="analysisTabs", selected="panel_completers")
      updateTabsetPanel(session, inputId="analysisTabs", selected="panel_composite")
      updateTabsetPanel(session, inputId="analysisTabs", selected="panel_LOCF")
    }
  })#end observe
})
shinyApp(ui = ui, server = server)
