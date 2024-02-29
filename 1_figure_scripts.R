# This is the code that produced the figures and tables for the paper "Temporal trends in test-seeking behaviour during the COVID-19 pandemic"
# The code relies on input data that is not publicly available (lines 24-31) and so will crash if you try to run it.

library(ggplot2)
library(mgcv)
library(epitools)
library(cowplot)
library(grid)
library(dplyr)
################################################################################
## Read in data ##
# This data is not publicly available
# This data is not publicly available
# This data is not publicly available
# This data is not publicly available
# This data is not publicly available
# This data is not publicly available
# This data is not publicly available
# This data is not publicly available

# Survey data
sdat <- read.csv('cleaned_data/survey_data_comb.csv')

# Case data
cases <- read.csv("data/cases data/local_cases_input_2023-10-05.csv")

# Flutracking data
flu <- read.csv('cleaned_data/flutrack_data.csv')

################################################################################
## Focus on specific period

# Set date limit
mindate <- as.Date("2021-11-02")
maxdate <- as.Date("2023-09-03")

# Remove data out of bounds of limits
cases <- cases[cases$date_onset >= mindate & cases$date_onset <= maxdate,]
sdat <- sdat[sdat$date >= mindate & sdat$date <= maxdate,]
flu <- flu[flu$SurveyWeek >= mindate & flu$SurveyWeek <=maxdate,]

################################################################################
week_starts <- seq(as.Date("2021-11-02")+27, as.Date("2023-09-03"), by=1)
sympt_cat <- "at least one core symptom"



sdat$response2 <- ifelse(sdat$test=="No", 0,
                         ifelse(sdat$test_symp=="Yes", 1, 0))

################################################################################
## Functions

reformat_data <- function(tdat, symptom_cat, weeks){
  
  tdat <- tdat[tdat$symptoms == symptom_cat & tdat$response=="Yes",]
  
  test <- table(tdat$date, tdat$response2)
  tab_df <- data.frame(neg = as.numeric(test[,1]),
                       pos= as.numeric(test[,2]),
                       date = rownames(test))
  
  
  df <- data.frame()
  for(i in seq_len(length(weeks))){
    print(i)
    temp<- tab_df[tab_df$date>=weeks[i]-27 & tab_df$date<=weeks[i],]
    
    
    neg <- sum(temp$neg)
    pos <- sum(temp$pos)
    if(pos+neg > 0){
      
      p <- binom.exact(pos, pos+neg)
      
    } else{
      p <- data.frame(proportion=NA,
                      lower=NA,
                      upper=NA)
    }
    
    row_df <- data.frame(p = p$proportion,
                         lwr = p$lower,
                         upr = p$upper,
                         week = i,
                         date = weeks[i])
    
    df <- rbind(df, row_df)
    
  }
  
  return(df)
}



get_tab_rtt <- function(tdat, sympt_cat){
  tdat <- tdat[tdat$symptoms == sympt_cat & tdat$response=="Yes",]
  
  tdat$response3 <- ifelse(tdat$test=="No", 0,
                           ifelse(tdat$test_contact=="Yes", 1, 0))
  test1 <- table(tdat$response3)
  
  tdat$response3 <- ifelse(tdat$test=="No", 0,
                           ifelse(tdat$test_other=="Yes", 1, 0))
  test2 <- table(tdat$response3)
  
  tdat$response3 <- ifelse(tdat$test=="No", 0,
                           ifelse(tdat$test_symp=="Yes", 1, 0))
  test3 <- table(tdat$response3)
  
  tdat$response3 <- ifelse(tdat$test=="No", 0,
                           ifelse(tdat$test_job=="Yes", 1, 0))
  test4 <- table(tdat$response3)
  
  
  tdat$response3 <- ifelse(tdat$test=="No", 0,
                           ifelse(tdat$test_job=="Yes", 1, 
                                  ifelse(tdat$test_symp=="Yes", 1, 
                                         ifelse(tdat$test_other=="Yes", 1, 
                                                ifelse(tdat$test_contact=="Yes", 1, 0)))))
  
  test5 <- table(tdat$response3)
  
  bn1 <- binom.exact(test1[2], test1[1]+test1[2])
  bn2 <- binom.exact(test2[2], test2[1]+test2[2])
  bn3 <- binom.exact(test3[2], test3[1]+test3[2])
  bn4 <- binom.exact(test4[2], test4[1]+test4[2])
  bn5 <- binom.exact(test5[2], test5[1]+test5[2])
  bn1$lab <- "Contact"
  bn2$lab <- "Other"
  bn3$lab <- "Sympt"
  bn4$lab <- "Job"
  bn5$lab <- "All"
  
  df <- rbind(bn5, bn1, bn2, bn3, bn4)
  return(df)
  
}



get_tab_rtt_RATPCR <- function(tdat, sympt_cat){
  tdat <- tdat[tdat$symptoms == sympt_cat & tdat$response=="Yes",]
  
  tdat$response3 <- ifelse(tdat$test=="No", 0,
                           ifelse(tdat$test_contact=="Yes", 1, 0))
  test1 <- table(tdat$response3)
  
  tdat$response3 <- ifelse(tdat$test=="No", 0,
                           ifelse(tdat$test_other=="Yes", 1, 0))
  test2 <- table(tdat$response3)
  
  tdat$response3 <- ifelse(tdat$test=="No", 0,
                           ifelse(tdat$test_symp=="Yes", 1, 0))
  test3 <- table(tdat$response2)
  
  tdat$response3 <- ifelse(tdat$test=="No", 0,
                           ifelse(tdat$test_job=="Yes", 1, 0))
  test4 <- table(tdat$response3)
  
  bn1 <- binom.exact(test1[2], test1[1]+test1[2])
  bn2 <- binom.exact(test2[2], test2[1]+test2[2])
  bn3 <- binom.exact(test3[2], test3[1]+test3[2])
  bn4 <- binom.exact(test4[2], test4[1]+test4[2])
  bn1$lab <- "Contact"
  bn2$lab <- "Other"
  bn3$lab <- "Sympt"
  bn4$lab <- "Job"
  
  df <- rbind(bn1, bn2, bn3, bn4)
  return(df)
  
}




reformat_data_rat_pos <- function(tdat, symptom_cat, weeks){
  
  #tdat <- tdat[tdat$symptoms == symptom_cat & tdat$response=="Yes" &tdat$test_symp=="Yes",]
  tdat <- tdat[tdat$symptoms == symptom_cat & tdat$response=="Yes",]
  #tdat <- tdat[tdat$test_symp=="Yes" & is.na(tdat$test_symp)==FALSE,]
  tdat$rat_response <- ifelse(tdat$rat_result=="The test showed I have COVID-19 (positive result)", 1,
                              ifelse(tdat$rat_result=="The test showed I do not have COVID-19 (negative result)", 0, NA))
  
  
  test <- table(tdat$date, tdat$rat_response)
  tab_df <- data.frame(neg = as.numeric(test[,1]),
                       pos= as.numeric(test[,2]),
                       date = rownames(test))
  
  neg<-sum(tab_df$neg)
  pos<-sum(tab_df$pos)
  p <- binom.exact(pos, pos+neg)
  row_df <- data.frame(p = p$proportion,
                       lwr = p$lower,
                       upr = p$upper,
                       week = "Overall",
                       date = as.Date("2000-01-01"))
  
  df <- data.frame()
  df <- rbind(df, row_df)
  
  

  for(i in seq_len(length(weeks))){
    print(i)
    temp<- tab_df[tab_df$date>=weeks[i]-27 & tab_df$date<=weeks[i],]
    
    
    neg <- sum(temp$neg)
    pos <- sum(temp$pos)
    if(pos+neg > 0){
      
      p <- binom.exact(pos, pos+neg)
      
    } else{
      p <- data.frame(proportion=NA,
                      lower=NA,
                      upper=NA)
    }
    
    row_df <- data.frame(p = p$proportion,
                         lwr = p$lower,
                         upr = p$upper,
                         week = i,
                         date = weeks[i])
    
    df <- rbind(df, row_df)
    
  }
  
  return(df)
}




reformat_data_pcr_pos <- function(tdat, symptom_cat, weeks){
  
  #tdat <- tdat[tdat$symptoms == symptom_cat & tdat$response=="Yes" &tdat$test_symp=="Yes",]
  tdat <- tdat[tdat$symptoms == symptom_cat & tdat$response=="Yes",]
  #tdat <- tdat[tdat$test_symp=="Yes" & is.na(tdat$test_symp)==FALSE,]
  tdat$pcr_response <- ifelse(tdat$pcr_result=="The test showed I have COVID-19 (positive result)", 1,
                              ifelse(tdat$pcr_result=="The test showed I do not have COVID-19 (negative result)", 0, NA))
  
  
  test <- table(tdat$date, tdat$pcr_response)
  tab_df <- data.frame(neg = as.numeric(test[,1]),
                       pos= as.numeric(test[,2]),
                       date = rownames(test))
  
  neg<-sum(tab_df$neg)
  pos<-sum(tab_df$pos)
  p <- binom.exact(pos, pos+neg)
  row_df <- data.frame(p = p$proportion,
                       lwr = p$lower,
                       upr = p$upper,
                       week = "Overall",
                       date = as.Date("2000-01-01"))
  
  df <- data.frame()
  df <- rbind(df, row_df)
  
  
  #df <- data.frame()
  for(i in seq_len(length(weeks))){
    print(i)
    temp<- tab_df[tab_df$date>=weeks[i]-27 & tab_df$date<=weeks[i],]
    
    
    neg <- sum(temp$neg)
    pos <- sum(temp$pos)
    if(pos+neg > 0){
      
      p <- binom.exact(pos, pos+neg)
      
    } else{
      p <- data.frame(proportion=NA,
                      lower=NA,
                      upper=NA)
    }
    
    row_df <- data.frame(p = p$proportion,
                         lwr = p$lower,
                         upr = p$upper,
                         week = i,
                         date = weeks[i])
    
    df <- rbind(df, row_df)
    
  }
  
  return(df)
}



reformat_data_other_reason_to_test <- function(tdat, symptom_cat, weeks, cat){
  
  tdat <- tdat[tdat$symptoms == symptom_cat & tdat$response=="Yes",]
  
  if(cat=="contact"){
    tdat$response3 <- ifelse(tdat$test=="No", 0,
                             ifelse(tdat$test_contact=="Yes", 1, 0))
  } else if(cat=="other"){
    tdat$response3 <- ifelse(tdat$test=="No", 0,
                             ifelse(tdat$test_other=="Yes", 1, 0))
  } else if(cat == "sympt"){
    tdat$response3 <- ifelse(tdat$test=="No", 0,
                             ifelse(tdat$test_symp=="Yes", 1, 0))
  } else if(cat == "job"){
    tdat$response3 <- ifelse(tdat$test=="No", 0,
                             ifelse(tdat$test_job=="Yes", 1, 0))
  }
  
  
  test <- table(tdat$date, tdat$response3)
  tab_df <- data.frame(neg = as.numeric(test[,1]),
                       pos= as.numeric(test[,2]),
                       date = rownames(test))
  
  
  df <- data.frame()
  for(i in seq_len(length(weeks))){
    print(i)
    temp<- tab_df[tab_df$date>=weeks[i]-27 & tab_df$date<=weeks[i],]
    
    
    neg <- sum(temp$neg)
    pos <- sum(temp$pos)
    if(pos+neg > 0){
      
      p <- binom.exact(pos, pos+neg)
      
    } else{
      p <- data.frame(proportion=NA,
                      lower=NA,
                      upper=NA)
    }
    
    row_df <- data.frame(p = p$proportion,
                         lwr = p$lower,
                         upr = p$upper,
                         week = i,
                         date = weeks[i])
    
    df <- rbind(df, row_df)
    
  }
  
  return(df)
}


reformat_data_flu_test <- function(tdat, weeks){
  
  tdat <- tdat[tdat$IsCase == 1,]
  
  test <- table(tdat$SurveyWeek, tdat$TestCovid)
  tab_df <- data.frame(neg = as.numeric(test[,1]),
                       pos= as.numeric(test[,2]),
                       date = rownames(test))
  
  
  df <- data.frame()
  for(i in seq_len(length(weeks))){
    print(i)
    temp<- tab_df[tab_df$date>=weeks[i]-27 & tab_df$date<=weeks[i],]
    
    
    neg <- sum(temp$neg)
    pos <- sum(temp$pos)
    if(pos+neg > 0){
      
      p <- binom.exact(pos, pos+neg)
      
    } else{
      p <- data.frame(proportion=NA,
                      lower=NA,
                      upper=NA)
    }
    
    row_df <- data.frame(p = p$proportion,
                         lwr = p$lower,
                         upr = p$upper,
                         week = i,
                         date = weeks[i])
    
    df <- rbind(df, row_df)
    
  }
  
  return(df)
}


reformat_data_flu_SMC <- function(tdat, weeks){
  
  tdat <- tdat[tdat$IsCase == 1,]
  
  test <- table(tdat$SurveyWeek, tdat$SeekMedics)
  tab_df <- data.frame(neg = as.numeric(test[,1]),
                       pos= as.numeric(test[,2]),
                       date = rownames(test))
  
  
  df <- data.frame()
  for(i in seq_len(length(weeks))){
    print(i)
    temp<- tab_df[tab_df$date>=weeks[i]-27 & tab_df$date<=weeks[i],]
    
    
    neg <- sum(temp$neg)
    pos <- sum(temp$pos)
    if(pos+neg > 0){
      
      p <- binom.exact(pos, pos+neg)
      
    } else{
      p <- data.frame(proportion=NA,
                      lower=NA,
                      upper=NA)
    }
    
    row_df <- data.frame(p = p$proportion,
                         lwr = p$lower,
                         upr = p$upper,
                         week = i,
                         date = weeks[i])
    
    df <- rbind(df, row_df)
    
  }
  
  return(df)
}


reformat_data_ratres <- function(tdat, symptom_cat, weeks){
  
  tdat <- tdat[tdat$symptoms == symptom_cat & tdat$response=="Yes",]
  
  test <- table(tdat$date, tdat$rat_report)
  tab_df <- data.frame(neg = as.numeric(test[,1]),
                       pos= as.numeric(test[,2]),
                       date = rownames(test))
  
  
  df <- data.frame()
  for(i in seq_len(length(weeks))){
    print(i)
    temp<- tab_df[tab_df$date>=weeks[i]-27 & tab_df$date<=weeks[i],]
    
    
    neg <- sum(temp$neg)
    pos <- sum(temp$pos)
    if(pos+neg > 0){
      
      p <- binom.exact(pos, pos+neg)
      
    } else{
      p <- data.frame(proportion=NA,
                      lower=NA,
                      upper=NA)
    }
    
    row_df <- data.frame(p = p$proportion,
                         lwr = p$lower,
                         upr = p$upper,
                         week = i,
                         date = weeks[i])
    
    df <- rbind(df, row_df)
    
  }
  
  return(df)
}




reformat_data_what_test <- function(tdat, symptom_cat, weeks, w_test){
  
  tdat$response2 <- ifelse(tdat$test=="No", 0,
                           ifelse(tdat$test_symp=="Yes" & tdat$what_test == w_test, 1, 0))
  
  tdat <- tdat[tdat$symptoms == symptom_cat & tdat$response=="Yes",]
  
  test <- table(tdat$date, tdat$response2)
  tab_df <- data.frame(neg = as.numeric(test[,1]),
                       pos= as.numeric(test[,2]),
                       date = rownames(test))
  
  
  df <- data.frame()
  for(i in seq_len(length(weeks))){
    print(i)
    temp<- tab_df[tab_df$date>=weeks[i]-27 & tab_df$date<=weeks[i],]
    
    
    neg <- sum(temp$neg)
    pos <- sum(temp$pos)
    if(pos+neg > 0){
      
      p <- binom.exact(pos, pos+neg)
      
    } else{
      p <- data.frame(proportion=NA,
                      lower=NA,
                      upper=NA)
    }
    
    row_df <- data.frame(p = p$proportion,
                         lwr = p$lower,
                         upr = p$upper,
                         week = i,
                         date = weeks[i])
    
    df <- rbind(df, row_df)
    
  }
  
  return(df)
}


################################################################################
sympt1 <- reformat_data(sdat, symptom_cat = "at least one core symptom", weeks = week_starts)
sympt2 <- reformat_data(sdat, symptom_cat = "at least one symptom", weeks = week_starts)
sympt3 <- reformat_data(sdat, symptom_cat = "at least two core symptoms", weeks = week_starts)
sympt4 <- reformat_data(sdat, symptom_cat = "at least fever and cough", weeks = week_starts)

sympt1$lab <- "at least one core symptom"
sympt2$lab <- "at least one symptom"
sympt3$lab <- "at least two core symptoms"
sympt4$lab <- "at least fever and cough"

sympt_df <- rbind(sympt1, sympt2, sympt3, sympt4)

sympt_df$lab <- factor(sympt_df$lab, levels=c("at least one symptom", "at least one core symptom", "at least two core symptoms", "at least fever and cough"))

sympt_df1 <- sympt_df[sympt_df$date < as.Date("2021-12-20"),]
sympt_df2 <- sympt_df[sympt_df$date > as.Date("2021-12-20")+28,]


fig1a<-ggplot()+
  geom_line(data = sympt_df1,aes(x=date, y= p, color=lab))+
  geom_ribbon(data = sympt_df1,aes(x=date, y= p, ymin=lwr, ymax=upr, fill=lab), alpha = 0.2)+
  geom_line(data = sympt_df2, aes(x=date, y= p, color=lab))+
  geom_ribbon(data = sympt_df2,aes(x=date, y= p, ymin=lwr, ymax=upr, fill=lab), alpha = 0.2)+
  theme_bw()+
  scale_color_brewer("Symptom status",palette = "Dark2")+
  scale_fill_brewer("Symptom status", palette = "Dark2")+
  scale_x_date(date_breaks = "2 month", date_labels =  "%b\n%Y")+
  scale_y_continuous(breaks=c(0,0.2,0.4,0.6,0.8,1.0))+
  ylab("Probability of testing due to symptoms\n")+
  xlab("Date")+
  coord_cartesian(ylim=c(0,1), xlim=c(as.Date("2021-12-15"), as.Date("2023-08-15")))+
  theme(legend.position = c(0.83,0.83),
        legend.background = element_rect(color='black'))

fig1a

ggsave('figures/Figure1.pdf', width=8, height=6)


sympt_df1$min_date <- sympt_df1$date-27
sympt_df2$min_date <- sympt_df2$date-27
sympt_df2[sympt_df2$min_date < as.Date("2022-01-11"),]$min_date <- as.Date("2022-01-11")
sympt_df <- rbind(sympt_df1, sympt_df2)

write.csv(sympt_df, 'figures/figure1a.csv')

#########################################################################################################################################
state1 <- reformat_data(sdat[sdat$state=="NSW",], symptom_cat = "at least one core symptom", weeks = week_starts)
state2 <- reformat_data(sdat[sdat$state=="VIC",], symptom_cat = "at least one core symptom", weeks = week_starts)
state3 <- reformat_data(sdat[sdat$state=="QLD",], symptom_cat = "at least one core symptom", weeks = week_starts)
state4 <- reformat_data(sdat[sdat$state=="TAS",], symptom_cat = "at least one core symptom", weeks = week_starts)
state5 <- reformat_data(sdat[sdat$state=="ACT",], symptom_cat = "at least one core symptom", weeks = week_starts)
state6 <- reformat_data(sdat[sdat$state=="WA",], symptom_cat = "at least one core symptom", weeks = week_starts)
state7 <- reformat_data(sdat[sdat$state=="NT",], symptom_cat = "at least one core symptom", weeks = week_starts)
state8 <- reformat_data(sdat[sdat$state=="SA",], symptom_cat = "at least one core symptom", weeks = week_starts)

state1$lab <- "NSW"
state2$lab <- "VIC"
state3$lab <- "QLD"
state4$lab <- "TAS"
state5$lab <- "ACT"
state6$lab <- "WA"
state7$lab <- "NT"
state8$lab <- "SA"

state_df <- rbind(state1,state2,state3,state4,
                  state5,state6,state7,state8)


state_df1 <- state_df[state_df$date < as.Date("2021-12-20"),]
state_df2 <- state_df[state_df$date > as.Date("2021-12-20")+28,]#as.Date("2022-01-11")


#########################################################################################################################################
# Cases for description of Aus COVID-19 situation
cases$lab <- cases$state
cases$rate <- cases$infections
cases[cases$state=="VIC",]$rate <- cases[cases$state=="VIC",]$rate*100/6766.6
cases[cases$state=="NSW",]$rate <- cases[cases$state=="NSW",]$rate*100/8294.0
cases[cases$state=="QLD",]$rate <- cases[cases$state=="QLD",]$rate*100/5418.5
cases[cases$state=="TAS",]$rate <- cases[cases$state=="TAS",]$rate*100/572.7
cases[cases$state=="NT",]$rate <- cases[cases$state=="NT",]$rate*100/251.7
cases[cases$state=="SA",]$rate <- cases[cases$state=="SA",]$rate*100/1844.6
cases[cases$state=="WA",]$rate <- cases[cases$state=="WA",]$rate*100/2855.6
cases[cases$state=="ACT",]$rate <- cases[cases$state=="ACT",]$rate*100/464.6


fig1b<-ggplot()+
  geom_rect(aes(ymin=0, ymax=Inf,xmin=as.Date("2021-01-01"), xmax= as.Date("2021-12-20")), alpha=0.1)+
  geom_rect(aes(ymin=0, ymax=Inf,xmin=as.Date("2021-12-20"), xmax= as.Date("2022-02-28")), alpha=0.2)+
  geom_rect(aes(ymin=0, ymax=Inf,xmin=as.Date("2022-02-28"), xmax= as.Date("2022-06-20")), alpha=0.1)+
  geom_rect(aes(ymin=0, ymax=Inf,xmin=as.Date("2022-06-20"), xmax= as.Date("2022-11-07")), alpha=0.2)+
  geom_rect(aes(ymin=0, ymax=Inf,xmin=as.Date("2022-11-07"), xmax= as.Date("2023-12-20")), alpha=0.1)+
  geom_label(label="Delta", aes(y=10500,x=as.Date("2021-12-02")), size=2)+
  geom_label(label="BA.1", aes(y=10500,x=as.Date("2022-01-24")), size=2)+
  geom_label(label="BA.2", aes(y=10500,x=as.Date("2022-04-24")), size=2)+
  geom_label(label="BA.4/BA.5", aes(y=10500,x=as.Date("2022-08-30")), size=2)+
  geom_label(label="Assorted Omicron variants", aes(y=10500,x=as.Date("2023-04-15")), size=2)+
  geom_line(data=cases, aes(x=as.Date(date_onset), y=rate, color=state), size=0.1)+
  geom_point(data=cases, aes(x=as.Date(date_onset), y=rate, color=state), size=0.5)+
  theme_bw()+
  scale_color_brewer("State/territory",palette="Dark2")+
  scale_x_date(date_breaks = "2 month", date_labels =  "%b\n%Y")+
  scale_y_log10()+
  ylab("Daily cases per 100,000 people")+
  xlab("Date")+
  coord_cartesian(ylim=c(1,10000), xlim=c(as.Date("2021-12-15"), as.Date("2023-08-15")))+
  theme(legend.position = c(0.81,0.78),
        legend.background = element_rect(color='black'),
        panel.grid.minor.y = element_blank())+
  guides(color=guide_legend(ncol=4,byrow=TRUE))

fig1b

ggsave('figures/Figure2.pdf', width=8, height=6)
write.csv(cases[c("date_onset", "state","infections", "rate")],'figures/figure1b.csv')

fig1a <- fig1a + labs(tag="A")+
  theme(plot.tag.position=c(0.01,0.99),
        axis.text.x=element_blank(),
        axis.title.x = element_blank())


fig1b <- fig1b + labs(tag="B")+
  theme(plot.tag.position=c(0.01,0.99))

plot_grid(fig1a, fig1b, nrow=2, rel_heights=c(1,0.8))

ggsave('figures/Figure1AB.pdf', width=8, height=9)

##############################################################################################


ggplot()+
  geom_line(data = state_df1,aes(x=date, y= p))+
  geom_ribbon(data = state_df1,aes(x=date, y= p, ymin=lwr, ymax=upr), alpha = 0.2)+
  geom_line(data = state_df2, aes(x=date, y= p))+
  geom_ribbon(data = state_df2,aes(x=date, y= p, ymin=lwr, ymax=upr), alpha = 0.2)+
  geom_line(data=cases, aes(x=as.Date(date_onset), y=log10(rate)/5),color='red3', size=0.1)+
  geom_point(data=cases, aes(x=as.Date(date_onset), y=log10(rate)/5),color='red3', size=0.02)+
  scale_y_continuous("Probability of testing due to symptoms\n(at least one core symptom)", sec.axis = sec_axis(name="Daily cases per 100,000 people",~.*5, breaks = c(0,1,2,3,4,5,6), labels = c(expression(10^0),expression(10^1),expression(10^2),expression(10^3),expression(10^4),expression(10^5),expression(10^6))))+
  theme_bw()+
  facet_wrap(.~lab, nrow=4)+
  scale_x_date(date_breaks = "3 month", date_labels =  "%b\n%Y")+
  ylab("Probability of test given symptoms")+
  xlab("Date")+
  coord_cartesian(ylim=c(0,0.75))+
  theme(legend.position = "none",
        legend.background = element_rect(color='black'),
        panel.grid.minor.x = element_blank(),
        axis.line.y.right = element_line(color='red3'))

ggsave('figures/SFigure2.pdf', width=8, height=10)


ggplot()+
  geom_line(data = state_df1,aes(x=date, y= p))+
  geom_ribbon(data = state_df1,aes(x=date, y= p, ymin=lwr, ymax=upr), alpha = 0.2)+
  geom_line(data = state_df2, aes(x=date, y= p))+
  geom_ribbon(data = state_df2,aes(x=date, y= p, ymin=lwr, ymax=upr), alpha = 0.2)+
  geom_line(data=cases, aes(x=as.Date(date_onset)+14, y=log10(rate)/5),color='red3', size=0.1)+
  geom_point(data=cases, aes(x=as.Date(date_onset)+14, y=log10(rate)/5),color='red3', size=0.02)+
  scale_y_continuous("Probability of testing because of symptoms\n(at least one core symptom)", sec.axis = sec_axis(name="Daily cases per 100,000 people",~.*5, breaks = c(0,1,2,3,4,5,6), labels = c(expression(10^0),expression(10^1),expression(10^2),expression(10^3),expression(10^4),expression(10^5),expression(10^6))))+
  theme_bw()+
  facet_wrap(.~lab, nrow=4)+
  scale_x_date(date_breaks = "3 month", date_labels =  "%b\n%Y")+
  ylab("Probability of test given symptoms")+
  xlab("Date")+
  coord_cartesian(ylim=c(0,0.75))+
  theme(legend.position = "none",
        legend.background = element_rect(color='black'),
        panel.grid.minor.x = element_blank(),
        axis.line.y.right = element_line(color = "red3"))

ggsave('figures/SFigure2_lagged.pdf', width=8, height=10)


state_df1$min_date <- state_df1$date-27
state_df2$min_date <- state_df2$date-27
state_df2[state_df2$min_date < as.Date("2022-01-11"),]$min_date <- as.Date("2022-01-11")
state_df <- rbind(state_df1, state_df2)

write.csv(state_df, 'figures/Sfigure2.csv')


#########################################################################################################################################
# States grouped

state_g1 <- reformat_data(sdat[sdat$state %in% c("NSW", "ACT", "VIC"),], symptom_cat = "at least one core symptom", weeks = week_starts)
state_g2 <- reformat_data(sdat[sdat$state%in%c("QLD", "NT", "SA", "TAS"),], symptom_cat = "at least one core symptom", weeks = week_starts)
state_g3 <- reformat_data(sdat[sdat$state=="WA",], symptom_cat = "at least one core symptom", weeks = week_starts)


state_g1$lab <- "ACT/NSW/VIC"
state_g2$lab <- "QLD/NT/SA/TAS"
state_g3$lab <- "WA"


state_gdf <- rbind(state_g1,state_g2,state_g3)


state_gdf1 <- state_gdf[state_gdf$date < as.Date("2021-12-20"),]
state_gdf2 <- state_gdf[state_gdf$date > as.Date("2021-12-20")+28,]

plt_as_a<-ggplot()+
  geom_line(data = state_gdf1,aes(x=date, y= p, color=lab))+
  geom_ribbon(data = state_gdf1,aes(x=date, y= p, ymin=lwr, ymax=upr, fill=lab), alpha = 0.2)+
  geom_line(data = state_gdf2, aes(x=date, y= p, color=lab))+
  geom_ribbon(data = state_gdf2,aes(x=date, y= p, ymin=lwr, ymax=upr, fill=lab), alpha = 0.2)+
  theme_bw()+
  scale_color_brewer("State/territory",palette = "Dark2")+
  scale_fill_brewer("State/territory", palette = "Dark2")+
  scale_x_date(date_breaks = "2 month", date_labels =  "%b\n%Y")+
  ylab("Probability of testing due to symptoms\n(at least one core symptom)")+
  coord_cartesian(ylim=c(0, 0.75), xlim=c(as.Date("2021-12-15"), as.Date("2023-08-15")))+
  theme(legend.position = c(0.9,0.8),
        legend.background = element_rect(color='black'))



plt_as_a2<- plt_as_a

#########################################################################################################################################
# ages grouped

age_1 <- reformat_data(sdat[sdat$age<30,], symptom_cat = "at least one core symptom", weeks = week_starts)
age_2 <- reformat_data(sdat[sdat$age<60&sdat$age>=30,], symptom_cat = "at least one core symptom", weeks = week_starts)
#age_3 <- reformat_data(sdat[sdat$age<60&sdat$age>=45,], symptom_cat = "at least one core symptom", weeks = week_starts)
age_3 <- reformat_data(sdat[sdat$age>=60,], symptom_cat = "at least one core symptom", weeks = week_starts)


age_1$lab <- "18-29"
age_2$lab <- "30-59"
#age_3$lab <- "45-59"
age_3$lab <- "60+"


age_df <- rbind(age_1, age_2, age_3)


age_df1 <- age_df[age_df$date < as.Date("2021-12-20"),]
age_df2 <- age_df[age_df$date > as.Date("2021-12-20")+28,]

plt_as_b<-ggplot()+
  geom_line(data = age_df1,aes(x=date, y= p, color=lab))+
  geom_ribbon(data = age_df1,aes(x=date, y= p, ymin=lwr, ymax=upr, fill=lab), alpha = 0.2)+
  geom_line(data = age_df2, aes(x=date, y= p, color=lab))+
  geom_ribbon(data = age_df2,aes(x=date, y= p, ymin=lwr, ymax=upr, fill=lab), alpha = 0.2)+
  theme_bw()+
  scale_color_brewer("Ages",palette = "Dark2")+
  scale_fill_brewer("Ages", palette = "Dark2")+
  scale_x_date(date_breaks = "2 month", date_labels =  "%b\n%Y")+
  ylab("Probability of testing due to symptoms\n(at least one core symptom)")+
  xlab("Date")+
  coord_cartesian(ylim=c(0, 0.75), xlim=c(as.Date("2021-12-15"), as.Date("2023-08-15")))+
  theme(legend.position = c(0.9,0.8),
        legend.background = element_rect(color='black'))


plt_as_a <- plt_as_a+labs(tag="A")+
  theme(plot.tag.position = c(0.02, 0.98),
        axis.title.x = element_blank(),
        axis.text.x = element_blank())

plt_as_b <- plt_as_b+labs(tag="B")+
  theme(plot.tag.position = c(0.02, 0.98))

plot_grid(plt_as_a, plt_as_b, nrow=2, rel_heights = c(1,1.1))
ggsave('figures/Figure3.pdf', width=10, height=8)


state_gdf1$min_date <- state_gdf1$date-27
state_gdf2$min_date <- state_gdf2$date-27
state_gdf2[state_gdf2$min_date < as.Date("2022-01-11"),]$min_date <- as.Date("2022-01-11")
state_gdf <- rbind(state_gdf1, state_gdf2)

write.csv(state_gdf, 'figures/figure3a.csv')

age_df1$min_date <- age_df1$date-27
age_df2$min_date <- age_df2$date-27
age_df2[age_df2$min_date < as.Date("2022-01-11"),]$min_date <- as.Date("2022-01-11")
age_df <- rbind(age_df1, age_df2)

write.csv(age_df, 'figures/figure3b.csv')

#########################################################################################################################################
# What test / rat reported?

test_1 <- reformat_data_what_test(sdat, symptom_cat = "at least one core symptom", weeks = week_starts, w_test="Both")
test_2 <- reformat_data_what_test(sdat, symptom_cat = "at least one core symptom", weeks = week_starts, w_test="PCR only")
test_3 <- reformat_data_what_test(sdat, symptom_cat = "at least one core symptom", weeks = week_starts, w_test="RAT only")


test_1$lab <- "Both"
test_2$lab <- "PCR only"
test_3$lab <- "RAT only"

test_1 <- test_1[test_1$date >= as.Date("2022-01-11"),]
test_3 <- test_3[test_3$date >= as.Date("2022-01-11"),]


test_df <- rbind(test_1, test_2, test_3)


test_df1 <- test_df[test_df$date < as.Date("2021-12-20"),]
test_df2 <- test_df[test_df$date > as.Date("2021-12-20")+28,]

test_plt_a <- ggplot()+
  geom_line(data = test_df1,aes(x=date, y= p, color=lab))+
  geom_ribbon(data = test_df1,aes(x=date, y= p, ymin=lwr, ymax=upr, fill=lab), alpha = 0.2)+
  geom_line(data = test_df2, aes(x=date, y= p, color=lab))+
  geom_ribbon(data = test_df2,aes(x=date, y= p, ymin=lwr, ymax=upr, fill=lab), alpha = 0.2)+
  theme_bw()+
  scale_color_brewer("Test type",palette = "Dark2")+
  scale_fill_brewer("Test type", palette = "Dark2")+
  scale_x_date(date_breaks = "2 month", date_labels =  "%b %Y")+
  ylab("Probability of testing due to symptoms\n(at least one core symptom)")+
  xlab("Date")+
  coord_cartesian(ylim=c(0, 0.4), xlim=c(as.Date("2021-12-15"), as.Date("2023-08-15")))+
  theme(legend.position = c(0.93,0.84),
        legend.background = element_rect(color='black'))


rat <- reformat_data_ratres(sdat, symptom_cat = "at least one core symptom", weeks = week_starts)
rat <- rat[rat$date >= as.Date("2022-02-09")+27,]


test_plt_b <- ggplot()+
  geom_line(data = rat[rat$date<as.Date("2023-06-01"),],aes(x=date, y= p))+
  geom_ribbon(data = rat[rat$date<as.Date("2023-06-01"),],aes(x=date, y= p, ymin=lwr, ymax=upr), alpha = 0.2)+
  theme_bw()+
  scale_x_date(date_breaks = "2 month", date_labels =  "%b\n%Y")+
  scale_y_continuous(breaks = c(0,0.2,0.4,0.6,0.8,1.0))+
  ylab("Probability of reporting\npositive RAT test")+
  xlab("Date")+
  coord_cartesian(ylim=c(0, 1), xlim=c(as.Date("2021-12-15"), as.Date("2023-08-15")))



test_plt_a <- test_plt_a + labs(tag = "A")+
  theme(plot.tag.position = c(0.01, 1),
        axis.text.x = element_blank(),
        axis.title.x = element_blank())

test_plt_b <- test_plt_b + labs(tag = "B")+
  theme(plot.tag.position = c(0.01, 1))

plot_grid(test_plt_a, test_plt_b,  nrow=2, rel_heights = c(1,0.5))

ggsave('figures/Figure4.pdf', width=10, height=8)

test_df1$min_date <- test_df1$date-27
test_df2$min_date <- test_df2$date-27
test_df2[test_df2$min_date < as.Date("2022-01-11"),]$min_date <- as.Date("2022-01-11")
test_df <- rbind(test_df1, test_df2)

write.csv(test_df, 'figures/figure4a.csv')

rat$min_date <- rat$date-27

write.csv(rat, 'figures/figure4b.csv')
#########################################################################################################################################
# Flu tracking vs national
nat1 <- reformat_data(sdat, symptom_cat = "at least fever and cough", weeks = week_starts)


sdat$response2 <- ifelse(sdat$test=="No", 0,
                         ifelse(sdat$test=="Yes – Influenza", 0, 
                                ifelse(sdat$test=="7", NA, 
                                       ifelse(sdat$test == "8", NA, 1))))
nat2 <- reformat_data(sdat, symptom_cat = "at least fever and cough", weeks = week_starts)

sdat$response2 <- ifelse(sdat$test=="No", 0,
                         ifelse(sdat$test_symp=="Yes", 1, 0))


flu <- flu[is.na(flu$RegionCode)==FALSE,]
flu1 <- reformat_data_flu_test(flu, weeks = week_starts)
flu2 <- reformat_data_flu_test(flu[flu$WorkWithPatients==1,], weeks = week_starts)
flu3 <- reformat_data_flu_test(flu[flu$WorkWithPatients!=1,], weeks = week_starts)


nat1$lab <- "National survey (test due to symptoms only)"
nat2$lab <- "National survey (test due to any reason)"
flu1$lab <- "Flutracking (all)"
flu2$lab <- "Flutracking (works with patients)"
flu3$lab <- "Flutracking (does not work with patients)"

nat <- rbind(nat1, nat2)

nat_1 <- nat[nat$date < as.Date("2021-12-20"),]
nat_2 <- nat[nat$date > as.Date("2021-12-20")+28,]


flu_df <- rbind( flu2, flu3)


ggplot()+
  geom_line(data = flu_df,aes(x=date, y= p, color=lab))+
  geom_ribbon(data = flu_df,aes(x=date, y= p, ymin=lwr, ymax=upr, fill=lab), alpha = 0.2)+
  geom_line(data = nat_1,aes(x=date, y= p, color=lab))+
  geom_ribbon(data = nat_1,aes(x=date, y= p, ymin=lwr, ymax=upr, fill=lab), alpha = 0.2)+
  geom_line(data = nat_2,aes(x=date, y= p, color=lab))+
  geom_ribbon(data = nat_2,aes(x=date, y= p, ymin=lwr, ymax=upr, fill=lab), alpha = 0.2)+
  theme_bw()+
  scale_color_brewer("Source of data",palette = "Dark2")+
  scale_fill_brewer("Source of data", palette = "Dark2")+
  scale_x_date(date_breaks = "2 month", date_labels =  "%b\n%Y")+
  ylab("Probability of testing\n(individuals with at least fever and cough)")+
  xlab("Date")+
  coord_cartesian(ylim=c(0, 1), xlim=c(as.Date("2021-12-15"), as.Date("2023-08-15")))+
  theme(legend.position = c(0.8,0.2),
        legend.background = element_rect(color='black'))

ggsave('figures/Figure5.pdf', width=10, height=6)


flu_df$min_date <- flu_df$date-27

nat_1$min_date <- nat_1$date-27
nat_2$min_date <- nat_2$date-27
nat_2[nat_2$min_date < as.Date("2022-01-11"),]$min_date <- as.Date("2022-01-11")
nat_df <- rbind(flu_df, nat_1, nat_2)

write.csv(nat_df, 'figures/figure5.csv')


#########################################################################################################################################
#Flu tracking vs national seek medical care as proxy
nat <- reformat_data(sdat, symptom_cat = "at least fever and cough", weeks = week_starts)

flu1 <- reformat_data_flu_SMC(flu, weeks = week_starts)
flu2 <- reformat_data_flu_SMC(flu[flu$WorkWithPatients==1,], weeks = week_starts)
flu3 <- reformat_data_flu_SMC(flu[flu$WorkWithPatients!=1,], weeks = week_starts)


nat$lab <- "National survey"
#flu1$lab <- "Flutracking (all)"
flu2$lab <- "Flutracking (works with patients)"
flu3$lab <- "Flutracking (does not work with patients)"



nat_1 <- nat[nat$date < as.Date("2021-12-20"),]
nat_2 <- nat[nat$date > as.Date("2021-12-20")+28,]


flu_df <- rbind( flu2, flu3)


ggplot()+
  geom_line(data = flu_df,aes(x=date, y= p, color=lab))+
  geom_ribbon(data = flu_df,aes(x=date, y= p, ymin=lwr, ymax=upr, fill=lab), alpha = 0.2)+
  geom_line(data = nat_1,aes(x=date, y= p, color=lab))+
  geom_ribbon(data = nat_1,aes(x=date, y= p, ymin=lwr, ymax=upr, fill=lab), alpha = 0.2)+
  geom_line(data = nat_2,aes(x=date, y= p, color=lab))+
  geom_ribbon(data = nat_2,aes(x=date, y= p, ymin=lwr, ymax=upr, fill=lab), alpha = 0.2)+
  theme_bw()+
  scale_color_brewer("Test type",palette = "Dark2")+
  scale_fill_brewer("Test type", palette = "Dark2")+
  scale_x_date(date_breaks = "2 month", date_labels =  "%b %Y")+
  ylab("Probability of test given symptoms")+
  xlab("Date")+
  coord_cartesian(ylim=c(0,1),  xlim=c(as.Date("2021-12-01"), as.Date("2023-09-01")))+
  theme(legend.position = c(0.9,0.9),
        legend.background = element_rect(color='black'))



#########################################################################################################################################
## Other reasons to test

rsn1a <- reformat_data_other_reason_to_test(sdat, symptom_cat = "at least one core symptom", weeks = week_starts, cat = "sympt")
rsn2a <- reformat_data_other_reason_to_test(sdat, symptom_cat = "at least one core symptom", weeks = week_starts, cat = "job")
rsn3a <- reformat_data_other_reason_to_test(sdat, symptom_cat = "at least one core symptom", weeks = week_starts, cat = "contact")
rsn4a <- reformat_data_other_reason_to_test(sdat, symptom_cat = "at least one core symptom", weeks = week_starts, cat = "other")

rsn1b <- reformat_data_other_reason_to_test(sdat, symptom_cat = "at least one symptom", weeks = week_starts, cat = "sympt")
rsn2b <- reformat_data_other_reason_to_test(sdat, symptom_cat = "at least one symptom", weeks = week_starts, cat = "job")
rsn3b <- reformat_data_other_reason_to_test(sdat, symptom_cat = "at least one symptom", weeks = week_starts, cat = "contact")
rsn4b <- reformat_data_other_reason_to_test(sdat, symptom_cat = "at least one symptom", weeks = week_starts, cat = "other")

rsn1c <- reformat_data_other_reason_to_test(sdat, symptom_cat = "at least two core symptoms", weeks = week_starts, cat = "sympt")
rsn2c <- reformat_data_other_reason_to_test(sdat, symptom_cat = "at least two core symptoms", weeks = week_starts, cat = "job")
rsn3c <- reformat_data_other_reason_to_test(sdat, symptom_cat = "at least two core symptoms", weeks = week_starts, cat = "contact")
rsn4c <- reformat_data_other_reason_to_test(sdat, symptom_cat = "at least two core symptoms", weeks = week_starts, cat = "other")

rsn1d <- reformat_data_other_reason_to_test(sdat, symptom_cat = "at least fever and cough", weeks = week_starts, cat = "sympt")
rsn2d <- reformat_data_other_reason_to_test(sdat, symptom_cat = "at least fever and cough", weeks = week_starts, cat = "job")
rsn3d <- reformat_data_other_reason_to_test(sdat, symptom_cat = "at least fever and cough", weeks = week_starts, cat = "contact")
rsn4d <- reformat_data_other_reason_to_test(sdat, symptom_cat = "at least fever and cough", weeks = week_starts, cat = "other")


rsn1a$lab <- "Symptoms"
rsn2a$lab <- "Job"
rsn3a$lab <- "Contact"
rsn4a$lab <- "Other"

rsn1b$lab <- "Symptoms"
rsn2b$lab <- "Job"
rsn3b$lab <- "Contact"
rsn4b$lab <- "Other"

rsn1c$lab <- "Symptoms"
rsn2c$lab <- "Job"
rsn3c$lab <- "Contact"
rsn4c$lab <- "Other"

rsn1d$lab <- "Symptoms"
rsn2d$lab <- "Job"
rsn3d$lab <- "Contact"
rsn4d$lab <- "Other"

rsn_dfa <- rbind(rsn1a, rsn2a, rsn3a, rsn4a)
rsn_dfb <- rbind(rsn1b, rsn2b, rsn3b, rsn4b)
rsn_dfc <- rbind(rsn1c, rsn2c, rsn3c, rsn4c)
rsn_dfd <- rbind(rsn1d, rsn2d, rsn3d, rsn4d)

rsn_dfa$symp <- "at least one core symptom"
rsn_dfb$symp <- "at least one symptom"
rsn_dfc$symp <- "at least two core symptoms"
rsn_dfd$symp <- "at least fever and cough"

rsn_df <- rbind(rsn_dfa, rsn_dfb, rsn_dfc, rsn_dfd)
rsn_df1 <- rsn_df[rsn_df$date < as.Date("2021-12-20"),]
rsn_df2 <-rsn_df[rsn_df$date > as.Date("2021-12-20")+28,]


ggplot()+
  geom_line(data = rsn_df1,aes(x=date, y= p, color=lab))+
  geom_ribbon(data = rsn_df1,aes(x=date, y= p, ymin=lwr, ymax=upr, fill=lab), alpha = 0.2)+
  geom_line(data = rsn_df2, aes(x=date, y= p, color=lab))+
  geom_ribbon(data = rsn_df2,aes(x=date, y= p, ymin=lwr, ymax=upr, fill=lab), alpha = 0.2)+
  theme_bw()+
  facet_wrap(.~symp, nrow=2)+
  scale_color_brewer("Reason to test",palette = "Dark2")+
  scale_fill_brewer("Reason to test", palette = "Dark2")+
  scale_x_date(date_breaks = "2 month", date_labels =  "%b\n%Y")+
  ylab("Probability of testing")+
  xlab("Date")+
  coord_cartesian(ylim=c(0, 1), xlim=c(as.Date("2021-12-15"), as.Date("2023-08-15")))+
  theme(legend.position = c(0.9,0.9),
        legend.background = element_rect(color='black'))

ggsave('figures/SFigure1.pdf', width=10, height=8)


rsn_df1$min_date <- rsn_df1$date-27
rsn_df2$min_date <- rsn_df2$date-27
rsn_df2[rsn_df2$min_date < as.Date("2022-01-11"),]$min_date <- as.Date("2022-01-11")
rsn_df <- rbind(rsn_df1, rsn_df2)

write.csv(rsn_df, 'figures/Sfigure1.csv')




######################################################################################################
# Supplementary Tables

rtt1 <- get_tab_rtt(sdat, sympt_cat = "at least one symptom")
rtt2 <- get_tab_rtt(sdat, "at least one core symptom")
rtt3 <- get_tab_rtt(sdat, "at least fever and cough")
rtt4 <- get_tab_rtt(sdat, "at least two core symptoms")
rtt1$symp <- "at least one symptom"
rtt2$symp <- "at least one core symptom"
rtt3$symp <- "at least fever and cough"
rtt4$symp <- "at least two core symptoms"

rtt_df <- rbind(rtt1, rtt2, rtt3, rtt4)
write.csv(rtt_df, 'figures/rtt_sup.csv')


########################################################################################################
st1 <- get_tab_rtt(sdat[sdat$state %in% c("NSW","VIC","ACT"),], "at least one core symptom")
st2 <- get_tab_rtt(sdat[sdat$state %in% c("QLD","SA","TAS","NT"),], "at least one core symptom")
st3 <- get_tab_rtt(sdat[sdat$state %in% c("WA"),], "at least one core symptom")

st1$lab2 <- "NSW/VIC/ACT"
st2$lab2 <- "QLD/SA/TAS/NT"
st3$lab2 <- "WA"

ag1 <- get_tab_rtt(sdat[sdat$age < 30,], "at least one core symptom")
ag2 <- get_tab_rtt(sdat[sdat$age >= 30 & sdat$age < 60,], "at least one core symptom")
ag3 <- get_tab_rtt(sdat[sdat$age >= 60,], "at least one core symptom")

ag1$lab2 <- "18-29"
ag2$lab2 <- "30-59"
ag3$lab2 <- "60+"

st_ag_df <- rbind(st1, st2, st3,
                  ag1, ag2, ag3)

write.csv(st_ag_df, 'figures/state_age_average_sup.csv')

###################################################################################################
# Supplementary table or figure
# By state
st1 <- get_tab_rtt(sdat[sdat$state %in% c("NSW"),], "at least one core symptom")
st2 <- get_tab_rtt(sdat[sdat$state %in% c("VIC"),], "at least one core symptom")
st3 <- get_tab_rtt(sdat[sdat$state %in% c("ACT"),], "at least one core symptom")
st4 <- get_tab_rtt(sdat[sdat$state %in% c("QLD"),], "at least one core symptom")
st5 <- get_tab_rtt(sdat[sdat$state %in% c("SA"),], "at least one core symptom")
st6 <- get_tab_rtt(sdat[sdat$state %in% c("NT"),], "at least one core symptom")
st7 <- get_tab_rtt(sdat[sdat$state %in% c("TAS"),], "at least one core symptom")
st8 <- get_tab_rtt(sdat[sdat$state %in% c("WA"),], "at least one core symptom")

st1$lab2 <- "NSW"
st2$lab2 <- "VIC"
st3$lab2 <- "ACT"
st4$lab2 <- "QLD"
st5$lab2 <- "SA"
st6$lab2 <- "NT"
st7$lab2 <- "TAS"
st8$lab2 <- "WA"

st_sup_df <- rbind(st1, st2, st3, st4,
                   st5, st6, st7, st8)
df <- st_sup_df[st_sup_df$lab =="Sympt",]
df$pvalue <- "ref"

for(i in 2:nrow(df)){
  row1 <- data.frame(p = df$x[1],
                     n = df$n[1]-df$x[1])
  row2 <- data.frame(p = df$x[i],
                     n = df$n[i]-df$x[i])
  tdf <- rbind(row1, row2)
  
  df$pvalue[i] <- fisher.test(tdf)$p.value
  
}

df_state_all <- df


#By age
st1 <- get_tab_rtt(sdat[sdat$age<25,], "at least one core symptom")
st2 <- get_tab_rtt(sdat[sdat$age>=25 & sdat$age <30,], "at least one core symptom")
st3 <- get_tab_rtt(sdat[sdat$age>=30 & sdat$age <35,], "at least one core symptom")
st4 <- get_tab_rtt(sdat[sdat$age>=35 & sdat$age <40,], "at least one core symptom")
st5 <- get_tab_rtt(sdat[sdat$age>=40 & sdat$age <45,], "at least one core symptom")
st6 <- get_tab_rtt(sdat[sdat$age>=45 & sdat$age <50,], "at least one core symptom")
st7 <- get_tab_rtt(sdat[sdat$age>=50 & sdat$age <55,], "at least one core symptom")
st8 <- get_tab_rtt(sdat[sdat$age>=55 & sdat$age <60,], "at least one core symptom")
st9 <- get_tab_rtt(sdat[sdat$age>=60 & sdat$age <65,], "at least one core symptom")
st10 <- get_tab_rtt(sdat[sdat$age>=65,], "at least one core symptom")

st1$lab2 <- "18-24"
st2$lab2 <- "25-29"
st3$lab2 <- "30-34"
st4$lab2 <- "35-39"
st5$lab2 <- "40-44"
st6$lab2 <- "45-49"
st7$lab2 <- "50-54"
st8$lab2 <- "55-59"
st9$lab2 <- "60-64"
st10$lab2 <- "65+"



st_sup_df <- rbind(st1, st2, st3, st4,
                   st5, st6, st7, st8,
                   st9, st10)
df <- st_sup_df[st_sup_df$lab =="Sympt",]
df$pvalue <- "ref"

for(i in 2:nrow(df)){
  row1 <- data.frame(p = df$x[1],
                     n = df$n[1]-df$x[1])
  row2 <- data.frame(p = df$x[i],
                     n = df$n[i]-df$x[i])
  tdf <- rbind(row1, row2)
  
  df$pvalue[i] <- fisher.test(tdf)$p.value
  
}

df_age_all <- df

write.csv(rbind(df_state_all, df_age_all), "figures/STab2.csv")


new1<-ggplot(df_age_all)+
  geom_point(aes(x=proportion, y=lab2))+
  geom_errorbar(aes(x=proportion, y=lab2, xmin=lower, xmax=upper), width=0.)+
  theme_bw()+
  coord_cartesian(xlim=c(0.28, 0.45))+
  ylab("Age-group")+
  xlab("Probability of testing due to symptoms\n(at least one core symptom)")

new2<-ggplot(df_state_all)+
  geom_point(aes(x=proportion, y=lab2))+
  geom_errorbar(aes(x=proportion, y=lab2, xmin=lower, xmax=upper), width=0.)+
  theme_bw()+
  coord_cartesian(xlim=c(0.28, 0.45))+
  ylab("State/territory")+
  xlab("Probability of testing due to symptoms\n(at least one core symptom)")


new1 <- new1+labs(tag="D")+
  theme(plot.tag.position = c(0.01,0.99),
        panel.grid = element_blank())

new2 <- new2+labs(tag="C")+
  theme(plot.tag.position = c(0.01,0.99),
        axis.title.x = element_blank(),
        panel.grid  = element_blank())

plot_grid(new2, new1, nrow=2, rel_heights = c(1,1))



plt_as_a <- plt_as_a+labs(tag="A")+
  theme(plot.tag.position = c(0.02, 0.98),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        legend.position = c(0.85,0.8))

plt_as_b <- plt_as_b+labs(tag="B")+
  theme(plot.tag.position = c(0.02, 0.98))

grd1<-plot_grid(plt_as_a, plt_as_b, nrow=2, rel_heights = c(1,1.1))
grd2<-plot_grid(new2, new1, nrow=2, rel_heights = c(1,1))

plot_grid(grd1, grd2, nrow=1, rel_widths = c(1,0.5))
ggsave('figures/Figure3NEW.pdf', width=10, height=8)


sdat$response2 <- ifelse(sdat$test=="No", 0,
                         ifelse(sdat$test_symp=="Yes" & sdat$what_test == "PCR only", 1, 0))
st1 <- get_tab_rtt_RATPCR(sdat[sdat$date>=as.Date("2022-01-11"),], "at least one core symptom")

sdat$response2 <- ifelse(sdat$test=="No", 0,
                         ifelse(sdat$test_symp=="Yes" & sdat$what_test == "RAT only", 1, 0))
st2 <- get_tab_rtt_RATPCR(sdat[sdat$date>=as.Date("2022-01-11"),], "at least one core symptom")

sdat$response2 <- ifelse(sdat$test=="No", 0,
                         ifelse(sdat$test_symp=="Yes" & sdat$what_test == "Both", 1, 0))
st3 <- get_tab_rtt_RATPCR(sdat[sdat$date>=as.Date("2022-01-11"),], "at least one core symptom")

sdat$response2 <- ifelse(sdat$test=="No", 0,
                         ifelse(sdat$test_symp=="Yes", 1, 0))


st1$lab2 <- "PCR only"
st2$lab2 <- "RAT only"
st3$lab2 <- "Both"

st_sup_df <- rbind(st1, st2, st3)
df <- st_sup_df[st_sup_df$lab =="Sympt",]
df$pvalue <- "ref"

for(i in 2:nrow(df)){
  row1 <- data.frame(p = df$x[1],
                     n = df$n[1]-df$x[1])
  row2 <- data.frame(p = df$x[i],
                     n = df$n[i]-df$x[i])
  tdf <- rbind(row1, row2)
  
  df$pvalue[i] <- fisher.test(tdf)$p.value
  
}
df


tab <- table(sdat$rat_report)
tab <- table(sdat[sdat$date<as.Date("2023-06-01"),]$rat_report)


tmp<-binom.exact(tab[[2]],tab[[1]]+ tab[[2]])
tmp$lab <- "Rat report"
tmp$lab2 <- "Rat report"
tmp$pvalue <- "Rat report"
write.csv(rbind(df,tmp),"figures/What_test_sup.csv")


#################################################################################################################################################
################################################################################
rat_pos1 <- reformat_data_rat_pos(sdat, symptom_cat = "at least one core symptom", weeks = week_starts)
rat_pos2 <- reformat_data_rat_pos(sdat, symptom_cat = "at least one symptom", weeks = week_starts)
rat_pos3 <- reformat_data_rat_pos(sdat, symptom_cat = "at least two core symptoms", weeks = week_starts)
rat_pos4 <- reformat_data_rat_pos(sdat, symptom_cat = "at least fever and cough", weeks = week_starts)

rat_pos1$lab <- "at least one core symptom"
rat_pos2$lab <- "at least one symptom"
rat_pos3$lab <- "at least two core symptoms"
rat_pos4$lab <- "at least fever and cough"

rat_pos_df <- rbind(rat_pos1, rat_pos2, rat_pos3, rat_pos4)


rat_pos_df$lab <- factor(rat_pos_df$lab, levels=c("at least one symptom", "at least one core symptom", "at least two core symptoms", "at least fever and cough"))

rat_pos_df1 <- rat_pos_df[rat_pos_df$date < as.Date("2021-12-20"),]
rat_pos_df2 <- rat_pos_df[rat_pos_df$date >= as.Date("2022-01-11")+27,]


plt_pos_a<-ggplot()+
  geom_line(data = rat_pos_df1,aes(x=date, y= p, color=lab))+
  geom_ribbon(data = rat_pos_df1,aes(x=date, y= p, ymin=lwr, ymax=upr, fill=lab), alpha = 0.2)+
  geom_line(data = rat_pos_df2, aes(x=date, y= p, color=lab))+
  geom_ribbon(data = rat_pos_df2,aes(x=date, y= p, ymin=lwr, ymax=upr, fill=lab), alpha = 0.2)+
  theme_bw()+
  scale_color_brewer("Symptom status",palette = "Dark2")+
  scale_fill_brewer("Symptom status", palette = "Dark2")+
  scale_x_date(date_breaks = "2 month", date_labels =  "%b\n%Y")+
  ylab("Probability of positive\nRAT test given symptoms")+
  coord_cartesian(ylim=c(0, 1), xlim=c(as.Date("2022-02-20"), as.Date("2023-08-15")))+
  theme(legend.position = c(0.85,0.8),
        legend.background = element_rect(color='black'))


#################################################################################################################################################
################################################################################
pcr_pos1 <- reformat_data_pcr_pos(sdat, symptom_cat = "at least one core symptom", weeks = week_starts)
pcr_pos2 <- reformat_data_pcr_pos(sdat, symptom_cat = "at least one symptom", weeks = week_starts)
pcr_pos3 <- reformat_data_pcr_pos(sdat, symptom_cat = "at least two core symptoms", weeks = week_starts)
pcr_pos4 <- reformat_data_pcr_pos(sdat, symptom_cat = "at least fever and cough", weeks = week_starts)

pcr_pos1$lab <- "at least one core symptom"
pcr_pos2$lab <- "at least one symptom"
pcr_pos3$lab <- "at least two core symptoms"
pcr_pos4$lab <- "at least fever and cough"

pcr_pos_df <- rbind(pcr_pos1, pcr_pos2, pcr_pos3, pcr_pos4)

pcr_pos_df$lab <- factor(pcr_pos_df$lab, levels=c("at least one symptom", "at least one core symptom", "at least two core symptoms", "at least fever and cough"))

pcr_pos_df1 <- pcr_pos_df[pcr_pos_df$date < as.Date("2021-12-20"),]
pcr_pos_df2 <- pcr_pos_df[pcr_pos_df$date >= as.Date("2022-01-11")+27,]


plt_pos_b<-ggplot()+
  geom_line(data = pcr_pos_df1,aes(x=date, y= p, color=lab))+
  geom_ribbon(data = pcr_pos_df1,aes(x=date, y= p, ymin=lwr, ymax=upr, fill=lab), alpha = 0.2)+
  geom_line(data = pcr_pos_df2, aes(x=date, y= p, color=lab))+
  geom_ribbon(data = pcr_pos_df2,aes(x=date, y= p, ymin=lwr, ymax=upr, fill=lab), alpha = 0.2)+
  theme_bw()+
  scale_color_brewer("Symptom status",palette = "Dark2")+
  scale_fill_brewer("Symptom status", palette = "Dark2")+
  scale_x_date(date_breaks = "2 month", date_labels =  "%b\n%Y")+
  ylab("Probability of positive\nPCR test given symptoms")+
  xlab("Date")+
  coord_cartesian(ylim=c(0, 1), xlim=c(as.Date("2022-02-20"), as.Date("2023-08-15")))+
  theme(legend.position = c(0.9,0.9),
        legend.background = element_rect(color='black'))



plt_pos_a <- plt_pos_a + labs(tag="A")+
  theme(plot.tag.position = c(0.02,0.98),
        axis.title.x = element_blank(),
        axis.text.x = element_blank())

plt_pos_b <- plt_pos_b + labs(tag="B")+
  theme(plot.tag.position = c(0.02,0.98),
        legend.position="none")
plot_grid(plt_pos_a, plt_pos_b, nrow=2, rel_heights=c(1,1.1))

ggsave('figures/SFigure3.pdf', width=8, height=8.5)


rat_pos_df1$min_date <- rat_pos_df1$date-27
rat_pos_df2$min_date <- rat_pos_df2$date-27
rat_pos_df2[rat_pos_df2$min_date < as.Date("2022-01-11"),]$min_date <- as.Date("2022-01-11")
rat_pos_df <- rbind(rat_pos_df1, rat_pos_df2)

write.csv(rat_pos_df, 'figures/Sfigure3a.csv')

pcr_pos_df1$min_date <- pcr_pos_df1$date-27
pcr_pos_df2$min_date <- pcr_pos_df2$date-27
pcr_pos_df2[pcr_pos_df2$min_date < as.Date("2022-01-11"),]$min_date <- as.Date("2022-01-11")
pcr_pos_df <- rbind(pcr_pos_df1, pcr_pos_df2)

write.csv(pcr_pos_df, 'figures/Sfigure3b.csv')

ST41<-pcr_pos_df1[pcr_pos_df1$week=="Overall",]
ST42<-rat_pos_df1[rat_pos_df1$week=="Overall",]

write.csv(rbind(ST41, ST42), 'figures/Stab4.csv')


tdat <- sdat[sdat$symptoms == sympt_cat & sdat$response=="Yes",]

#####################################################################################


tdat <- flu[flu$IsCase == 1 & flu$WorkWithPatients==1,]

test <- table(tdat$SurveyWeek, tdat$TestCovid)
tab_df <- data.frame(neg = as.numeric(test[,1]),
                     pos= as.numeric(test[,2]),
                     date = rownames(test))


neg <- sum(tab_df$neg)
pos <- sum(tab_df$pos)

binom.exact(pos, pos+neg)


tdat <- flu[flu$IsCase == 1 & flu$WorkWithPatients!=1,]

test <- table(tdat$SurveyWeek, tdat$TestCovid)
tab_df <- data.frame(neg = as.numeric(test[,1]),
                     pos= as.numeric(test[,2]),
                     date = rownames(test))


neg <- sum(tab_df$neg)
pos <- sum(tab_df$pos)

binom.exact(pos, pos+neg)

###############################################################################################################################
flu <- mutate(flu, age_group = ifelse(Age<5, "0-4",
                                      ifelse(Age<5, "0-4",
                                             ifelse(Age<10, "5-9",
                                                    ifelse(Age<15, "10-14",
                                                           ifelse(Age<20, "15-19",
                                                                  ifelse(Age<25, "20-24",
                                                                         ifelse(Age<30, "25-39",
                                                                                ifelse(Age<35, "30-34",
                                                                                       ifelse(Age<40, "35-39",
                                                                                              ifelse(Age<45, "40-44",
                                                                                                     ifelse(Age<50, "45-49",
                                                                                                            ifelse(Age<55, "50-54",
                                                                                                                   ifelse(Age<60, "55-59",
                                                                                                                          ifelse(Age<65, "60-64",
                                                                                                                                 ifelse(Age<70, "65-69",
                                                                                                                                        ifelse(Age<75, "70-74","75+")))))))))))))))))

sdat <- mutate(sdat, age_group = ifelse(age<5, "0-4",
                                      ifelse(age<5, "0-4",
                                             ifelse(age<10, "5-9",
                                                    ifelse(age<15, "10-14",
                                                           ifelse(age<20, "15-19",
                                                                  ifelse(age<25, "20-24",
                                                                         ifelse(age<30, "25-39",
                                                                                ifelse(age<35, "30-34",
                                                                                       ifelse(age<40, "35-39",
                                                                                              ifelse(age<45, "40-44",
                                                                                                     ifelse(age<50, "45-49",
                                                                                                            ifelse(age<55, "50-54",
                                                                                                                   ifelse(age<60, "55-59",
                                                                                                                          ifelse(age<65, "60-64",
                                                                                                                                 ifelse(age<70, "65-69",
                                                                                                                                        ifelse(age<75, "70-74","75+")))))))))))))))))


tabF_R<- table(flu$RegionCode, exclude=NULL)
tabF_A<- table(flu$age_group, exclude=NULL)
tabF_G<- table(flu$Gender, exclude=NULL)

tabS_R<- table(sdat[sdat$symptoms%in%c("at least one core symptom"),]$state)
tabS_A<- table(sdat[sdat$symptoms%in%c("at least one core symptom"),]$age_group)
tabS_G<- table(sdat[sdat$symptoms%in%c("at least one core symptom"),]$gender)


dfF_R <- data.frame( per=100*tabF_R/sum(tabF_R),
                     num = tabF_R)
dfF_A <- data.frame( per=100*tabF_A/sum(tabF_A),
                     num = tabF_A)
dfF_G <- data.frame( per=100*tabF_G/sum(tabF_G),
                     num = tabF_G)

dfS_R <- data.frame( per=100*tabS_R/sum(tabS_R),
                     num = tabS_R)
dfS_A <- data.frame( per=100*tabS_A/sum(tabS_A),
                     num = tabS_A)
dfS_G <- data.frame( per=100*tabS_G/sum(tabS_G),
                     num = tabS_G)


df <- rbind(dfF_R,
            dfF_A,
            dfF_G,
            dfS_R,
            dfS_A,
            dfS_G)

write.csv(df, 'figures/representativeness_tab.csv')



####################################################################################################################
# Alternate analysis
date_brk1 <- as.Date("2022-06-20")
date_brk2 <- as.Date("2022-11-07")
###################################################################################################
### By state
# First time period

st1 <- get_tab_rtt(sdat[sdat$state %in% c("NSW") & sdat$date <= date_brk1 ,], "at least one core symptom")
st2 <- get_tab_rtt(sdat[sdat$state %in% c("VIC") & sdat$date <= date_brk1 ,], "at least one core symptom")
st3 <- get_tab_rtt(sdat[sdat$state %in% c("ACT") & sdat$date <= date_brk1 ,], "at least one core symptom")
st4 <- get_tab_rtt(sdat[sdat$state %in% c("QLD") & sdat$date <= date_brk1 ,], "at least one core symptom")
st5 <- get_tab_rtt(sdat[sdat$state %in% c("SA") & sdat$date <= date_brk1 ,], "at least one core symptom")
st6 <- get_tab_rtt(sdat[sdat$state %in% c("NT") & sdat$date <= date_brk1 ,], "at least one core symptom")
st7 <- get_tab_rtt(sdat[sdat$state %in% c("TAS") & sdat$date <= date_brk1 ,], "at least one core symptom")
st8 <- get_tab_rtt(sdat[sdat$state %in% c("WA") & sdat$date <= date_brk1 ,], "at least one core symptom")

st1$lab2 <- "NSW"
st2$lab2 <- "VIC"
st3$lab2 <- "ACT"
st4$lab2 <- "QLD"
st5$lab2 <- "SA"
st6$lab2 <- "NT"
st7$lab2 <- "TAS"
st8$lab2 <- "WA"

st_sup_df <- rbind(st1, st2, st3, st4,
                   st5, st6, st7, st8)
df <- st_sup_df[st_sup_df$lab =="Sympt",]
df$pvalue <- "ref"

for(i in 2:nrow(df)){
  row1 <- data.frame(p = df$x[1],
                     n = df$n[1]-df$x[1])
  row2 <- data.frame(p = df$x[i],
                     n = df$n[i]-df$x[i])
  tdf <- rbind(row1, row2)
  
  df$pvalue[i] <- fisher.test(tdf)$p.value
  
}

df_state_all_a <- df


# Second time period
st1 <- get_tab_rtt(sdat[sdat$state %in% c("NSW") & sdat$date <= date_brk2 & sdat$date > date_brk1,], "at least one core symptom")
st2 <- get_tab_rtt(sdat[sdat$state %in% c("VIC") & sdat$date <= date_brk2 & sdat$date > date_brk1,], "at least one core symptom")
st3 <- get_tab_rtt(sdat[sdat$state %in% c("ACT") & sdat$date <= date_brk2 & sdat$date > date_brk1,], "at least one core symptom")
st4 <- get_tab_rtt(sdat[sdat$state %in% c("QLD") & sdat$date <= date_brk2 & sdat$date > date_brk1,], "at least one core symptom")
st5 <- get_tab_rtt(sdat[sdat$state %in% c("SA") & sdat$date <= date_brk2 & sdat$date > date_brk1,], "at least one core symptom")
st6 <- get_tab_rtt(sdat[sdat$state %in% c("NT") & sdat$date <= date_brk2 & sdat$date > date_brk1,], "at least one core symptom")
st7 <- get_tab_rtt(sdat[sdat$state %in% c("TAS") & sdat$date <= date_brk2 & sdat$date > date_brk1,], "at least one core symptom")
st8 <- get_tab_rtt(sdat[sdat$state %in% c("WA") & sdat$date <= date_brk2 & sdat$date > date_brk1,], "at least one core symptom")

st1$lab2 <- "NSW"
st2$lab2 <- "VIC"
st3$lab2 <- "ACT"
st4$lab2 <- "QLD"
st5$lab2 <- "SA"
st6$lab2 <- "NT"
st7$lab2 <- "TAS"
st8$lab2 <- "WA"

st_sup_df <- rbind(st1, st2, st3, st4,
                   st5, st6, st7, st8)
df <- st_sup_df[st_sup_df$lab =="Sympt",]
df$pvalue <- "ref"

for(i in 2:nrow(df)){
  row1 <- data.frame(p = df$x[1],
                     n = df$n[1]-df$x[1])
  row2 <- data.frame(p = df$x[i],
                     n = df$n[i]-df$x[i])
  tdf <- rbind(row1, row2)
  
  df$pvalue[i] <- fisher.test(tdf)$p.value
  
}

df_state_all_b <- df

# Third time period
st1 <- get_tab_rtt(sdat[sdat$state %in% c("NSW") & sdat$date > date_brk2,], "at least one core symptom")
st2 <- get_tab_rtt(sdat[sdat$state %in% c("VIC") & sdat$date > date_brk2,], "at least one core symptom")
st3 <- get_tab_rtt(sdat[sdat$state %in% c("ACT") & sdat$date > date_brk2,], "at least one core symptom")
st4 <- get_tab_rtt(sdat[sdat$state %in% c("QLD") & sdat$date > date_brk2,], "at least one core symptom")
st5 <- get_tab_rtt(sdat[sdat$state %in% c("SA") & sdat$date > date_brk2,], "at least one core symptom")
st6 <- get_tab_rtt(sdat[sdat$state %in% c("NT") & sdat$date > date_brk2,], "at least one core symptom")
st7 <- get_tab_rtt(sdat[sdat$state %in% c("TAS") & sdat$date > date_brk2,], "at least one core symptom")
st8 <- get_tab_rtt(sdat[sdat$state %in% c("WA") & sdat$date > date_brk2,], "at least one core symptom")

st1$lab2 <- "NSW"
st2$lab2 <- "VIC"
st3$lab2 <- "ACT"
st4$lab2 <- "QLD"
st5$lab2 <- "SA"
st6$lab2 <- "NT"
st7$lab2 <- "TAS"
st8$lab2 <- "WA"

st_sup_df <- rbind(st1, st2, st3, st4,
                   st5, st6, st7, st8)
df <- st_sup_df[st_sup_df$lab =="Sympt",]
df$pvalue <- "ref"

for(i in 2:nrow(df)){
  row1 <- data.frame(p = df$x[1],
                     n = df$n[1]-df$x[1])
  row2 <- data.frame(p = df$x[i],
                     n = df$n[i]-df$x[i])
  tdf <- rbind(row1, row2)
  
  df$pvalue[i] <- fisher.test(tdf)$p.value
  
}

df_state_all_c<- df


### By age
# First time period
st1 <- get_tab_rtt(sdat[sdat$age<25 & sdat$date <= date_brk1 ,], "at least one core symptom")
st2 <- get_tab_rtt(sdat[sdat$age>=25 & sdat$age <30 & sdat$date <= date_brk1 ,], "at least one core symptom")
st3 <- get_tab_rtt(sdat[sdat$age>=30 & sdat$age <35 & sdat$date <= date_brk1 ,], "at least one core symptom")
st4 <- get_tab_rtt(sdat[sdat$age>=35 & sdat$age <40 & sdat$date <= date_brk1 ,], "at least one core symptom")
st5 <- get_tab_rtt(sdat[sdat$age>=40 & sdat$age <45 & sdat$date <= date_brk1 ,], "at least one core symptom")
st6 <- get_tab_rtt(sdat[sdat$age>=45 & sdat$age <50 & sdat$date <= date_brk1 ,], "at least one core symptom")
st7 <- get_tab_rtt(sdat[sdat$age>=50 & sdat$age <55 & sdat$date <= date_brk1 ,], "at least one core symptom")
st8 <- get_tab_rtt(sdat[sdat$age>=55 & sdat$age <60 & sdat$date <= date_brk1 ,], "at least one core symptom")
st9 <- get_tab_rtt(sdat[sdat$age>=60 & sdat$age <65 & sdat$date <= date_brk1 ,], "at least one core symptom")
st10 <- get_tab_rtt(sdat[sdat$age>=65 & sdat$date <= date_brk1,], "at least one core symptom")

st1$lab2 <- "18-24"
st2$lab2 <- "25-29"
st3$lab2 <- "30-34"
st4$lab2 <- "35-39"
st5$lab2 <- "40-44"
st6$lab2 <- "45-49"
st7$lab2 <- "50-54"
st8$lab2 <- "55-59"
st9$lab2 <- "60-64"
st10$lab2 <- "65+"

st_sup_df <- rbind(st1, st2, st3, st4,
                   st5, st6, st7, st8,
                   st9, st10)
df <- st_sup_df[st_sup_df$lab =="Sympt",]
df$pvalue <- "ref"

for(i in 2:nrow(df)){
  row1 <- data.frame(p = df$x[1],
                     n = df$n[1]-df$x[1])
  row2 <- data.frame(p = df$x[i],
                     n = df$n[i]-df$x[i])
  tdf <- rbind(row1, row2)
  
  df$pvalue[i] <- fisher.test(tdf)$p.value
  
}

df_age_all_a <- df


# Second time period
st1 <- get_tab_rtt(sdat[sdat$age<25 & sdat$date <= date_brk2 & sdat$date > date_brk1,], "at least one core symptom")
st2 <- get_tab_rtt(sdat[sdat$age>=25 & sdat$age <30 & sdat$date <= date_brk2 & sdat$date > date_brk1,], "at least one core symptom")
st3 <- get_tab_rtt(sdat[sdat$age>=30 & sdat$age <35 & sdat$date <= date_brk2 & sdat$date > date_brk1,], "at least one core symptom")
st4 <- get_tab_rtt(sdat[sdat$age>=35 & sdat$age <40 & sdat$date <= date_brk2 & sdat$date > date_brk1,], "at least one core symptom")
st5 <- get_tab_rtt(sdat[sdat$age>=40 & sdat$age <45 & sdat$date <= date_brk2 & sdat$date > date_brk1,], "at least one core symptom")
st6 <- get_tab_rtt(sdat[sdat$age>=45 & sdat$age <50 & sdat$date <= date_brk2 & sdat$date > date_brk1,], "at least one core symptom")
st7 <- get_tab_rtt(sdat[sdat$age>=50 & sdat$age <55 & sdat$date <= date_brk2 & sdat$date > date_brk1,], "at least one core symptom")
st8 <- get_tab_rtt(sdat[sdat$age>=55 & sdat$age <60 & sdat$date <= date_brk2 & sdat$date > date_brk1,], "at least one core symptom")
st9 <- get_tab_rtt(sdat[sdat$age>=60 & sdat$age <65 & sdat$date <= date_brk2 & sdat$date > date_brk1,], "at least one core symptom")
st10 <- get_tab_rtt(sdat[sdat$age>=65 & sdat$date <= date_brk2 & sdat$date > date_brk1,], "at least one core symptom")

st1$lab2 <- "18-24"
st2$lab2 <- "25-29"
st3$lab2 <- "30-34"
st4$lab2 <- "35-39"
st5$lab2 <- "40-44"
st6$lab2 <- "45-49"
st7$lab2 <- "50-54"
st8$lab2 <- "55-59"
st9$lab2 <- "60-64"
st10$lab2 <- "65+"



st_sup_df <- rbind(st1, st2, st3, st4,
                   st5, st6, st7, st8,
                   st9, st10)
df <- st_sup_df[st_sup_df$lab =="Sympt",]
df$pvalue <- "ref"

for(i in 2:nrow(df)){
  row1 <- data.frame(p = df$x[1],
                     n = df$n[1]-df$x[1])
  row2 <- data.frame(p = df$x[i],
                     n = df$n[i]-df$x[i])
  tdf <- rbind(row1, row2)
  
  df$pvalue[i] <- fisher.test(tdf)$p.value
  
}

df_age_all_b <- df

# Third time period
st1 <- get_tab_rtt(sdat[sdat$age<25 & sdat$date > date_brk2,], "at least one core symptom")
st2 <- get_tab_rtt(sdat[sdat$age>=25 & sdat$age <30 & sdat$date > date_brk2,], "at least one core symptom")
st3 <- get_tab_rtt(sdat[sdat$age>=30 & sdat$age <35 & sdat$date > date_brk2,], "at least one core symptom")
st4 <- get_tab_rtt(sdat[sdat$age>=35 & sdat$age <40 & sdat$date > date_brk2,], "at least one core symptom")
st5 <- get_tab_rtt(sdat[sdat$age>=40 & sdat$age <45 & sdat$date > date_brk2,], "at least one core symptom")
st6 <- get_tab_rtt(sdat[sdat$age>=45 & sdat$age <50 & sdat$date > date_brk2,], "at least one core symptom")
st7 <- get_tab_rtt(sdat[sdat$age>=50 & sdat$age <55 & sdat$date > date_brk2,], "at least one core symptom")
st8 <- get_tab_rtt(sdat[sdat$age>=55 & sdat$age <60 & sdat$date > date_brk2,], "at least one core symptom")
st9 <- get_tab_rtt(sdat[sdat$age>=60 & sdat$age <65 & sdat$date > date_brk2,], "at least one core symptom")
st10 <- get_tab_rtt(sdat[sdat$age>=65 & sdat$date > date_brk2,], "at least one core symptom")

st1$lab2 <- "18-24"
st2$lab2 <- "25-29"
st3$lab2 <- "30-34"
st4$lab2 <- "35-39"
st5$lab2 <- "40-44"
st6$lab2 <- "45-49"
st7$lab2 <- "50-54"
st8$lab2 <- "55-59"
st9$lab2 <- "60-64"
st10$lab2 <- "65+"


st_sup_df <- rbind(st1, st2, st3, st4,
                   st5, st6, st7, st8,
                   st9, st10)
df <- st_sup_df[st_sup_df$lab =="Sympt",]
df$pvalue <- "ref"

for(i in 2:nrow(df)){
  row1 <- data.frame(p = df$x[1],
                     n = df$n[1]-df$x[1])
  row2 <- data.frame(p = df$x[i],
                     n = df$n[i]-df$x[i])
  tdf <- rbind(row1, row2)
  
  df$pvalue[i] <- fisher.test(tdf)$p.value
  
}

df_age_all_c <- df


write.csv(rbind(df_state_all_a, df_state_all_b, df_state_all_c, 
                df_age_all_a, df_age_all_b, df_age_all_c), "figures/STab2_alt.csv")


df_age_all_a$facet_lab <- "BA.1/BA.2"
df_age_all_b$facet_lab <- "BA.4/BA.5"
df_age_all_c$facet_lab <- "Assorted Omicron variants"
df_age_all_comb <- rbind(df_age_all_a,
                         df_age_all_b,
                         df_age_all_c)

df_state_all_a$facet_lab <- "BA.1/BA.2"
df_state_all_b$facet_lab <- "BA.4/BA.5"
df_state_all_c$facet_lab <- "Assorted Omicron variants"
df_state_all_comb <- rbind(df_state_all_a,
                           df_state_all_b,
                           df_state_all_c)

df_state_all_comb$facet_lab <- factor(df_state_all_comb$facet_lab, levels = c("BA.1/BA.2","BA.4/BA.5","Assorted Omicron variants"))
df_age_all_comb$facet_lab <- factor(df_age_all_comb$facet_lab, levels = c("BA.1/BA.2","BA.4/BA.5","Assorted Omicron variants"))


new1<-ggplot(df_age_all_comb)+
  geom_point(aes(x=proportion, y=lab2, color=facet_lab), position = position_dodge(.5))+
  geom_errorbar(aes(x=proportion, y=lab2, xmin=lower, xmax=upper, color=facet_lab), width=0.,  position = position_dodge(.5))+
  theme_bw()+
  scale_color_brewer("Time period",palette = "Set1")+
  ylab("Age-group")+
  xlab("Probability of testing due to symptoms\n(at least one core symptom)")+
  theme(legend.background = element_rect(color = "black"),
        legend.position = "bottom")+
  guides(color=guide_legend(ncol=1,byrow=TRUE, title.position = "top"))

new2<-ggplot(df_state_all_comb)+
  geom_point(aes(x=proportion, y=lab2, color=facet_lab), position = position_dodge(.5))+
  geom_errorbar(aes(x=proportion, y=lab2, xmin=lower, xmax=upper, color=facet_lab), width=0.,  position = position_dodge(.5))+
  theme_bw()+
  scale_color_brewer("Time period",palette = "Set1")+
  ylab("State/territory")+
  xlab("Probability of testing due to symptoms\n(at least one core symptom)")+
  theme(legend.background = element_rect(color = "black"),
        legend.position = "bottom")+
  guides(color=guide_legend(ncol=1,byrow=TRUE, title.position = "top"))


new1 <- new1+labs(tag="B")+
  theme(plot.tag.position = c(0.01,0.99),
        panel.grid = element_blank())

new2 <- new2+labs(tag="B")+
  theme(plot.tag.position = c(0.01,0.99),
        panel.grid  = element_blank())

plt_as_b <- plt_as_b+labs(tag="A")+
  geom_vline(xintercept = as.Date("2022-06-20"), linetype="dashed")+
  geom_vline(xintercept = as.Date("2022-11-07"), linetype="dashed")

plt_as_a <- plt_as_a2+labs(tag="A")+
  theme(plot.tag.position = c(0.02, 0.98),
        legend.position = c(0.85,0.8))+
  geom_vline(xintercept = as.Date("2022-06-20"), linetype="dashed")+
  geom_vline(xintercept = as.Date("2022-11-07"), linetype="dashed")


plot_grid(plt_as_a, new2, nrow=1, rel_widths = c(1,0.5))
ggsave('figures/Figure3_alternate_State.pdf', width=10, height=8)

plot_grid(plt_as_b, new1, nrow=1, rel_widths = c(1,0.5))
ggsave('figures/Figure3_alternate_Age.pdf', width=10, height=8)

