######### Rework the .sum file to mesh with the filtering algorithm
library(tidyverse)
setwd("C:/Users/gsinger/Dropbox/GitHubRepos/TeknoFilter")

#read in known tag list
tags <- read.csv("./taglist/FriantTaglist.csv")
#set counter
counter <- 1:12
#read in .SUM file
dat <- read.csv("./raw/2015-6003170601804.SUM", skip = 8, header=T)  
#rename columns to match db
names(dat)<- c("Detection", "RecSN", "dtf", "Hex", "Tilt", "Volt", "Temp", "Pres", "Amp",
               "Freq", "Thresh", "nbw", "snr", "Valid")              

#change object format to tbl
dat <- as.tbl(dat)                                                   
#drop the barker code and the CRC from tagid field
dat$Hex <- substr(dat$Hex, 4, 7)    
#name and print new tbl w/ bad detect lines removed
(dat <- dat %>%
  filter(Detection != "      -"))                                    
#filter receiver file by known taglist
dat<- dat[dat$Hex %in% tags$Tag.ID..hex., ]                          
#set nPRI (Nominal PRI) for the tag 
dat$nPRI<- 5     
#convert to POSIXct note: fractional secons will no longer print, but they are there
dat$dtf<- as.POSIXct(dat$dtf, format = "%m/%d/%Y %H:%M:%OS", tz="Etc/GMT-8")                                 
#run the next line to verify that you haven't lost your frac seconds
#(strftime(dat$dtf, format = "%m/%d/%Y %H:%M:%OS6"))

dat2 <- dat

dat2<- arrange(dat2, Hex, dtf)#sort by TagID and then dtf, Frac Second
#calculate tdiff, then remove multipath
dat3 <- data.frame(dat2, tdiff=c(NA, difftime(dat2$dtf[-1], dat2$dtf[-nrow(dat2)])), 
                   winmax=dat2$dtf+((dat2$nPRI*1.3*max(counter))+1))
dat4 <- data.frame(dat3, crazy=c(NA,dat3$Hex[-nrow(dat3)]==dat3$Hex[-1]))
dat4$tdiff[dat4$crazy==0] <- NA
dat5 <- dat4[,-18]
dat5 <- dat5[dat5$tdiff>0.2 | is.na(dat5$tdiff),]


###Figure out what got dropped between dat4 and dat5
head(dat4)
head(dat5)
dropped <- dat4[!paste(strftime(dat4$dtf, format = "%m/%d/%Y %H:%M:%OS6"), dat4$Hex) %in% 
                     paste(strftime(dat5$dtf, format = "%m/%d/%Y %H:%M:%OS6"), dat5$Hex),]
#opened 'dropped' and dat4 to manually inspect the lines that were dropped due to multipath detections.
#everything appears to check out!  So, now I think that the .SUM file is ready to be fed through the filtering alorithm

