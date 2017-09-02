####################################################################################################################################
#                                                                                                                                  #
#                         Tag Filter for Lotek Receiver Files converted from CBR description                                       #
#                           Written by: Gabe Singer, Damien Caillaud     On: 05/16/2017                                            #
#                                            Last Updated: 06/22/2017                                                              #
#                   Issues: I think that I have some circular logic that is making the function keep dups                          #
####################################################################################################################################
#1)Rewrite script to run by Hex instead of Dec, no Dec in Tekno files - check
#2)Cycle through each unique tagID in the file - check
#3)End Prouct 2 df, one with all of the isAccecpted==T, an one with the isAccepted==F   - check
#4)Write output with original file name and filteredAccepted/filteredRejected '20176011_FilteredAccepted' 20176011_FilteredRejected - check
#5)Cycle through each file in a folder - check

###Eventually, turn this into a shiny app that lets you change some input parameters and see how your filtered data is affected

#setwd to whatever directory the receiver files are in

install.load <- function(package.name)
{
  if (!require(package.name, character.only=T)) install.packages(package.name)
  library(package.name, character.only=T)
}


install.load('tidyverse')
install.load('lubridate')
install.load('plotly')

###############################Data loading and cleaning#########################################################
###Working directory setup
#***folder names are non-negotiable
#In your wd you should have the following folders: 1)'raw' - holds raw .jst files
#                                                  2)'cleaned' - where the code will place your cleaned raw files
#                                                  3)'accepted' - where the code will place your good detections
#                                                  4)'rejected' - where the code will place you detections that were filtered out
#                                                  5)'taglist' - holds a .csv of known Tag IDs to filter by

setwd('C:/users/gabri/Dropbox/GitHubRepos/TeknoFilter/') #set to whatever wd you will work out of
tags<- read.csv("./taglist/FriantTaglist.csv", header = T) #list of known Tag IDs
dat<- read.csv("2017-6008DF.csv", header=F)#not adding extra stuff for now, just dealing with detections,
#also this file is already filtered for known tag ids***

dat<- read.csv("./raw/2015-6003170601804.jst", header=F) #load data 
names(dat)<- c("Filename", "RecSN", "DT", "FracSec", "Hex", "Junk", "valid", "junk2", "junk3") #rename headers
dat$nPRI<- 5 #nPRI = 'Nominal PRI' (ie set the pulse rate for the tag)
head(dat) #check changes in file

#combine the DT and FracSec columns into a single time column and convert to POSIXct
dat2<- dat #copy dat, so if we screw up we on't have to load data again
dat2$dtf<- paste0(dat2$DT, substring(dat2$FracSec,2)) #paste the fractional seconds to the end of the DT in a new column
dat2$dtf<- as.POSIXct(dat2$dtf, format = "%m/%d/%Y %H:%M:%OS", tz="Etc/GMT-8") #convert to POSIXct beware this may change value of 0.0000X
#head(strftime(dat2$dtf, format = "%m/%d/%Y %H:%M:%OS5")) #verify that although the fractional seconds don't print, they are indeed there
head(dat2) #check changes made to data frame
class(dat2$dtf)

#sort by TagID and then dtf, Frac Second
dat2<- as.tbl(dat2) #change from data.frame to tbl (tidyverse/dplyr package) for easy manipulation
dat2$Hex <- as.character(dat2$Hex) #change the class of the tag id from hexmode to character
dat2<- arrange(dat2, Hex, dtf) #sort tbl by TagID and then by datetime with frac seconds
dat2 #check changes



# ###Loop for cleaning raw .jst files
# for(i in list.files("./raw")){    #./raw is whatever folder your files are in
#       dat <- read.csv(paste0("./raw/", i), header=F)        #read in each file
#       names(dat)<- c("Filename", "RecSN", "DT", "FracSec", "Hex", "Junk", "valid", "junk2", "junk3") #rename headers
#       dat$nPRI<- 5 #nPRI = 'Nominal PRI' (ie set the pulse rate for the tag)
#       
#       #combine the DT and FracSec columns into a single time column and convert to POSIXct
#       dat2<- dat #copy dat, so if we screw up we on't have to load data again
#       dat2$dtf<- paste0(dat2$DT, substring(dat2$FracSec,2)) #paste the fractional seconds to the end of the DT in a new column
#       dat2$dtf<- as.POSIXct(dat2$dtf, format = "%m/%d/%Y %H:%M:%OS", tz="Etc/GMT-8") #convert to POSIXct beware this may change value of 0.0000X
#       #head(strftime(dat2$dtf, format = "%m/%d/%Y %H:%M:%OS5")) #verify that although the fractional seconds don't print, they are indeed there
#       
#       #sort tbl
#       dat2<- as.tbl(dat2) #change from data.frame to tbl (tidyverse/dplyr package) for easy manipulation
#       dat2$Hex <- as.character(dat2$Hex) #change the class of the tag id from hexmode to character
#       dat2<- arrange(dat2, Hex, dtf) #sort tbl by TagID and then by datetime with frac seconds
#       write.csv(dat2, paste0("./cleaned/", dat2$RecSN[1], "_cleaned.csv"), row.names=F)
}
#**********add code to filter by known TagID
###################Begin writing filtering algorithm based on the CBR description#####################################
#later, add filtering criteria that AA has, release date based, after expected tag failure date, etc

###create an index of time between subsequent detections (name this vector tdiff) assuming that subsequent tagIDs match
#also add a winmax column identifing the window search window.
###dplyr is way more efficient than a loop
#NominalPRI<- 5 #set value of nomimal PRI


# dat3 <- dat2 %>%
# group_by(Dec) %>%
# mutate(tdiff = difftime(dtf, lag(dtf))) %>%
# mutate(winmax = dtf+((NominalPRI*1.3*12)+1))

###write function to find the mode (or multiple modes)
mode <- function(x, i){
  ta <- table(x)
  tam <- max(ta)
  if (all(ta == tam))
    mod <- NA
  else
    if(is.numeric(x))
      mod <- as.numeric(names(ta)[ta == tam])
  else
    mod <- names(ta)[ta == tam]
  #print(i) 
  return(mod)
}
#  head(dat4)


#function has 3 inputs, 'dat' = your data file, 'tagHex' = Hex tagID that you are working with, and 'counter' = the intervals for hits (e.g 12)
magicFunc <- function(dat, tagHex, counter){
  dat3 <- data.frame(dat, tdiff=c(NA, difftime(dat2$dtf[-1], dat2$dtf[-nrow(dat2)])), winmax=dat2$dtf+((dat2$nPRI*1.3*max(counter))+1))
  ##get rid off crazy values when tagID changes
  head(dat3)
  dat4 <- data.frame(dat3, crazy=c(NA,dat3$Hex[-nrow(dat3)]==dat3$Hex[-1]))
  
  dat4$tdiff[dat4$crazy==0] <- NA
  dat5 <- dat4[,-14]
  
  dat5 <- dat5[dat5$tdiff>0.2 | is.na(dat5$tdiff),]
  
  tagid <- dat5[dat5$Hex==tagHex,]##use 2892 as an example
  logTable <- data.frame(hitRowNb=numeric(0), initialHitRowNb=numeric(0), isAccepted=logical(0), nbAcceptedHitsForThisInitialHit =logical(0))
  for(i in 1:nrow(tagid)){
    hits<- tagid[tagid$dtf>=tagid$dtf[i] & tagid$dtf<=tagid$winmax[i], ] #hits within window
    #head(tagid[tagid$dtf>=tagid$dtf[i], ], 100)
    #if nrow(hits>=4) then... (need to write code for this qualifier, if_else???)
    if(nrow(hits)>1){
      retained<- numeric(0) #make empty df to hold detections retained
      # indices<- numeric(0) #make empty df to hold detections retained
      
      for(j in 2:nrow(hits)){
        candidates<- round(as.numeric(hits$dtf[j]-hits$dtf[1])/counter, digits = 2) #subsubtract time of each det in window from initial and divide by 1:12
        candidates<- candidates[candidates>= tagid$nPRI[i]*0.651 & candidates<= tagid$nPRI[i]*1.3] #constrain values
        retained<- c(retained, candidates)#bind each batch of retained candiate PRIs to a single list
        # indices <- c(indices, rep(j, length(candidates)))
      }
      if(length(retained)==0){
        logTable <- rbind(logTable, data.frame(hitRowNb=NA, initialHitRowNb=i, isAccepted=FALSE, nbAcceptedHitsForThisInitialHit=0)) 
      } else {
        #plot(table(retained))
        ePRI<- min(mode(retained, i))
        if(is.na(ePRI)) ePRI <- min(retained)
        
        nbHits <- 1
        for(j in 2:nrow(hits)){	
          ii <- round(as.numeric(hits$dtf[j]-hits$dtf[1])/ePRI)
          uppb <- ii*ePRI+hits$dtf[1]+0.006+ii*0.006
          lowb <- ii*ePRI+hits$dtf[1]-(0.006+ii*0.006)
          nbHits <- nbHits+(hits$dtf[j]>= lowb & hits$dtf[j]<= uppb)
          logTable <- rbind(logTable, data.frame(hitRowNb=i+j-1, initialHitRowNb=i, isAccepted=(hits$dtf[j]>= lowb & hits$dtf[j]<= uppb), nbAcceptedHitsForThisInitialHit=NA))
        }
        logTable$nbAcceptedHitsForThisInitialHit[logTable$initialHitRowNb ==i] <- nbHits
      }
    } else {
      logTable <- rbind(logTable, data.frame(hitRowNb=NA, initialHitRowNb=i, isAccepted=FALSE, nbAcceptedHitsForThisInitialHit=0))     	
    }
    #print(paste("row", i,"done"))
    #readline("next\n")
  }
  return(logTable) 
}

res <- magicFunc(dat=dat2, tagHex='00CA', counter=1:12)
head(res, 100)
hist(dat5$tdiff[dat5$tdiff<20])

###create a loop that goes through all possible tagHex values, uses the function above
###'filterthresh' = the number of accepted hits in a window to record an event (eg 4)
dataFilter <- function(dat, filterthresh, counter){
  res <- dat[numeric(0),]
  timer <- 0
  for(i in unique(dat$Hex)){
    ans <- magicFunc(dat, tagHex=i, counter=1:12)
    ans[!is.na(ans$hitRowNb),]
    ans2 <- ans[ans$nbAcceptedHitsForThisInitialHit>= filterthresh,]
    keep <- c(ans2$hitRowNb[ans2$isAccepted], ans2$initialHitRowNb[ans2$isAccepted])
    keep <- keep[!duplicated(keep)]
    ans3 <- dat[dat$Hex==i,][keep,]
    res <- rbind(res, ans3)
    timer <- timer+1
    print(timer/length(unique(dat$Hex)))
  }
  return(res)
}
myResults <- dataFilter(dat=dat2, filterthresh=4, counter=1:12)
head(myResults)

dim(dat2)
dim(myResults)

#####which rows have we discarded?
rejecteds <- dat2[!paste(dat2$dtf, dat2$Hex) %in% paste(myResults$dtf, myResults$Hex),]
dim(dat2)
dim(myResults)
dim(rejecteds)

###save results
write.csv(myResults, paste0(myResults$RecSN[1], "_accepted.csv"), row.names=F)
write.csv(rejecteds, paste0(rejecteds$RecSN[1], "_rejected.csv"), row.names=F)

###cycling through file names: ****ISSUE - need to clean and sort data before running through the loops****
#list.files()
#list.files(recursive=T)
#list.files("./test")


for(i in list.files("./raw")){
  dat <- read.csv(paste0("./raw/", i), header=F)        #read in each file
  names(dat)<- c("Filename", "RecSN", "DT", "FracSec", "Hex", "Junk", "valid", "junk2", "junk3") #rename columns
  dat$nPRI<- 5   # set nPRI (Nominal PRI) for the tag 
  #combine the DT and FracSec columns into a single time column and convert to POSIXct
  dat2<- dat #copy dat, so if we screw up we on't have to load data again
  dat2$dtf<- paste0(dat2$DT, substring(dat2$FracSec,2)) #paste the fractional seconds to the end of the DT in a new column
  dat2$dtf<- as.POSIXct(dat2$dtf, format = "%m/%d/%Y %H:%M:%OS", tz="Etc/GMT-8") #convert to POSIXct beware this may change value of 0.0000X
  #head(strftime(dat2$dtf, format = "%m/%d/%Y %H:%M:%OS5")) #verify that although the fractional seconds don't print, they are indeed there
  dat2<- arrange(dat2, Hex, dtf)#sort by TagID and then dtf, Frac Second
  dat2<- as.tbl(dat2)
  dat2$Hex <- as.character(dat2$Hex)
  write.csv(dat2, paste0("./cleaned/", dat2$RecSN[1], "_cleaned.csv"), row.names=F)
  myResults <- dataFilter(dat=dat2, filterthresh=4, counter=1:12)
  rejecteds <- dat2[!paste(strftime(dat2$dtf, format = "%m/%d/%Y %H:%M:%OS5"), dat2$Hex) %in% 
                      paste(strftime(myResults$dtf, format = "%m/%d/%Y %H:%M:%OS5"), myResults$Hex),]
  write.csv(myResults, paste0("./accepted/", myResults$RecSN[1], "_accepted.csv"), row.names=F)
  write.csv(rejecteds, paste0("./rejected/", rejecteds$RecSN[1], "_rejected.csv"), row.names=F)
} #above loop runs, but has ups in the accepted file



dim(myResults)
dim(rejecteds)
dim(dat2)
dim(myResults)
dim(rejecteds)
##rerun at home
####Maybe copy cleaned up code over to new file to make it run smoother

dat<- read.csv("./raw/2015-6034170371415.jst", header=F)
rejecteds <- dat2[!paste(strftime(dat2$dtf, format = "%m/%d/%Y %H:%M:%OS5"), dat2$Hex) %in% 
                    paste(strftime(myResults$dtf, format = "%m/%d/%Y %H:%M:%OS5"), myResults$Hex),]
