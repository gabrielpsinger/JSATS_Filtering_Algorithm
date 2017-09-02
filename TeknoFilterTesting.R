#####why isn't it removing dups
#####Changing magicFunction, changed 

magicFunc <- function(datos, tagHex, counter){
  dat3 <- data.frame(datos, tdiff=c(NA, difftime(datos$dtf[-1], datos$dtf[-nrow(datos)])), winmax=datos$dtf+((datos$nPRI*1.3*max(counter))+1))
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

dataFilter <- function(datos, filterthresh, counter){
  res <- datos[numeric(0),]
  timer <- 0
  for(i in unique(dat$Hex)){
    ans <- magicFunc(datos, tagHex=i, counter=1:12)
    ans[!is.na(ans$hitRowNb),]
    ans2 <- ans[ans$nbAcceptedHitsForThisInitialHit>= filterthresh,]
    keep <- c(ans2$hitRowNb[ans2$isAccepted], ans2$initialHitRowNb[ans2$isAccepted])
    keep <- keep[!duplicated(keep)]
    ans3 <- datos[datos$Hex==i,][keep,]
    res <- rbind(res, ans3)
    timer <- timer+1
    print(timer/length(unique(datos$Hex)))
  }
  return(res)
}
myResults <- dataFilter(datos=dat2, filterthresh=4, counter=1:12)



###############################################################
#ggplotly test
set.seed(100)
d <- diamonds[sample(nrow(diamonds), 1000), ]
plot_ly(d, x = ~carat, y = ~price, color = ~carat,
        size = ~carat, text = ~paste("Clarity: ", clarity))
p <- ggplot(data = d, aes(x = carat, y = price)) +
  geom_point(aes(text = paste("Clarity:", clarity))) +
  geom_smooth(aes(colour = cut, fill = cut)) + facet_wrap(~ cut)
p
ggplotly(p)

# testing commits process
# testing commits process 2.2
# CAnt seem to get it to merge

