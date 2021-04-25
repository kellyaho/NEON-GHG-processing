###This script calculates three different ATM mixing ratios to use for the function

### LOOP for different ATM inputs (paired, smoothed, median) ####
x <- read.csv("input.csv") #input file should have sample ID, site, timestamp, and the raw NEON ATM mixing ratios 


predList = list()
for (site in unique(x$site)){ 
  
  pred <- data.frame(Sample.ID = x$Sample.ID[x$site==site],
                     timestamp = x$collectDate[x$site==site],
                     HS.mCH4.before = x$HS.mCH4.before[x$site==site],
                     HS.mCO2.before = x$HS.mCO2.before[x$site==site],
                     HS.mN2O.before = x$HS.mN2O.before[x$site==site],
                     smoothed.HS.mCH4.before = predict(loess(HS.mCH4.before ~ as.numeric(collectDate), 
                                                             data = x[which(x$site==site),], 
                                                             span=0.25, se = F), 
                                                       x[which(x$site==site),], se=F),
                     smoothed.HS.mCO2.before = predict(loess(HS.mCO2.before ~ as.numeric(collectDate), 
                                                             data = x[which(x$site==site),], 
                                                             span=0.25, se = F), 
                                                       x[which(x$site==site),], se=F),
                     smoothed.HS.mN2O.before = predict(loess(HS.mN2O.before ~ as.numeric(collectDate), 
                                                             data = x[which(x$site==site),], 
                                                             span=0.25, se = F), 
                                                       x[which(x$site==site),], se=F),
                     median.HS.mCH4.before = rep(median(x$HS.mCH4.before[x$site==site])),
                     median.HS.mCO2.before = rep(median(x$HS.mCO2.before[x$site==site])),
                     median.HS.mN2O.before = rep(median(x$HS.mN2O.before[x$site==site])),
                     HS.mCH4.after = x$HS.mCH4.after[x$site==site],
                     HS.mCO2.after = x$HS.mCO2.after[x$site==site],
                     HS.mN2O.after = x$HS.mN2O.after[x$site==site],
                     Temp.insitu = x$Temp.insitu[x$site==site],
                     Temp.equil = x$Temp.equil[x$site==site],
                     Alk.measured = x$Alkalinity.measured[x$site==site],
                     Volume.gas = x$Volume.gas[x$site==site],
                     Volume.water = x$Volume.water[x$site==site],
                     Bar.pressure = x$Bar.pressure[x$site==site])
  pred$site <- site
  predList[[site]] <- pred # add it to your list
}

write.csv(predList, "input_for_function.csv")