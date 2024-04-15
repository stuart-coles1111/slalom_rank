library(stringr)
library(dplyr)
library(ggplot2)
load("~/slalom_rank/slalom.RData")
files <- c('~/slalom_rank/data/gurgl_2023_1.rtf','~/slalom_rank/data/madonna_di_campiglio_2023_2.rtf','~/slalom_rank/data/adelboden_2024_3.rtf','~/slalom_rank/data/wengen_2024_4.rtf','~/slalom_rank/data/kitzbuhel_2024_5.rtf','~/slalom_rank/data/schladming_2024_6.rtf',
           '~/slalom_rank/data/chamonix_2024_7.rtf','~/slalom_rank/data/Palisades_Tahoe_2024_8.rtf', '~/slalom_rank/data/aspen_2024_9.rtf', '~/slalom_rank/data/saalbach_2024_10.rtf')

races <- c('1: Gurgl','2: Madonna di Campiglio','3: Adelboden','4: Wengen','5: Kitzbuhel','6: Schladming','7: Chamonix','8: Palisades Tahoe', '9: Aspen', '10: Saalbach')

alldata <- c()

for(i in 1:10){ temp=read.csv(files[i])
temp$race=races[i]
temp$bib = nrow(temp) - rank(temp$H) +1
alldata=rbind(alldata,temp)
}
colnames(alldata)[1] <- 'A'

alldata$D <- str_trim(alldata$D)
alldata <- left_join(alldata,racers_df,join_by(D==racers))
#alldata$race <- factor(alldata$race, levels=1:10)

alldata$race=factor(alldata$race,levels=unique(alldata$race),order=T)

colnames(alldata) <- c('position','bib1','id','surname','name','dob','country','time1','rank1','time2','rank2','total_time','diff','race_points','ski','race','bib2','points')

temperature <- c(-4,0,-2,0,3,5,8,8,-2,11)
wind = c(rep('none',7),'med','med','none')
weather= c('sunny','cloudy','cloudy','cloudy','cloudy','rain','sunny','sunny','cloudy','cloudy')
snow = c('compact','compact','compact','hard','hard','compact','compact','hard','packed','spring')
