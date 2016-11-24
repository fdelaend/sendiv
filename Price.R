require(lattice)

i <- 8 #study considered
ResultsFolder <- "/Users/frederik/Documents/results/sendiv/" #results folder

#1. Parameters #####################################################################
Folder            <- c("~/Dropbox/sendiv/Hartgers/",
                       "~/Dropbox/sendiv/Brock/",
                       "~/Dropbox/sendiv/Roessink/",
                       "~/Dropbox/sendiv/VandenBrink/",
                       "~/Dropbox/sendiv/md2/",
                       "~/Dropbox/sendiv/md3/",
                       "~/Dropbox/sendiv/geest/",
                       "~/Dropbox/sendiv/pvdb/") #Folder with biovolume and count data
CountCols         <- c(4,4,7,4,7,6,6,5)                       #Counts in raw data start at this column
Time              <- c("Week", "week", "week", 
                       "week", "Days.p.a.", "Days.p.a.",
                       "Week", "Week")            #Name used to indicate time in count data
StartDates        <- c(1, 9, 1, 1, 1, 1, 1, 1)                      #Average densities are calculated for Time>= this time
EndDates          <- c(4, 22, 10, 11, 21, 24, 11, 4)                     #...and for Time<= this time
Treat             <- c("Treatment", "Treatment", "treatment", 
                       "Treatment", "Treatment", "Treatment",
                       "Treatment", "Treatment")                  #Name used to indicate trt level in raw data
Replicate         <- c("Replicate", "Replicate", "replica", 
                       "Replicate", "Replicate", "Replicate",
                       "Replicate", "Replicate")                  #Name used to indicate cosm in raw data
ControlTreat      <- c(1,1,0.000,1,1,1,1,1)                            #Code used for control treatment

#2. Function for partition #####################################################################
Price <- function(tox, ref)
{
  sc <- sum(tox*ref>0)                            #nr of species both sites have in common
  s  <- sum(ref>0)                                #nr of species at ref site
  s_ <- sum(tox>0)                                #nr of species at tox site
  z  <- mean(ref[which(ref>0)])                   #mean function per species at ref site. Only for species that are present.
  z_ <- mean(tox[which(tox>0)])                   #mean function per species at tox site. Only for species that are present.

  #re-order
  ind<- which((tox*ref)>0)                        #indices of species present at both sites
  oth<- c(1:length(ref))[-ind]                    #other species
  tox<- tox[c(ind,oth)]
  ref<- ref[c(ind,oth)]

  #now get zc and z_c
  zc <- mean(ref[c(1:sc)])
  z_c<- mean(tox[c(1:sc)])

  #PRICE EQUATIONS
  SREL <- (sc-s)*z
  SREG <- (s_-sc)*z_
  SCEL <- sc*(zc-z)
  SCEG <- -sc*(z_c-z_)
  CDE  <- sc*(z_c-zc)
  return(list(BDE = c(SREL=SREL,
                      SREG=SREG,
                      SCEL=SCEL,
                      SCEG=SCEG,
                      CDE=CDE)))
}

#3. Read count and BioV data #####################################################################
CountData  <- read.delim(paste(Folder[i],"Data/phytoplankton.txt",
                               sep=""),
                         header=TRUE)
Counts <- data.matrix(CountData[,c(CountCols[i]:ncol(CountData))]) #get counts from CountData
BioVols  <- read.delim(paste(Folder[i],"Biovolumes/phytoplankton/PVolumes.txt",
                               sep=""))
#Put the 'species' in BioVols in the same order as in Counts
BioVols   <- BioVols[match(colnames(Counts), BioVols[,"Sp"]),]

#4. Convert to Biov #####################################################################
Volumes <- Counts %*% diag(BioVols[,"Norm_Volume"])
colnames(Volumes) <- colnames(Counts)
VolumeData <- CountData[,c(Treat[i], Time[i], Replicate[i])]

#5. Make sum per genus #####################################################################
Genera <- unique(as.character(BioVols[,"Genus"])) #pull out different genera
Test <- sum(is.na(Genera)) #check for NAs
if (Test>0) {Genera <- Genera[-which(is.na(Genera))]} #throw them out
for (Genus in Genera)
{
  #the following species belong to genus
  SpeciesGenus <- as.character(BioVols[which(BioVols[,"Genus"]==Genus),"Sp"]) 
  #get their bioVs and sum if multiple species belong to genus
  VolumeGenus <- Volumes[,SpeciesGenus]
  if (length(SpeciesGenus) >1) {VolumeGenus <- rowSums(VolumeGenus)}
  #add to VolumeData
  VolumeData <- cbind(VolumeData, VolumeGenus)
}
#Now give proper names to VolumeData
colnames(VolumeData) <- c("Treat", "Time", "Replicate", Genera)

#6. Get out selected time period #####################################################################
#...and identify stress levels
Ind <- which((VolumeData[,"Time"]>=StartDates[i])&(VolumeData[,"Time"]<=EndDates[i]))
VolumeDataSelect <- VolumeData[Ind,]
#Make average densities over time
Key <- paste(VolumeDataSelect$Treat,
             VolumeDataSelect$Replicate,sep="") #to use to make averages - unique combi of cosm and treat
VolumeDataSelect <- aggregate(VolumeDataSelect, 
                              by=list(Key), FUN=mean)
Stress<- unique(VolumeDataSelect[,"Treat"])
#get control site
Ind <- which(VolumeDataSelect[,"Treat"]==ControlTreat[i])
VolumeDataSelectControl <- VolumeDataSelect[Ind,
                                            c(5:ncol(VolumeDataSelect))] 

#7. Loop over stress levels #####################################################################
Results <- NULL
for (l in Stress[2:length(Stress)]) #1 being the control
{
  Ind      <- which(VolumeDataSelect[,"Treat"]==l)
  VolumeDataSelectStress <- VolumeDataSelect[Ind,c(5:ncol(VolumeDataSelect))]
  #all possible combinations
  Combs <- expand.grid(REF=c(1:nrow(VolumeDataSelectControl)),
                       TOX=c(1:nrow(VolumeDataSelectStress)))
  #loop over them
  for (r in c(1:nrow(Combs)))
  {
    Ref <- VolumeDataSelectControl[Combs[r,"REF"],]
    Tox <- VolumeDataSelectStress[Combs[r,"TOX"],]
    #only keep species that are present in at least one site
    Keep<- which(((Ref>0)+(Tox>0))>0)
    Ref <- as.numeric(Ref[Keep])
    Tox <- as.numeric(Tox[Keep])
    Del <- sum(Tox) - sum(Ref) #deltaT
    #for fun, track proportion occupied by species lost in tox site
    Lost <- sum(Ref[which(Tox==0)])/sum(Ref) 
    #for fun, track proportion occupied by species gained in tox site
    Gain <- sum(Tox[which(Ref==0)])/sum(Tox) 
    Results <- rbind(Results, 
                     c(l,Delta=Del,
                       Price(Tox, Ref)[[1]],
                       Lost=Lost, 
                       Gain=Gain))
    
  }
}

colnames(Results)[1] <- "Stress"
Results <- as.data.frame(Results)

quartz("",6,6,type="pdf",
       file=paste(ResultsFolder,i,"Test.pdf",sep=""))
plot(0,0,xlim=c(min(Stress),max(Stress)+0.5), xlab="Stress intensity",
     ylab="Effect",
     ylim=c(min(Results),max(Results)+20),col="white")
col  <- c("black","blue","brown","orange","red","green")
pch <- c(17,rep(1,5))
sapply(2:7,function(x) points(Results[,1]+0.00001,Results[,x],
                              col=col[x-1],pch=pch[x-1])) 
legend("topleft",c("deltaT","SREL","SREG", "SCEL", "SCEG",
                   "CDE"),col=col,pch=pch, horiz=TRUE)
dev.off()

quartz("",6,6,type="pdf",
       file=paste(ResultsFolder,i,"Test2.pdf",sep=""))
plot(Results$CDE,Results$Delta)
abline(0,1)
dev.off()

#small additional analysis showing that stress increasingly filtered out species doing bad in ref site
plot(Results$Stress,Results$Lost)
#small additional analysis showing that stress increasingly led to appearance of species doing good in tox site
plot(Results$Stress,Results$Gain)



