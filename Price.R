require(lattice)
#-1. Parameters #####################################################################
Data <- "~/Dropbox/sendiv/Hartgers/Data/phytoplankton.txt" #raw data
CountCols <- 4 #counts in raw data start at this column
print(paste("ATTENTION: Data in",Data,"are assumed to be counts from column",CountCols,"to last column"))
VolumeData <- "~/Dropbox/sendiv/Hartgers et al/Biovolumes/phytoplankton/PVolumes.txt" #volume or weight per individual (cell). Need to have exactly same names as in raw data
Time <- "Week" #Name used to indicate time in raw data
Time_s<- 5.0 #selected time point for analysis
Treat <- "Treatment" #Name used to indicate trt level in raw data
ControlTreat <- 1 #Code used for control treatment
#END Parameters #####################################################################

#0. Function for partition
Price <- function(tox, ref)
{
  sc <- sum(tox*ref>0) #nr of species both sites have in common
  s  <- sum(ref>0)     #nr of species at ref site
  s_ <- sum(tox>0)     #nr of species at tox site
  z  <- mean(ref[which(ref>0)]) #mean function per species at ref site. Only for species that are present.
  z_ <- mean(tox[which(tox>0)]) #mean function per species at tox site. Only for species that are present.

  #re-order
  ind<- which((tox*ref)>0) #indices of species present at both sites
  oth<- c(1:length(ref))[-ind]#other species
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

#1. Read data
Data  <- read.delim(Data)

#2. Convert to biov: the taxon order in Data and 
#..in biovolumes are exactly the same
Counts <- data.matrix(Data[,c(CountCols:ncol(Data))])
#the biovolume data
vols   <- read.delim(VolumeData)
#the biovolume data, in the proper order
vols   <- vols[match(colnames(Counts), vols[,1]),]
Volumes <- Counts %*% diag(vols[,2])
Data[,c(CountCols:ncol(Data))] <- Volumes

#3. Get out time points and stress levels
Data_s<- Data[which(Data[,Time]==Time_s),]
Stress<- unique(Data_s[,Treat])
Data_s_1 <- Data_s[which(Data_s[,Treat]==ControlTreat),c(CountCols:ncol(Data_s))] #the control data (REF site)

#4. Loop over stress levels
Results <- NULL
for (l in Stress[2:length(Stress)]) #1 being the control
{
  ind      <- which(Data_s[,Treat]==l)
  Data_s_l <- Data_s[ind,c(CountCols:ncol(Data_s))]
  #all possible combinations
  Combs <- expand.grid(REF=c(1:nrow(Data_s_1)),
                       TOX=c(1:nrow(Data_s_l)))
  #loop over them
  for (r in c(1:nrow(Combs)))
  {
    Ref <- Data_s_1[Combs[r,"REF"],]
    Tox <- Data_s_l[Combs[r,"TOX"],]
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

quartz("",6,6,type="pdf",file="Test.pdf")
plot(0,0,xlim=c(min(Stress),max(Stress)+0.5), xlab="Stress intensity",
     ylab="Effect",
     ylim=c(min(Results),max(Results)+20),col="white")
col  <- c("black","blue","brown","orange","red","green")
pch <- c(17,rep(1,5))
sapply(2:7,function(x) points(Results[,1]+(x-3)/12,Results[,x],
                              col=col[x-1],pch=pch[x-1])) 
legend("topleft",c("deltaT","SREL","SREG", "SCEL", "SCEG",
                   "CDE"),col=col,pch=pch, horiz=TRUE)
dev.off()

#small additional analysis showing that stress increasingly filtered out species doing bad in ref site
plot(Results$Stress,Results$Lost)
#small additional analysis showing that stress increasingly led to appearance of species doing good in tox site
plot(Results$Stress,Results$Gain)





