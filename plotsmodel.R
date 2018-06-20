#################################################
# Source: SEIR Network Model                      
# Author: Kimia Ameri (hameri@unomaha.edu)      
# Date:   Unknown             
# Edited: 06.19.2018 by Kate Cooper kmcooper@unomaha.edu                      
#################################################

#Load required libraries
library("EpiDynamics")
library("deSolve")
library("forecast")
library("igraph")

#######################################
###     User-defined Variables      ###
#######################################
workingDir <- ("/Users/katedempsey/git/seir/SEIR/")
year <- 2015
reportsFile <- "/Users/katedempsey/git/seir/SEIR/reports.csv"
populationFile <-"/Users/katedempsey/git/seir/SEIR/population.csv"
m <- 13     #The constant out-degree of variables in the AB network model
power <- 1  #The power of the AB network model

#######################################
###    Initialize other  Variables  ###
#######################################
I<-0  #I is the number of infected individuals 

#######################################
###     Input data                  ###
#######################################
input<- read.table(reportsFile,sep=",",header=TRUE)
input<- as.matrix(input)
Nebraska_population_allyears<- as.matrix(read.table(populationFile,sep=",",header=TRUE),ncol=2)

#What is going on here?
#Why 1:18?
#Does Ntemp mean the population? Can we rename it? I thought it was temperature.
for(i in 1:18){
  if(Nebraska_population_allyears[i,1]==year)
  {
    infected_population_temp<- as.numeric(input[i,2:54])
    Nebraska_population_temp<- as.numeric(Nebraska_population[i,2])
  }
}


########################################
### Make a Barabasi network model    ###
########################################
net<-barabasi.game(
  n=Nebraska_population_temp, # The number of nodes is the Nebraska population for that year
  m=m,                     # The constant out-degree of variables
  directed = FALSE,
  power=power,
)

x<-c()
for (i in 1:1000)
{
  r<-as.integer(runif(1,1,Ntemp))
  ne<-neighbors(net,r,mode ="all")
  x<-c(x,length(ne))
  
}

(E<-as.integer(mean(x)))
  
  ########################################
  ###            set parameters        ###
  ########################################
  I<-inputtemp[1] 
  R<-0 
  S<-(Ntemp-I-R)* 0.90
  N<- as.integer(S+E+I+R)
  mu<-0.0000085   # population rate
  beta<- 0.9      # susceptible-infected transmission to exposure
  sigma <- 0.5    # transmission from exposed to infected
  gamma<- 0.2     # transmission from infected to recovered
  init       <- c(S=S, E=E, I=I, R=R)
  parameters <- c(beta=beta, gamma=gamma, mu=mu, sigma=sigma)
  times      <- seq(0, 52, by = 1)
  #########################################
  ###   Create an SEIR function         ###
  #########################################
  NET.SEIR <- function(time, state, parameters) {
    
    with(as.list(c(state, parameters)), {
      
      dS=(mu*(N-S))-(beta*((S*I)/N))-S
      dE=(beta*(S*I)/N)-((mu+sigma)*E)
      dI=(sigma*E)-((mu+gamma)*I)
      dR=(gamma*I)-(mu*R)+S
      
      return(list(c(dS, dE, dI, dR)))
    })
  }
  #########################################
  ###      calculate accurary           ###
  #########################################
  simulation.net<-as.data.frame(lsoda(init,times,NET.SEIR,parameters))
  predictednet <-as.integer(simulation.net$I)
  accresnet<-accuracy(predictednet,inputtemp)
 # write.table(accresnet[2], file = "RMSEvariant.csv",append = TRUE, row.names=FALSE, col.names=FALSE, sep=",")
#}
#########################################
###    Create Scalefree SEIR plot     ###
#########################################
pdf(paste("NB-SEIR for.pdf",year))
plot(simulation.net$I,main=paste("NB-SEIR Model for",year),xlab="Time (weeks)", ylab="Number of Infected",type="l",col="red",ylim=c(0,max(simulation.net)))
lines(simulation.net$S,type="l",col="green",ylim=c(0,max(simulation.net)))
lines(simulation.net$E,type="l",col="orange",ylim=c(0,max(simulation.net)))
lines(simulation.net$R,type="l",col="blue",ylim=c(0,max(simulation.net)))
legend('topright', legend=c('Susceptible','Exposed', 'Infected', 'Recovered'), col=c('green','orange','red','blue'), lwd=2.5)
dev.off()
#########################################
###              SEIR                 ###
#########################################
I<-inputtemp[1] 
R<-0 
S<-(Ntemp-I-R)* 0.90
E<-(Ntemp-I-R)* 0.10
N<- as.integer(S+E+I+R)
mu<-0.0000085   # population rate
beta<- 0.9      # susceptible-infected transmission to exposure
sigma <- 0.5    # transmission from exposed to infected
gamma<- 0.2     # transmission from infected to recovered
init       <- c(S=S, E=E, I=I, R=R)
parameters <- c(beta=beta, gamma=gamma, mu=mu, sigma=sigma)
times      <- seq(0, 52, by = 1)
#########################################
###   Create an SEIR function         ###
#########################################
seir <- function(time, state, parameters) {
  
  with(as.list(c(state, parameters)), {
    
    dS=(mu*(N-S))-(beta*((S*I)/N))-S
    dE=(beta*(S*I)/N)-((mu+sigma)*E)
    dI=(sigma*E)-((mu+gamma)*I)
    dR=(gamma*I)-(mu*R)+S
    
    return(list(c(dS, dE, dI, dR)))
  })
}
#########################################
###      calculate accurary           ###
#########################################
simulation<-as.data.frame(lsoda(init,times,seir,parameters))
predictedk <-as.integer(simulation$I)
accres<-accuracy(predictedk,inputtemp)
#########################################
###          Create SEIR plot         ###
#########################################
pdf(paste("SEIR Plot for.pdf",year))
plot(simulation$I,main=paste("SEIR Model for",year),xlab="Time (weeks)", ylab="Number of Infected",type="l",col="red",ylim=c(0,max(simulation)))
lines(simulation$S,type="l",col="green",ylim=c(0,max(simulation)))
lines(simulation$E,type="l",col="orange",ylim=c(0,max(simulation)))
lines(simulation$R,type="l",col="blue",ylim=c(0,max(simulation)))
legend('topright', legend=c('Susceptible','Exposed', 'Infected', 'Recovered'), col=c('green','orange','red','blue'), lwd=2.5)
dev.off()

#########################################
###          Create LOGplot           ###
#########################################
pdf("logplot.pdf")
plot(as.integer(log(simulation$I+1,base=10)),main=paste("SEIR Vs. NB-SEIR Vs. RealData for",year),xlab="Time (weeks)", ylab="Number of Infected ",type="l",col="red",lty=1,ylim=c(0,max(log(simulation$I))))
lines(as.integer(log(simulation.net$I+1,base=10)),type="l",col="darkgreen",ylim=c(0,max(log(simulation$I+1))))
lines(as.integer(log(inputtemp+1,base=10)),type="l",col="blue",ylim=c(0,N),lty=1)
legend('topright', legend=c('SEIR Model','NB-SEIR', 'Real data'), col=c('darkgreen','red','blue'),lty = c(1,1,1))
dev.off()

pdf(paste("LOGplots/LOG comparition for.pdf",year))
plot(as.integer(log(simulation$I+1,base=10)),main=paste("SEIR Vs. Scale-free Network Vs. RealData for",year),xlab="Time (weeks)", ylab="Number of Infected ",type="l",col="red",lty=1,ylim=c(0,max(log(simulation$I))))
lines(as.integer(log(simulation.net$I+1,base=10)),type="l",col="darkgreen",ylim=c(0,max(log(simulation$I+1))))
lines(as.integer(log(inputtemp+1,base=10)),type="p",col="blue",ylim=c(0,N),pch=1)
legend('topright', legend=c('SEIR Model','ScaleFree Network', 'Real data'), col=c('darkgreen','red','blue'),lty = c(1,1,3))
dev.off()
#########################################
###    actual Vs Scale-Free  plot     ###
#########################################
pdf(paste("actual Vs NB-RMSE network for.pdf",year))
plot(as.integer(inputtemp),main=paste("Actual infected Vs NB-RMSE prediction for",year),xlab="weeks",ylab="Number of Infected",type="l",col="blue",lty=1,ylim=c(0,max(simulation.net$I)))
lines(as.integer(simulation.net$I),type="l",las=2,col="red",ylim=c(0,N),pch=1)
legend('topright', legend=c('Actual','NB-SEIR'), col=c('blue','red'),lty = c(1,1))
dev.off()

##########################################
###             RMSE variance          ###
##########################################
RMSE<- as.matrix(read.table("totalRMSE.csv",sep=",",header=TRUE,check.names = F))
RMSEsorted <- RMSE[,order(colnames(RMSE))]
pdf("boxplot.pdf")
boxplot(RMSEsorted,las=2,xlab="Year", ylab="NB-SEIR",main= "SEIR Variance for NB-SEIR for 100 iteration")
dev.off()
#########################################
outputs<- read.table("results.csv",sep=",",header=TRUE)
