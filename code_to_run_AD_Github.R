# Use the code below for the real-time simulation of the turning process
# Last edited on 20/11/2018
# Yuanzhi Huang

# Clear workspace and load in packages ------------------------------------
rm(list=ls())
require(ggplot2)

# Load data ---------------------------------------------------------------
DATBASE=readRDS(file="Turn_DATBASE.rds")

# Source in functions -----------------------------------------------------
source("RealTimeMachine.R") # the main function for the test
source("runGESD.R") # function for the GESD test

# Run "RealTimeMachine" on the training data ------------------------------
# Define the candidate weights
wset=seq(0.8,3,0.2)
# Define the threshold-based rule
minfval=0.5;mindiff=minfval*0.1
ntrain=12 # number of runs for training
ntests= 4 # number of runs for testing
nbroke= 5 # number of runs with broken tool
ntotal=ntrain+ntests+nbroke

# Create a new csv file to save all the results
newcsv=TRUE
# Turn on the training mode (runmode=1) and take data from the 12 runs of progressive wear
runmode=1
denominator=ntrain
# Choose the detection method (first digit of Methodset=1,2,3) and the data (second digit of Methodset=1,2) to use
Methodset=c(11,12,21,22)
# Choose a team streaming scheme
for(RW in c( 0,60)){
# If we use the rolling window scheme, it is feasible to run the GESD test
if (RW!=0){Methodset=c(Methodset,31,32)}
for(Method in Methodset){
Result=RealTimeMachine(DATBASE,Method,wset,RW,minfval,mindiff,runmode)
out=Result[[1]]
# Calculate the mean of the tool wear prediction scores
out[,c(4,7)]=round(out[,c(4,7)]/denominator,2)
tosave=cbind(out[,1],out[,4],out[,7],Method%/%10,RW,Method%%10)
if (newcsv){
colnames(tosave)=c("W","Score1","Score2","Method","Stream","Difference")
write.table(tosave,"out.csv",F,F,",",row.names=F,col.names=T)
newcsv=FALSE
}else{
write.table(tosave,"out.csv",T,F,",",row.names=F,col.names=F)
}
# save the results:
# 1. weight
# 2. mean score in the   univariate case
# 3. mean score in the multivariate case based on the Euclidean distance
# 4. detection method (1 for statistical process control; 2 for chi-squared test; 3 for GESD test)
# 5. time streaming (0 for rolling origin; 60 for rolling window of 60 seconds each)
# 6. data (1 for first-order difference; 2 for minimum successive difference)
}}

# Test performance --------------------------------------------------------
# Use the rolling origin scheme
RW=0
# Use statistical process control method for the first-order difference data
Method=11 
# Turn on the test mode
runmode=2
denominator=ntests
Result=RealTimeMachine(DATBASE,Method,wset,RW,minfval,mindiff,runmode)
out=Result[[1]]
out[,c(4,7)]=round(out[,c(4,7)]/denominator,2)
write.table(cbind(out[,1],out[,4],out[,7],Method%/%10,RW,Method%%10),"out_tests.csv",T,F,",",row.names=F,col.names=F)

# Run a second test on the broken tool cases (OPTIONAL)
runmode=3
denominator=nbroke
Result=RealTimeMachine(DATBASE,Method,wset,RW,minfval,mindiff,runmode)
out=Result[[1]]
out[,c(4,7)]=round(out[,c(4,7)]/denominator,2)
write.table(cbind(out[,1],out[,4],out[,7],Method%/%10,RW,Method%%10),"out_broke.csv",T,F,",",row.names=F,col.names=F)

# Plot results ---------------------------------------------------------
# FIG. 5 of the paper
# Load the results to plot the mean scores across different scenarios
allscores=read.csv("out.csv",T,)
allscores[,4]=factor(allscores[,4])
allscores[,5]=factor(allscores[,5])
allscores[,6]=factor(allscores[,6])
levels(allscores[,4])[1]= "Control"
levels(allscores[,4])[2]= "Chi2 Test"
levels(allscores[,4])[3]= "GESD Test"
levels(allscores[,5])[1]= "RO"
levels(allscores[,5])[2]= "RW=60"
levels(allscores[,6])[1]= "FOD"
levels(allscores[,6])[2]= "MSD"
# Mean scores in the   univariate case
qplot(W    ,Score1,data=allscores,color=Method,shape=Method,facets=Stream~Difference,geom="point",
xlab="Weight",ylab="Mean Performance Score",ylim=c(0,0.75))+geom_line(aes(linetype=Method))+ 
theme_bw()
# Mean scores in the multivariate case
qplot(W*0.5,Score2,data=allscores,color=Method,shape=Method,facets=Stream~Difference,geom="point",
xlab="Weight",ylab="Mean Performance Score",ylim=c(0,0.75))+geom_line(aes(linetype=Method))+ 
theme_bw()

# FIG. 6 of the paper
# Below we plot and examine how accurate our predictions are
# Define the weight
wset=1.2
# Turn on the mode when data of all the 21 runs are fed to "RealTimeMachine"
runmode=4
# The vector below contains the time when tool wear has reached 150 micrometres
brokeset=c(86.9,115.9,625.9,165.6,45.0,107.8,538.3,383.4,75.0,48.7,68.4,193.8,20.0,65.7,351.6,361.5,30.0,39.5,142.5,144.6,55.0)
# The vector below indicates the order of the 21 runs, which was chosen in random
randruns=c(1:4,8,10:12,15:16,18:19,6:7,14,20,5,9,13,17,21)
brokeset=brokeset[randruns]
Result=RealTimeMachine(DATBASE,Method,wset,RW,minfval,mindiff,runmode)
Test=as.factor(c(0,0,0,0,2,1,1,0,2,0,0,0,2,1,0,0,2,0,0,1,2))
Test=Test[randruns]
levels(Test)=c("No","Yes","Broken")
Score=round(Result[[2]][-nrow(Result[[2]]),8],2)
qplot(brokeset,Result[[2]][-nrow(Result[[2]]),4],color=Score,shape=Test,xlab="Tool Wear Time (Seconds)",ylab="Predicted Time (Seconds)",xlim=c(0,650),ylim=c(0,650))+ 
geom_abline(intercept=0,slope=1,lty=2,col="red")+
geom_point(size=2)+ 
theme_bw()