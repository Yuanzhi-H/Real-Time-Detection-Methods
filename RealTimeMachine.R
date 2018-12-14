RealTimeMachine <- function(DATBASE,Method,wset,RW,minfval,mindiff,runmode) {
# The order of the 16 progressive wear runs was SAMPLED in random and the last five runs suffered from broken tool 
randruns=c(1:4,8,10:12,15:16,18:19,6:7,14,20,5,9,13,17,21)
ntrain=12
ntests=4
nbroke=5
# Check the mode and select all the runs we use for training/testing
      if (runmode==1){randruns=randruns[1:ntrain]
}else if (runmode==2){randruns=randruns[(ntrain+1):(ntrain+ntests)]
}else if (runmode==3){randruns=randruns[(ntrain+ntests+1):(ntrain+ntests+nbroke)]}

nvar=3        # Number of forces
Hz=10         # There are 10 observations (90% quantiles) per second
maxobs=30*Hz; # We need at least 300 observations to start real-time detection
multicomb=as.logical(c(1,0,1)) # Choose the 1st and 3rd force data streams for the multivariate version
lifecoef1=5; lifecoef2=30;     # Two parameters for the tool wear prediction score function
mscore=4

# V: spindle speed. F: feed rate. 
V=c(40,40,30,35,55,45,25,28,48,58,38,33,53,43,23,24,44,54,34,39,59)
F=c(40.00,30.00,50.00,35.00,55.00,25.00,45.00,32.50,52.50,22.50,42.50,27.50,47.50,37.50,57.50,38.75,58.75,28.75,48.75,23.75,43.75)
# We remove some initial observations in each new session
burninset=round(120/(V+F),1)
# Define the time intervals when no cutting was performed
invalid=list(c(0,6.5),c(0,11.9),c(0,5.4,198,209.5,388.8,399.6,566.1,576.1,728,741),c(0,5.4),c(0,9.6),c(0,12.1),
c(0,5.4,216.9,229.9,426.5,440.1,622.6,633.8),c(0,10.8,213.3,228.4),c(0,4.4,63.8,73.5),c(0,7.2),c(0,5.8),c(0,7.5),c(0,5.5),c(0,7),
c(0,6.5,97.5,114.9,170.6,187,409.6,424.4),c(0,6.2,277.8,288.5),c(0,9.4),c(0,5.8),c(0,4.4,106,115.7),c(0,7.2),c(0,8.6,55.6,71.6))
# Time to finish the whole process and tool wear measured for each run
finalset    =c(126,163,879,236,52,141,767,570,82,63,116,272,42,95,629,526,50,48,203,190,159)
toolwear    =c(223.00,217.95,211.24,215.88,NA,202.01,214.42,228.27,NA,201.53,264.09,212.99,NA,224.93,177.65,221.59,NA,187.68,220.63,199.62,NA)

# Estimate the time when tool wear reached 150 micrometers 
toolwearthreshold=150
noncontact=numeric(length(finalset))
for(run in 1:length(finalset)){noncontact[run]=sum(diff(invalid[[run]])[c(TRUE,FALSE)])}
brokeset=round((finalset-noncontact)*toolwearthreshold/toolwear,1)
brokeset[15]=round((409-sum(diff(invalid[[15]][invalid[[15]]<409])[c(TRUE,FALSE)]))*toolwearthreshold/toolwear[15],1)
for(run in which(!is.na(brokeset))){delimit=invalid[[run]]
while (brokeset[run]>delimit[1]){brokeset[run]=brokeset[run]+delimit[2]-delimit[1];delimit=delimit[c(-1,-2)];if (length(delimit)==0) break}}
brokeset[is.na(brokeset)]=c(45,75,20,30,55) # five cases with broken tool

# Save the results: weight, number of runs where anomalies were detected, mean difference in distance, and mean score
out=matrix(NA,length(wset),1+6)
colnames(out)=c("W","Num1","Difference1","Score1","Num2","Difference2","Score2")
out[,1]=wset
nrun=length(randruns)

for(wi in 1:length(wset)){
tab=matrix(NA,nrun+1,nvar+1+4)
for(run in 1:nrun){
burninlimit=burninset[randruns[run]]
ntime=finalset[randruns[run]]*Hz-1
DAT=DATBASE[[randruns[run]]][1:(ntime+1),]
MD=numeric(ntime);MM=numeric(ntime)

for(rvar in 1:(nvar+1)){
noutliers=0
w=wset[wi]
if (rvar!=nvar+1){
rdat=DAT[,rvar]
# Calculate the first-order difference and minimum successive difference
M=diff(rdat)
D=c(0,M[-ntime])
M=cbind(D,-M)
M=M[cbind(1:ntime,max.col(-abs(M)))]
if (multicomb[rvar]){MD=MD+abs(D);MM=MM+abs(M)}
}else            {
rdat=rowMeans(DAT[,multicomb])
# Calculate the Euclidean Distance D in the multivariate case
M=sqrt(rowSums(diff(as.matrix(DAT[,multicomb]))^2))
D=c(0,M[-ntime])
M=pmin(D,M)
# Use half the weight if we combine two forces
w=w/sum(multicomb)
}

delimit=invalid[[randruns[run]]]
use=numeric(0)
nobs=nwarns=0
WARN=FALSE
for(i in 1:floor(ntime/Hz)){
if (sum(delimit)>0){if (i==ceiling(delimit[2])+1){correctdiff=median(rdat[max(delimit[1]*Hz+1,delimit[2]*Hz-3*Hz+1):(delimit[2]*Hz)])
rdat[(delimit[2]*Hz+1):(ntime+1)]=rdat[(delimit[2]*Hz+1):(ntime+1)]-correctdiff}}
if (WARN) {nwarns=nwarns+1;WARN=FALSE} else nwarns=max(0,nwarns-0)
if (sum(delimit)>0){if (i> delimit[1]-burninlimit && i<=delimit[2]+burninlimit) next else if (i> delimit[2]+burninlimit) delimit=delimit[c(-1,-2)]}

use=c(use,(i*Hz-Hz+1):(i*Hz))
nobs=nobs+Hz
if (nobs<30) next

MS=M[use]; DS=D[use]
      if (Method%/%10==1&&Method%%10==1){stat=abs(DS)
}else if (Method%/%10==1&&Method%%10==2){stat=abs(MS)
}else if (Method%%10==1) {stat=DS
}else if (Method%%10==2) {stat=MS}

if (RW>0){stat=tail(stat,min(nobs,RW*Hz));w=1*w}
new=tail(stat,Hz)
new=new[order(abs(new),decreasing=T)]
score=0
vec=round(c(0.1*Hz,0.2*Hz,0.4*Hz,0.6*Hz))
if (max(abs(tail(rdat[use],Hz)))>minfval || abs(new[vec[1]])>mindiff){
if (Method%/%10==1&& new[vec[1]]>(1.128+w*3*0.8525)/1.128*mean(stat)) score=score+1
if (Method%/%10==1&& new[vec[2]]>(1.128+w*2*0.8525)/1.128*mean(stat)) score=score+1
if (Method%/%10==1&& new[vec[3]]>(1.128+w*1*0.8525)/1.128*mean(stat)) score=score+1
if (Method%/%10==1&& new[vec[4]]>(1.128+w*0*0.8525)/1.128*mean(stat)) score=score+1
if (Method%/%10==2&&(mean(stat)-new[vec[1]])^2/var(stat)>(w*3)^2) score=score+2
if (Method%/%10==2&&(mean(stat)-new[vec[2]])^2/var(stat)>(w*2)^2) score=score+1
if (Method%/%10==2&&(mean(stat)-new[vec[3]])^2/var(stat)>(w*1)^2) score=score+1
if (Method%/%10==3&& GESD(stat,3,vec[1],Hz,RW,w)) score=score+2
if (Method%/%10==3&& GESD(stat,2,vec[2],Hz,RW,w)) score=score+1
if (Method%/%10==3&& GESD(stat,1,vec[3],Hz,RW,w)) score=score+1
}
if (score>=mscore){
noutliers=noutliers+1
if (nwarns>=0) {tab[run,rvar]=i;break}
WARN=TRUE}
}}

# Calculate results in the   univariate case below
if (sum(!is.na(tab[run,1:nvar]))>0){
delimit=invalid[[randruns[run]]]
if (brokeset[randruns[run]]>=sort(tab[run,1:nvar])[1]){
      delimit=delimit[delimit<brokeset[randruns[run]]];delimit= delimit[delimit>sort(tab[run,1:nvar])[1]]
}else{delimit=delimit[delimit>brokeset[randruns[run]]];delimit=-delimit[delimit<sort(tab[run,1:nvar])[1]]}
if (sum(delimit)!=0){tab[run,nvar+2]=brokeset[randruns[run]]-sort(tab[run,1:nvar])[1]-sum(diff(delimit)[c(TRUE,FALSE)])
}else tab[run,nvar+2]=brokeset[randruns[run]]-sort(tab[run,1:nvar])[1]
tab[run,nvar+2]=V[run]*tab[run,nvar+2]/60
x=tab[run,nvar+2]
      if (x< -lifecoef1){tab[run,nvar+3]=max(-(x+lifecoef1)^2/(lifecoef2-lifecoef1)^2+1,0)
}else if (x>  lifecoef1){tab[run,nvar+3]=max(-(x-lifecoef1)^2/(lifecoef2-lifecoef1)^2+1,0)
}else                    tab[run,nvar+3]=1
}
# Calculate results in the multivariate case below
if (    !is.na(tab[run,nvar+1])   ){
delimit=invalid[[randruns[run]]]
if (brokeset[randruns[run]]>=     tab[run,nvar+1]    ){
      delimit=delimit[delimit<brokeset[randruns[run]]];delimit= delimit[delimit>     tab[run,nvar+1]    ]   
}else{delimit=delimit[delimit>brokeset[randruns[run]]];delimit=-delimit[delimit<     tab[run,nvar+1]    ]}
if (sum(delimit)!=0){tab[run,nvar+4]=brokeset[randruns[run]]-     tab[run,nvar+1]    -sum(diff(delimit)[c(TRUE,FALSE)])
}else tab[run,nvar+4]=brokeset[randruns[run]]-     tab[run,nvar+1]    
tab[run,nvar+4]=V[run]*tab[run,nvar+4]/60
x=tab[run,nvar+4]
      if (x< -lifecoef1){tab[run,nvar+5]=max(-(x+lifecoef1)^2/(lifecoef2-lifecoef1)^2+1,0)
}else if (x>  lifecoef1){tab[run,nvar+5]=max(-(x-lifecoef1)^2/(lifecoef2-lifecoef1)^2+1,0)
}else                    tab[run,nvar+5]=1
}
}

tab[nrun+1,nvar+3]=sum(tab[,nvar+3],na.rm=T);tab[,nvar+3]=round(tab[,nvar+3]*10000)/10000
tab[nrun+1,nvar+5]=sum(tab[,nvar+5],na.rm=T);tab[,nvar+5]=round(tab[,nvar+5]*10000)/10000
out[wi,2]=sum(!is.na(tab[1:nrun,nvar+3]));out[wi,3]=mean(tab[1:nrun,nvar+2],na.rm=T)
out[wi,4]=tab[nrun+1,nvar+3]
out[wi,5]=sum(!is.na(tab[1:nrun,nvar+5]));out[wi,6]=mean(tab[1:nrun,nvar+4],na.rm=T)
out[wi,7]=tab[nrun+1,nvar+5]
}
return(list(out,tab))}