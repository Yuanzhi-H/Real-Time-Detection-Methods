GESD <- function(stat,coef,cutoff,Hz,RW,w){
n=min(max(Hz*10[RW==0],Hz*RW),length(stat))
# Standardise all the historical observations
stdstat=abs(stat-mean(stat))/sd(stat)
# Find the cutoff value in the data subset of size Hz
outlier=sort(tail(stdstat,Hz),TRUE)[cutoff]
# Calculate the critical value of the t distribution
tv=qt(1-(1-pt(coef,n))/(n-cutoff+1),df=(n-cutoff-1))
tv=w*tv
# See if the value is a remarkable outlier
test=outlier>(n-cutoff)*tv/sqrt((n-cutoff-1+tv^2)*(n-cutoff+1))
return(test)
}