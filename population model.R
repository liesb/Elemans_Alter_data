

require(deSolve)
rm(list=ls())



### Parameter values (defaults) ###
N0 <- 1.97E8  			#size US population - default=1.97E8
p <- 0.5					  #fraction of KIR2DL2-pos people in US-population - default=0.5
mu <- 0.0132				#death rate HIV- people (yr-1) - default=0.0132
B <- 5.5E6					#birth rate US population (yr-1) - default=5.5E6

beta.WT.preHAART <- 0.506    #transmission probability rate WT strain (yr-1) preHAART - default=0.506
beta.WT.postHAART <- 0.087	 #transmission probability rate WT strain (yr-1) postHAART - default=0.087
beta.V.preHAART <- 0.506		 #transmission probability rate V strain (yr-1) preHAART - default=0.506
beta.V.postHAART <- 0.087		 #transmission probability rate V strain (yr-1) postHAART - default=0.087

alpha.preHAART <- 0.1		     #death rate HIV+ people (yr-1) preHAART - default=0.1 
alpha.preHAART.PWT <- 0.04	 #death rate HIV+ people (yr-1) preHAART - HLA-matched, WT infected - default=0.04
alpha.postHAART <- 0.04	     #death rate HIV+ people (yr-1) postHAART - HLA-matched - default=0.04
alpha.postHAART.PWT <- 0.04	 #death rate HIV+ people (yr-1) postHAART - HLA-matched, WT infected - default=0.04


t0 <- 1966					   #time of introduction of HIV in US
now <- 2010					   #end time of the simulation
HAART <- 1996-t0			 #introduction of HAART after start infection
Sim.time <- now-t0		 #simulation time in years
actual.time <- seq(t0,t0+Sim.time,by=0.1)
time.out <- actual.time-t0
###____________________________________________________________________________




### ODE model ###
pop.ODE <- function(t,y,parms)
{  with(as.list(c(y, parms)), 
{	N=HP+HM+HX+P.WT+P.V+M.WT+M.V+X.WT+X.V
  lamda.WT=(t<HAART)*(beta.WT.preHAART/N)*(P.WT+M.WT+X.WT)+
    (t>=HAART)*(beta.WT.postHAART/N)*(P.WT+M.WT+X.WT)
  lamda.V=(t<HAART)*(beta.V.preHAART/N)*(P.V+M.V+X.V)+
    (t>=HAART)*(beta.V.postHAART/N)*(P.V+M.V+X.V)
  
  dHP=b*p*B-(lamda.WT+lamda.V+mu)*HP
  dHM=(1-b)*p*B-(lamda.WT+lamda.V+mu)*HM
  dHX=(1-p)*B-(lamda.WT+lamda.V+mu)*HX
  
  dP.WT=lamda.WT*HP-phi*P.WT-(t<HAART)*alpha.preHAART.PWT*P.WT-
    (t>=HAART)*alpha.postHAART.PWT*P.WT
  dP.V=lamda.V*HP+phi*P.WT-(t<HAART)*alpha.preHAART*P.V-
    (t>=HAART)*alpha.postHAART*P.V
  dM.WT=lamda.WT*HM+psi*M.V-(t<HAART)*alpha.preHAART*M.WT-
    (t>=HAART)*alpha.postHAART*M.WT
  dM.V=lamda.V*HM-psi*M.V-(t<HAART)*alpha.preHAART*M.V-
    (t>=HAART)*alpha.postHAART*M.V
  dX.WT=lamda.WT*HX+psi*X.V-(t<HAART)*alpha.preHAART*X.WT-
    (t>=HAART)*alpha.postHAART*X.WT
  dX.V=lamda.V*HX-psi*X.V-(t<HAART)*alpha.preHAART*X.V-
    (t>=HAART)*alpha.postHAART*X.V
  
  res<-c(dHP,dHM,dHX,dP.WT,dP.V,dM.WT,dM.V,dX.WT,dX.V)
  return(list(res))
})
}




### PLOT fraction V over time ###
phipsi.combi <- matrix(c(0.1,1,0.1,0,0,0.05), nrow=3,ncol=2)
rownames(phipsi.combi) <- c('run1','run2','run3')
colnames(phipsi.combi) <- c('phi','psi')

b <- 0.5

yinit <- c(HP=b*p*N0,HM=(1-b)*p*N0,HX=(1-p)*N0,P.WT=1,P.V=0,M.WT=1,M.V=0,X.WT=1,X.V=0)
names(yinit)<- c('HP','HM','HX','P.WT','P.V','M.WT','M.V','X.WT','X.V')

for(i in 1:dim(phipsi.combi)[1])
{
  P <- c(phi=phipsi.combi[i,'phi'],psi=phipsi.combi[i,'psi'])
  pop.pred <- as.data.frame(ode(func=pop.ODE,y=yinit,parms=P,times=time.out,method=rkMethod("rk45f")))
  freq.V <- rowSums(pop.pred[,c('P.V','M.V')])/rowSums(pop.pred[,c('P.WT','P.V','M.WT','M.V')])
  freq.Vmm <- pop.pred[,c('X.V')]/rowSums(pop.pred[,c('X.WT','X.V')])

  plot.mat <- cbind(time.out+1966,freq.V,freq.Vmm)
  colnames(plot.mat) <- c('time','fV-poss','fV-negs')

  write.table(plot.mat,paste0('C:\\Users\\lies\\Dropbox\\Marjet_Alter\\Death_rate_impacts\\plot_',i,'.txt'),row.names=F)
}

win.graph()
plot(c(1966,2012),c(0,1),col='white',xlim=c(1966,2012),ylim=c(0,1),xlab='time (yrs)',ylab='fraction infected with variant strain',bty='l')
for(i in 1:dim(phipsi.combi)[1])
{
  plot.mat<- read.table(paste('C:\\Users\\lies\\Dropbox\\Marjet_Alter\\Death_rate_impacts\\plot_',i,'.txt',sep=''),header=T)
  lines(plot.mat[,1:2],lty=i,col='black')
  lines(plot.mat[,c(1,3)],lty=i,col='grey')  
}


