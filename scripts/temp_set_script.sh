#!/bin/bash
# script to optimize the temperature schedule (to get nearly constant acceptance rate)
# based on the energy measurements from a short run.
	/usr/bin/R --slave --args $@ << "EOF"
args <- commandArgs(trailingOnly=T)
	d <- read.table(args[1],header=T);
	s <- smooth.spline(d$T,d$e);

nT <- length(d$T)
	if (length(args)>=2) {
		nT <-as.numeric(args[2])
	}
rT <- range(d$T)
	if (length(args)>=4) {
		rT <-as.numeric(args[3:4])
	}

Tset <- function(step_){
	cost <- function(T_) {
		dE <- predict(s,T)$y-predict(s,T_)$y
			dB <- 1./T-1./T_
			return(dE*dB-step_)
	}
	T <- min(rT)
		Tset_ <- 1:nT
		Tset_[1] <- T
		for(i in 2:nT){
			Tset_[i] <- T <- uniroot(cost,c(T,10*T),tol=1e-15)$root
		}
	return(Tset_)
}

f <- function(step) max(rT)-max(Tset(step))
	step <- uniroot(f,c(-1,-1e-15),tol=1e-15)$root
	write.table(t(sprintf("%.5f",Tset(step))),quote=F,col.names=F,row.names=F,sep="\n")

EOF
