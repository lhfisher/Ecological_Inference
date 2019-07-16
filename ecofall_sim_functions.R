
##----------------------------------------------------------------------------##
#### Functions for simulations in the absence of vaccinnation (Section 4.1) ####
##----------------------------------------------------------------------------##

##########################################################
#### Function to simulate one unvaccinated population
##########################################################
sim.pop=function(N, Nwks, AR, EN, y0=1, burnin=4, report_all=F){
        lt=Yt=St=rep(NA, Nwks+burnin)
        Yt[1]=y0
        St[1]=N-y0

        for(t in 2:(Nwks+burnin)){
                if(length(EN)==1){
                        lt[t] = exp(AR)*Yt[t-1]/N + exp(EN)
                }else{
                        lt[t] = exp(AR)*Yt[t-1]/N + exp(EN[t])
                }
                
                Yt[t]=rbinom(1, size=St[t-1], prob=1-exp(-lt[t]))
                
                St[t] = St[t-1]-Yt[t]
                
        }
        if(report_all){
                return(list(S=St, Y=Yt, lambda=lt))
        }else{
         return(list(S=St[(burnin+1):(Nwks+burnin)], Y=Yt[(burnin+1):(Nwks+burnin)], lambda=lt[(burnin+1):(Nwks+burnin)]))
        }
        
        
}



##########################################################
#### Likelihood functions in the absence of vaccination
##########################################################
## These are the 8 log-likelihoods fit to simulated data in the absence of vaccination 

# Model 1:  Y_{t+1}|lambda_t ~ Bin(S_t, 1-exp(-lambda_t))
loglik1=function(param, y, s, n){
        ar=exp(param[1])
        en=exp(param[2])
        l=0
        nwks=length(y)
        for(t in 2:nwks){
                lam=ar*y[t-1]/n + en
                l=l-dbinom(y[t], size=s[t-1], prob=1-exp(-lam), log=T)
        }
        return(l)
}

gr1=function(param, y, s, n){
        ar=exp(param[1])
        en=exp(param[2])
        d.ar=d.en=0
        nwks=length(y)
        for(t in 2:nwks){
                lam=ar*y[t-1]/n + en
                d.ar=d.ar+(y[t]*exp(lam)/(1-exp(-lam)) - (s[t-1]-y[t]))*y[t-1]/n
                d.en=d.en+(y[t]*exp(lam)/(1-exp(-lam)) - (s[t-1]-y[t]))
        }
        return(-c(d.ar*ar, d.en*en))
}

# Model 2:  Y_{t+1}|lambda_t ~ Bin(N, 1-exp(-lambda_t))
loglik2=function(param, y, n){
        ar=exp(param[1])
        en=exp(param[2])
        l=0
        nwks=length(y)
        for(t in 2:nwks){
                lam=ar*y[t-1]/n + en
                l=l-dbinom(y[t], size=n, prob=1-exp(-lam), log=T)
        }
        return(l)
}

# Model 3:  Y_{t+1}|lambda_t ~ Bin(S_t, lambda_t)
loglik3=function(param, y, s, n){
        ar=exp(param[1])
        en=exp(param[2])
        l=0
        nwks=length(y)
        for(t in 2:nwks){
                lam=ar*y[t-1]/n + en
                l=l-dbinom(y[t], size=s[t-1], prob=lam, log=T)
        }
        return(l)
}

# Model 4:  Y_{t+1}|lambda_t ~ Bin(N, lambda_t)
loglik4=function(param, y, n){
        ar=exp(param[1])
        en=exp(param[2])
        l=0
        nwks=length(y)
        for(t in 2:nwks){
                lam=ar*y[t-1]/n + en
                l=l-dbinom(y[t], size=n, prob=lam, log=T)
        }
        return(l)
}

# Model 5:  Y_{t+1}|lambda_t ~ Poi(S_t (1-exp(-lambda_t)))
loglik5=function(param, y, s, n){
        ar=exp(param[1])
        en=exp(param[2])
        l=0
        nwks=length(y)
        for(t in 2:nwks){
                lam=ar*y[t-1]/n + en
                l=l-dpois(y[t], lambda=s[t-1]*(1-exp(-lam)), log=T)
        }
        return(l)
}

# Model 6:  Y_{t+1}|lambda_t ~ Poi(N (1-exp(-lambda_t)))
loglik6=function(param, y, n){
        ar=exp(param[1])
        en=exp(param[2])
        l=0
        nwks=length(y)
        for(t in 2:nwks){
                lam=ar*y[t-1]/n + en
                l=l-dpois(y[t], lambda=n*(1-exp(-lam)), log=T)
        }
        return(l)
}

# Model 7:  Y_{t+1}|lambda_t ~ Poi(S_t lambda_t)
loglik7=function(param, y, s, n){
        ar=exp(param[1])
        en=exp(param[2])
        l=0
        nwks=length(y)
        for(t in 2:nwks){
                lam=ar*y[t-1]/n + en
                l=l-dpois(y[t], lambda=s[t-1]*lam, log=T)
        }
        return(l)
}

# Model 8:  Y_{t+1}|lambda_t ~ Poi(N lambda_t)
loglik8=function(param, y, n){
        ar=exp(param[1])
        en=exp(param[2])
        l=0
        nwks=length(y)
        for(t in 2:nwks){
                lam=ar*y[t-1]/n + en
                l=l-dpois(y[t], lambda=n*lam, log=T)
        }
        return(l)
}

###########################################
#### Function to simulate one totally unvaccinated population, and fit the 8 models to the data to assess simplifying assumptions in the ecological vaccine model
###########################################



fit.mod=function(Func, inits1=c(1, -5), inits2=c(0, -12), Gr=NULL, ...){
        tmp=optim(par=inits1, fn=Func, hessian=T, ...)
        
        ## check to see if there are 0's on the diag in the hessian
        if(sum(diag(tmp$hessian)==0)>0){
                mod=optim(par=inits1, fn=Func, hessian=T, method='BFGS', ...) # , control=list(reltol=1e-9)
        }else{
                mod = tmp}
        
        # mod$hessian = hessian(func=Func, x=mod$par, ...)
        
        if(sum(diag(solve(mod$hessian))<0)>0){
                mod2=optim(par=inits2, fn=Func, hessian=T, ...)
        }else{
                mod2=mod
        }
        
        
        
        return(mod2)
}


sim.simple=function(N, Nwks, alpha_ar, alpha_en, ...){
        obs=sim.pop(N=N, Nwks=Nwks, AR=alpha_ar, EN=alpha_en, ...)
        # fit each
        mod1=fit.mod(Func=loglik1, y=obs$Y, s=obs$S, n=N)
        
        mod2=fit.mod(Func=loglik2, y=obs$Y, n=N)
        
        mod3=fit.mod(Func=loglik3, y=obs$Y, s=obs$S, n=N)
        
        mod4=fit.mod(Func=loglik4, y=obs$Y, n=N)
        
        mod5=fit.mod(Func=loglik5, y=obs$Y, s=obs$S, n=N)
        
        mod6=fit.mod(Func=loglik6, y=obs$Y, n=N)
        
        mod7=fit.mod(Func=loglik7, y=obs$Y, s=obs$S, n=N)
        
        mod8=fit.mod(Func=loglik8, y=obs$Y, n=N)
        
        res=list(dat=obs, mod1=mod1, mod2=mod2, mod3=mod3, mod4=mod4, mod5=mod5, mod6=mod6, mod7=mod7, mod8=mod8)
        return(res)
        
}


##---------------------------------------------------------------------------------##
#### Functions for simulations in partially vaccinated populations (Section 4.2) ####
##---------------------------------------------------------------------------------##

#### functions to simulate one all-or-none vaccine population
sim.AoN.pop=function(Npop, nwks, AR, EN, Phi, X, y0=c(0, 0), Ntot, burn=4){
        lt=Yvtot=Yutot=Sutot=Svtot=Stot=Ytot=rep(NA, nwks+burn)
        Yutot[1]=y0[1]
        Yvtot[1]=y0[2]
        Ytot[1]=Yutot[1]+Yvtot[1]
        
        Sutot[1]=as.integer((1-X)*Npop -y0[1])
        Svtot[1]=as.integer((1-Phi)*X*Npop - y0[2])
        Stot[1]=Sutot[1]+Svtot[1]
        for(t in 2:(nwks+burn)){
                lt[t]=exp(AR)*Ytot[t-1]/Npop + exp(EN-log(Ntot))
                
                Yutot[t]=rbinom(1, size=Sutot[t-1], prob=1-exp(-lt[t]))
                Yvtot[t]=rbinom(1, size=Svtot[t-1], prob=1-exp(-lt[t]))
                
                up_su=Sutot[t-1]-Yutot[t]
                up_sv=Svtot[t-1]-Yvtot[t]
                
                Svtot[t] = ifelse(up_sv<0, 0, up_sv)
                Sutot[t] = ifelse(up_su<0, 0, up_su)
                Stot[t] = Sutot[t] + Svtot[t]
                Ytot[t] = Yutot[t] + Yvtot[t]
                
                Stot[t]=Stot[t-1]-Ytot[t]
        }
        return(list(Y=Ytot[(burn+1):(nwks+burn)], S=Stot[(burn+1):(nwks+burn)], Yu=Yutot[(burn+1):(nwks+burn)], Yv=Yvtot[(burn+1):(nwks+burn)], Su=Sutot[(burn+1):(nwks+burn)], Sv=Svtot[(burn+1):(nwks+burn)], lu=lt[(burn+1):(nwks+burn)], lv=lt[(burn+1):(nwks+burn)] ))
}

## simulate multiple areas at once
sim.AoN.pops=function(Npops, nwks, AR, EN, Phi, X, ...){
        lumat=lvmat=svmat=yvmat=sumat=yumat=smat=ymat=matrix(NA, nrow=length(Npops), ncol=nwks)
        for(i in 1:length(Npops)){
                tmp=sim.AoN.pop(Npop=Npops[i], nwks=nwks, AR=AR, EN=EN, Phi=Phi, X=X[i], Ntot=sum(Npops), ...)
                ymat[i, ]=tmp$Y
                smat[i, ]=tmp$S
                yumat[i, ]=tmp$Yu
                yvmat[i, ]=tmp$Yv
                sumat[i, ]=tmp$Su
                svmat[i, ]=tmp$Sv
                lumat[i, ]=tmp$lu
                lvmat[i, ]=tmp$lv
        }
        return(list(Y=ymat, S=smat, Yu=yumat, Yv=yvmat, Su=sumat, Sv=svmat, lu=lumat, lv=lvmat, params=list(alpha_ar=AR, alpha_en=EN, phi=Phi, N=Npops, Nwks=nwks, x=X, vacc='All-or-none')))
}



#### function to simulate leaky vaccine population ####
sim.lk.pop=function(Npop, nwks, AR, EN, Phi, X, y0=c(0, 0), Ntot){
        lvt=lut=Yvtot=Yutot=Sutot=Svtot=Stot=Ytot=rep(NA, nwks+4)
        Yutot[1]=y0[1]
        Yvtot[1]=y0[2]
        Ytot[1]=Yutot[1]+Yvtot[1]
        
        Sutot[1]=as.integer((1-X)*Npop-y0[1])
        Svtot[1]=as.integer(X*Npop - y0[2])
        Stot[1]=Sutot[1]+Svtot[1]
        for(t in 2:(nwks+4)){
                lut[t]=exp(AR)*Ytot[t-1]/Npop + exp(EN-log(Ntot))
                lvt[t]=(1-Phi)*lut[t]
                
                Yutot[t]=rbinom(1, size=Sutot[t-1], prob=1-exp(-lut[t]))
                Yvtot[t]=rbinom(1, size=Svtot[t-1], prob=1-exp(-lvt[t]))
                
                up_su=Sutot[t-1]-Yutot[t]
                up_sv=Svtot[t-1]-Yvtot[t]
                
                Svtot[t] = ifelse(up_sv<0,0,  up_sv)
                Sutot[t] = ifelse(up_su<0, 0,  up_su)
                Stot[t] = Sutot[t] + Svtot[t]
                Ytot[t] = Yutot[t] + Yvtot[t]
                
                Stot[t]=Stot[t-1]-Ytot[t]
        }
        return(list(Y=Ytot[5:(nwks+4)], S=Stot[5:(nwks+4)], Yu=Yutot[5:(nwks+4)], Yv=Yvtot[5:(nwks+4)], Su=Sutot[5:(nwks+4)], Sv=Svtot[5:(nwks+4)], lu=lut[5:(nwks+4)], lv=lvt[5:(nwks+4)] ))
}

## simulate multiple areas at once
sim.lk.pops=function(Npops, nwks, AR, EN, Phi, X, ...){
        lumat=lvmat=svmat=yvmat=sumat=yumat=smat=ymat=matrix(NA, nrow=length(Npops), ncol=nwks)
        for(i in 1:length(Npops)){
                tmp=sim.lk.pop(Npop=Npops[i], nwks=nwks, AR=AR, EN=EN, Phi=Phi, X=X[i], Ntot=sum(Npops), ...)
                ymat[i, ]=tmp$Y
                smat[i, ]=tmp$S
                yumat[i, ]=tmp$Yu
                yvmat[i, ]=tmp$Yv
                sumat[i, ]=tmp$Su
                svmat[i, ]=tmp$Sv
                lumat[i, ]=tmp$lu
                lvmat[i, ]=tmp$lv
        }
        return(list(Y=ymat, S=smat, Yu=yumat, Yv=yvmat, Su=sumat, Sv=svmat, lu=lumat, lv=lvmat, 
                    params=list(alpha_ar=AR, alpha_en=EN, phi=Phi, N=Npops, Nwks=nwks, x=X, vacc='Leaky')))
}



##################################################
#### Likelihood functions 
#### for partially vaccinated populations
##################################################
## These are the 4 log-likelihoods fit to simulated data in partially vaccinated populations

## Model 1: Completely observed binomial, leaky
loglik_bin_lk=function(param, yu, yv, xi, n){
        ar=exp(param[1])
        en=exp(param[2])
        v=expit(param[3])
        l=0
        if(length(xi)==1){
                nwks=length(yu)
                for(t in 2:nwks){
                        su=(1-xi)*n - sum(yu[1:(t-1)])
                        sv=xi*n - sum(yv[1:(t-1)])
                        lam_u=ar*(yu[t-1]+yv[t-1])/n + en/n
                        lam_v=(ar*(yu[t-1]+yv[t-1])/n + en/n)*(1-v)
                        l=l-dbinom(yu[t], size=su, prob=1-exp(-lam_u), log=T) -dbinom(yv[t], size=sv, prob=1-exp(-lam_v), log=T)
                }
        }else{
                nwks=dim(yu)[2]
                ni=dim(yu)[1]
                for(i in 1:ni){
                        for(t in 2:nwks){
                                su=(1-xi[i])*n[i] - sum(yu[i, 1:(t-1)])
                                sv=xi[i]*n[i] - sum(yv[i, 1:(t-1)])
                                lam_u=ar*(yu[i, t-1]+yv[i, t-1])/n[i] + en/sum(n)
                                lam_v=(ar*(yu[i, t-1]+yv[i, t-1])/n[i] + en/sum(n))*(1-v)
                                l=l-dbinom(yu[i, t], size=su, prob=1-exp(-lam_u), log=T) -dbinom(yv[i, t], size=sv, prob=1-exp(-lam_v), log=T)
                        }
                }
        }
        return(l)
}

## Model 2: Completely observed binomial all-or-none
loglik_bin_aon=function(param, yu, yv, xi, n){
        ar=exp(param[1])
        en=exp(param[2])
        v=expit(param[3])
        l=0
        nwks=dim(yu)[2]
        ni=dim(yu)[1]
        for(i in 1:ni){
                for(t in 2:nwks){
                        # st=n[i]*(1-v*xi[i])-sum(y[i, 1:(t-1)]) # S_{i, t-1}
                        # lam=ar*y[i, t-1]/n[i] + en/n[i]
                        # l=l-dbinom(y[i, t], size=as.integer(st), prob=1-exp(-lam), log=T)
                        su=(1-xi[i])*n[i] - sum(yu[i, 1:(t-1)])
                        sv=(1-v)*xi[i]*n[i] - sum(yv[i, 1:(t-1)])
                        su=ifelse(su<0, 0, su)
                        sv=ifelse(sv<0, 0, sv)
                        lam=ar*(yu[i, t-1]+yv[i, t-1])/n[i] + en/sum(n)
                        l=l-dbinom(yu[i, t], size=as.integer(su), prob=1-exp(-lam), log=T) - dbinom(yv[i, t], size=as.integer(sv), prob=1-exp(-lam), log=T)
                }
        }
        return(l)
}


## Model 3: Epidemic-endemic from Herzog et al. 2011
loglik_herzog=function(param, y, xi, n){
        ar=exp(param[1])
        en=exp(param[2])
        b=param[3]
        
        nwks=dim(y)[2]
        ni=dim(y)[1]
        l=0
        for(i in 1:ni){
                for(t in 2:nwks){
                        mu=(ar*(1-xi[i])^b)*y[i, t-1] + en*(n[i]/sum(n))
                        l=l-dpois(y[i, t], lambda=mu, log=T)
                }
        }
        return(l)
}


## Model 4: Ecological vaccine model
loglik_agg=function(param, y, xi, n){
        ar=exp(param[1])
        en=exp(param[2])
        v=expit(param[3])
        l=0
        nwks=dim(y)[2]
        ni=dim(y)[1]
        for(i in 1:ni){
                for(t in 2:nwks){
                        
                        lam=(ar*y[i, t-1]/n[i] + en/sum(n) )
                        l=l-dpois(y[i, t], lambda=n[i]*(1-v*xi[i])*lam, log=T)
                }
        }
        return(l)
}


############################################
### function to simulate one population, 
###     and then fit all the estimates
############################################

sim.one=function(vacc='AoN', N, Nwks, alpha_ar, alpha_en, phi, x, sts=c(0.5, 0, 1)){
        # Simulate data either all-or-none, or leaky
        if(vacc=='AoN'){
                obs=sim.AoN.pops(Npops=N, nwks=Nwks, AR=alpha_ar, EN=alpha_en, Phi=phi, X=x, y0=c(1, 0)) 
        }else{
                obs=sim.lk.pops(Npops=N, nwks=Nwks, AR=alpha_ar, EN=alpha_en, Phi=phi, X=x, y0=c(1, 0))
        }
        
        # Model 1: Fully-observed all-or-none model
        mod1=optim(par=sts, fn=loglik_bin_aon, yu=obs$Yu, yv=obs$Yv, xi=x, n=N, hessian=T, method='BFGS')

        
        # Model 2: Fully-observed leaky model
        mod2=optim(par=sts, fn=loglik_bin_lk, yu=obs$Yu, yv=obs$Yv, xi=x, n=N, hessian=T, method='BFGS')
        
        # Model 3: ecological model
        mod3=optim(par=sts, fn=loglik_agg, y=obs$Y, xi=x, n=N, hessian=T, method='BFGS')
        
        
        # Model 4: Epidemic-endemic model
        mod4=optim(par=sts, fn=loglik_herzog, y=obs$Y, xi=x, n=N, hessian=T, method='BFGS')
        
        
        ## return everything
        res=list(dat=obs, mod1=mod1, mod2=mod2, mod3=mod3, mod4=mod4)
        return(res)
}


###################################
## Code to demonstrate ecologial bias in the binary outcome setting
#################################

##----------------------------------------------------------------------------------------##
#### Functions for simulations for ecological bias from a binary covariate (Appendix A) ####
##----------------------------------------------------------------------------------------##


#### function to simulate simple population, to demonstrate ecological bias ####
sim.bias.pop=function(Npop, nwks, a0, a1, a2, X, y0=c(1, 0), Ntot){
        lvt=lut=Yvtot=Yutot=Sutot=Svtot=Stot=Ytot=rep(NA, nwks+4)
        Yutot[1]=y0[1]
        Yvtot[1]=y0[2]
        Ytot[1]=Yutot[1]+Yvtot[1]
        
        Sutot[1]=as.integer((1-X)*Npop-y0[1])
        Svtot[1]=as.integer(X*Npop - y0[2])
        Stot[1]=Sutot[1]+Svtot[1]
        for(t in 2:(nwks+4)){
                lut[t]=exp(a0)*Ytot[t-1]/Npop + exp(a2)
                lvt[t]=exp(a0 + a1)*Ytot[t-1]/Npop + exp(a2)
                
                Yutot[t]=rbinom(1, size=Sutot[t-1], prob=lut[t]) # 1-exp(-lut[t]))
                Yvtot[t]=rbinom(1, size=Svtot[t-1], prob=lvt[t]) # 1-exp(-lvt[t]))
                
                Sutot[t]=Sutot[t-1]-Yutot[t]
                Svtot[t]=Svtot[t-1]-Yvtot[t]
                
                Stot[t] = Sutot[t] + Svtot[t]
                Ytot[t] = Yutot[t] + Yvtot[t]
                
                Stot[t]=Stot[t-1]-Ytot[t]
        }
        return(list(Y=Ytot[5:(nwks+4)], S=Stot[5:(nwks+4)], Yu=Yutot[5:(nwks+4)], Yv=Yvtot[5:(nwks+4)], Su=Sutot[5:(nwks+4)], Sv=Svtot[5:(nwks+4)], lu=lut[5:(nwks+4)], lv=lvt[5:(nwks+4)] ))
}


## simulate multiple areas at once
sim.bias.pops=function(Npops, nwks, a0, a1, a2, X, ...){
        lumat=lvmat=svmat=yvmat=sumat=yumat=smat=ymat=matrix(NA, nrow=length(Npops), ncol=nwks)
        for(i in 1:length(Npops)){
                tmp=sim.bias.pop(Npop=Npops[i], nwks=nwks, a0=a0, a1=a1, a2=a2, X=X[i], Ntot=sum(Npops), ...)
                ymat[i, ]=tmp$Y
                smat[i, ]=tmp$S
                yumat[i, ]=tmp$Yu
                yvmat[i, ]=tmp$Yv
                sumat[i, ]=tmp$Su
                svmat[i, ]=tmp$Sv
                lumat[i, ]=tmp$lu
                lvmat[i, ]=tmp$lv
        }
        return(list(Y=ymat, S=smat, Yu=yumat, Yv=yvmat, Su=sumat, Sv=svmat, lu=lumat, lv=lvmat, 
                    params=list(a0=a0, a1=a1, a2=a2, N=Npops, Nwks=nwks, x=X)))
}



## Naive ecological regression (A6)
loglik_naive = function(param, y, xi, n){
        l=0
        nwks=dim(y)[2]
        ni=dim(y)[1]
        for(i in 1:ni){
                for(t in 2:nwks){
                        lam=exp(param[1] + param[2]*xi[i])*y[i, t-1] + exp(param[3])*n[i]
                        l = l-dpois(y[i, t], lambda = lam, log=T)
                }
        }
        return(l)
}

## Implied aggregate risk model (A5)
loglik_implied_agg = function(param, y, xi, n){
        l=0
        nwks=dim(y)[2]
        ni=dim(y)[1]
        for(i in 1:ni){
                for(t in 2:nwks){
                        lam=exp(param[1])*((1-xi[i]) + xi[i]*exp(param[2]))*y[i, t-1] + exp(param[3])*n[i]
                        l = l-dpois(y[i, t], lambda = lam, log=T)
                }
        }
        return(l)
}



############################################
### function to simulate one population, 
###     and then fit all the implied aggregate
###     and naive models
############################################
sim.eco.bias = function(N, Nwks, x, a0, a1, a2, sts=c(1, 1, -1)){
        dt = sim.bias.pops(Npop=N, nwks=Nwks, a0=a0, a1=a1, a2=a2, X=x, y0=c(1, 0))
        
        agg=optim(par=sts, fn=loglik_implied_agg, y=dt$Y, xi=x, n=N)
        # agg=optim(par=agg1$par, fn=loglik_implied_agg, y=dt$Y, xi=x, n=N)
        nai=optim(par=sts, fn=loglik_naive, y=dt$Y, xi=x, n=N)
        # nai=optim(par=nai1$par, fn=loglik_naive, y=dt$Y, xi=x, n=N)
        
        res=list(dat=dt, agg=agg, naive=nai)
        return(res)
}

