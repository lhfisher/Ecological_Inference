
mod_ENre='
data{
int<lower=1> I;          // number of areas
int<lower=1> Ts;         // number of times
int<lower=1> Nobs;    // number of total observations
int y[I, Ts] ;  // IxTs matrix of observations
real logNi[I];          // log population proportions 
vector[I] X;            // vector of Xi vacc cover
vector[Ts] sin1;           // vector of sin
vector[Ts] cos1;           // vector of cos
vector[2] a; //vector of parameters for the beta prior on phi
}
parameters{
real alpha_ar;  // AR component - intercept
real alpha_en;  // EN component - intercept
real gamma;     // sin component
real delta;     // cosine component 
real<lower=0, upper=1> phi;  // vaccine effect
vector[I] re_en; // EN random effects
vector[I] re_ar; // AR random effects
real<lower=0> sigma_ar; 
real<lower=0> sigma_en; 
}
model{
for(i in 1:I){
re_ar[i]~normal(0, sigma_ar);
re_en[i]~normal(0, sigma_en);
for(t in 2:Ts){
real muit;
muit = exp(alpha_ar+re_ar[i])*y[i, t-1] + exp(logNi[i] + alpha_en + re_en[i] + gamma*sin1[t] + delta*cos1[t]);
y[i, t] ~ poisson((1-phi*X[i])*muit);
}
}

//priors
alpha_ar ~ normal(0, 5);
alpha_en ~ normal(0, 5);
gamma ~ normal(0, 10);
delta ~ normal(0, 10);
phi ~ beta(a[1], a[2]);

}
generated quantities{
vector[Nobs-I] log_lik;
vector[I] ar;
vector[I] en;
vector[I] r; 
for(i in 1:I){
r[i] = exp(alpha_ar + re_ar[i])*(1-phi*X[i]); 
ar[i] = exp(alpha_ar + re_ar[i]);
en[i] = exp(alpha_en + re_en[i]);
for(t in 2:Ts){
real muit;
muit=exp(alpha_ar + re_ar[i])*y[i, t-1] + exp(logNi[i] + alpha_en + re_en[i] + gamma*sin1[t] + delta*cos1[t]);
log_lik[(i-1)*(Ts-1)+t-1 ] = poisson_lpmf(y[i, t] | (1-phi*X[i])*muit);
}
}
}
'

### fixed effects model - closer to the MLE 
mod_fe='
data{
int<lower=1> I;          // number of areas
int<lower=1> Ts;         // number of times
int<lower=1> Nobs;    // number of total observations
int y[I, Ts] ;  // IxTs matrix of observations
real logNi[I];          // log population sizes
vector[I] X;            // vector of Xi vacc cover
vector[Ts] sin1;           // vector of sin
vector[Ts] cos1;           // vector of cos
vector[2] a; //vector of parameters for the beta prior on phi
}
parameters{
real alpha_ar;  // AR component - intercept
real alpha_en;  // EN component - intercept
real gamma;     // sin component
real delta;     // cosine component 
real<lower=0, upper=1> phi;  // vaccine effect
}
model{
for(i in 1:I){
for(t in 2:Ts){
real muit;
muit = exp(alpha_ar)*y[i, t-1] + exp(logNi[i] + alpha_en + gamma*sin1[t] + delta*cos1[t]);
y[i, t] ~ poisson((1-phi*X[i])*muit);
}
}

//priors
alpha_ar ~ normal(0, 5);
alpha_en ~ normal(0, 5);
gamma ~ normal(0, 10);
delta ~ normal(0, 10);
phi ~ beta(a[1], a[2]);

}
generated quantities{
vector[Nobs-I] log_lik;
vector[I] ri; 
for(i in 1:I){
ri[i] = exp(alpha_ar)*(1-phi*X[i]); 
for(t in 2:Ts){
real muit;
muit=exp(alpha_ar)*y[i, t-1] + exp(logNi[i] + alpha_en + gamma*sin1[t] + delta*cos1[t]);
log_lik[(i-1)*(Ts-1)+t-1 ] = poisson_lpmf(y[i, t] | (1-phi*X[i])*muit);
}
}
}
'