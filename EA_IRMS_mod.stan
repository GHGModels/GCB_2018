data{

int<lower = 0> N; //#observations
int<lower = 0> K; // #pastures
real y[N];//insert response here 
int<lower = 0> G[N];// indicator for grazing 
int<lower = 0> Back[N];// indicator for yes/no background sample
int<lower = 0> Past[N]; //index for pasture 

//passed in priors 
real<lower = 0> sigma_scales; 
real<lower = 0> sigma_location; 

}

parameters{

real alpha;
real beta; 
real beta_back; // background sample? 
 
real<lower = 0> sigma_pas; 
real B[3]; // treatment coefficients
real Pas[K]; // plot coefficients

}


model {

// treatment effect priors 
B[1] ~ normal(0, 2*sigma_location);//effectively flat prior to identify location 
B[2] ~ normal(0, 5*sigma_scales); // weakly informative 
B[3] ~ normal(0, 5*sigma_scales); // weakly informative 

// priors for variance model 
alpha ~ normal(0,5);
beta ~ normal(0,5); 
beta_back ~ normal(0,5); 

// pasture hierarchical variance term 
sigma_pas ~ student_t(4,0, sigma_scales); 

// pasture random effects 
for(k in 1:3)
Pas[k] ~ student_t(4,0,sigma_pas); 


// sampling distribution with model for both mean and variance

for(i in 1:N)
y[i] ~ normal(B[1] + B[2]*G[i] + B[3]*Back[i] + Pas[Past[i]],
exp(alpha + beta*G[i] + beta_back*Back[i]));

} 

generated quantities{

vector[N] preds; 
real mean_G;
real mean_EXC; 

for(i in 1:N){
preds[i] <- normal_rng(B[1] + B[2]*G[i] + B[3]*Back[i] + Pas[Past[i]],
exp(alpha + beta*G[i] + beta_back*Back[i]));
}

mean_G <- B[1];
mean_EXC <- B[1] + B[2];

}
