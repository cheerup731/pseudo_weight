
### please run the 'data.R' file for importing population ### 

# convergence test of the true variance  
# probability sample from stratified sampling (on x2) 
# nonprobability sample is selective at random (on x1 and x4)
# P(S*=0|S=1, x) != P(S*=0|x)
start_time = Sys.time()

library(sampling)

set.seed(14)

colnames(pop) = c('x1', 'x2', 'x3', 'x4', 'y') 
# x1 = register year, x2 = engine type, x3 = fuel type, x4 = owner age, y = distance

N = nrow(pop)
# pop = pop[order(pop$x2),] # for stratified sampling

scenario = expand.grid(P = c(0.01, 0.1),
                       NP= c(0.05, 0.3, 0.5))

M = 2000

i = m = 11

scen.fun <- function(i){

  icf.P = scenario[i, 1] ; icf.NP = scenario[i, 2] # inclusion fraction
  bias = mse = var.est = 0

for (m in 1:M) {
  
  # probability sample 
  x = pop$x1/30
  f = function(theta) sum(exp(theta + x) / (1 + exp(theta + x))) - icf.P*N
  theta = uniroot(f, c(-100, 100))$root
  includP = exp(theta + x) / (1 + exp(theta + x))
  S = as.logical(UPrandomsystematic(includP))
  sample.P = pop[S,-5]
  
  # nonprobability sample 
  x = pop$x2 - pop$x4/20 
  x[S] = x[S]-2
  f = function(theta) sum(exp(theta + x) / (1 + exp(theta + x))) - icf.NP*N
  theta = uniroot(f, c(-100, 100))$root
  includNP = exp(theta + x) / (1 + exp(theta + x))
  Sstar = as.logical(UPrandomsystematic(includNP))
  sample.NP = pop[Sstar,]
  
  # design weight
  W.P  = 1/includP[S]
  W.NP = 1/includP[Sstar]
  
  # B set
  B = pop[Sstar+S == 1, ]
  B$Z = Sstar[Sstar+S == 1]
  
  # estimate O
  glmO = glm(Z ~ x1+x2+x4, data = B, family = "binomial")
  O = exp(predict(glmO, newdata = sample.NP))
  
  # estimate L
  sample.P$Sstar = Sstar[S == 1]
  glmL = glm(Sstar ~ x1+x2+x4, data = sample.P, family = "binomial")
  L = 1 - predict(glmL, newdata = sample.NP, type = "response")
  
  # formula (11)
  f11 = (W.NP-1)/(O*L) 
  
  bias = (c(sum(sample.NP$y*f11)/sum(f11))- mean(pop$y))/M + bias
  mse = (c(sum(sample.NP$y*f11)/sum(f11)) - mean(pop$y))^2/M + mse
           
  
}

  return(c(icf.NP, icf.P, mse-bias^2))
 
}


                     
library(doRNG); library(doParallel)

registerDoParallel(cores = min(c(nrow(scenario), (detectCores()-1))))
registerDoRNG(13)

result <- foreach(i = 1:nrow(scenario),
                  .packages = c('sampling'),
                  .combine = rbind) %dopar% scen.fun(i)

colnames(result) = c('n.NP', 'n.P', 'true_var_dep')

capture.output(round(result, 4), file = "//cbsp.nl/Productie/secundair/MPOnderzoek/Werk/Combineren/Projecten/nonprobability/variance/true_var_dep.txt")

total.time = Sys.time() - start_time

save.image("//cbsp.nl/Productie/secundair/MPOnderzoek/Werk/Combineren/Projecten/nonprobability/variance/true_var_dep.RData")

