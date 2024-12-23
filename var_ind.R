
### please run the 'data.R' file for importing population ### 

# bootstrapping
# probability sample from stratified sampling (on x2) 
# nonprobability sample is selective at random (on x1 and x4)
# P(S*=0|S=1, x) = P(S*=0|x)

start_time = Sys.time()

library(sampling)

set.seed(14)

colnames(pop) = c('x1', 'x2', 'x3', 'x4', 'y') 
# x1 = register year, x2 = engine type, x3 = fuel type, x4 = owner age, y = distance

N = nrow(pop)

scenario = expand.grid(P = c(0.01, 0.1),
                       NP= c(0.05, 0.3, 0.5))
M = 1000; boot = 500; D = 1

i = m = 1

# scenario function
scen.fun <- function(i){
  
  icf.P = scenario[i, 1] ; icf.NP = scenario[i, 2] # inclusion fraction
  bias = mse = var.est = CI = 0

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
  
  
  # formula (18)
  f18 = 1+(W.NP-1)/O 
  x_bar = sum(sample.NP$y*f18)/sum(f18)
  mu = mean(pop$y)
  
  bias = (x_bar - mu)/M + bias
  mse = (x_bar - mu)^2/M + mse
  
 
  # estimation function
  estfun <- function(bsample.p, bsample.np){
    
    y = bsample.np$y
    bsample.np$Z = 1
    bsample.p$Z = 0
    indP <- row.names(bsample.p)
    indNP <- row.names(bsample.np)
    A = intersect(indP, indNP)
    B = rbind(bsample.p[!c(indP %in% A), ], bsample.np[!c(indNP %in% A), ])
    
    glmO = glm(Z ~ x1+x2+x4, data = B, family = "binomial")
    O = exp(predict(glmO, newdata = bsample.np))
    
    W.NP = 1/bsample.np$includP
    bf18 = 1+(W.NP-1)/O
    p18 = sum(y*bf18)/sum(bf18)
    return(p18)
  }
  
  # choose one of the two following lines as a basis for bootstrapping
  f.norm <- N * f18 / sum(f18)
  
  v = 0

  for (d in 1:D) {
  
  ## finite population bootstrap based on NP
  f.norm.rounded <- trunc(f.norm) + UPrandomsystematic(f.norm - trunc(f.norm)) # round weights using the method of Fellegi (1975)
  inflate.NP <- rep(1:nrow(sample.NP), times = f.norm.rounded)
  fakepop.NP <- sample.NP[inflate.NP, ]
  
  fakepop.NP$W.NP <- W.NP[inflate.NP]
  fakepop.NP$includP <- inclusionprobabilities(1 / fakepop.NP$W.NP, n = nrow(sample.P))
  fakepop.NP$includNP <- inclusionprobabilities((1 / f.norm)[inflate.NP], n = nrow(sample.NP))
  eps = 1e-6

  bsample.p <- replicate(boot, fakepop.NP[as.logical(UPrandomsystematic(fakepop.NP$includP, eps = eps/10)), ], simplify = F)
  bsample.np <- replicate(boot, fakepop.NP[as.logical(UPrandomsystematic(fakepop.NP$includNP, eps = eps/10)), ], simplify = F)
    
  estPNP = mapply(estfun, bsample.p = bsample.p, bsample.np = bsample.np)
  
  v[d] = var(estPNP)
  
  }
  
  var.est[m] = mean(v)
  
  CI = ( ((x_bar - 1.96*sqrt(v)) < mu) & (mu < (x_bar + 1.96*sqrt(v))) )/M + CI
}

  return(list(icf.NP, icf.P, var.est, mse-bias^2, CI))
  
}

# parallel 

library(doRNG); library(doParallel)

registerDoParallel(cores = min(c(nrow(scenario), (detectCores()-1))))
registerDoRNG(13)

result <- foreach(i = 1:nrow(scenario),
                  .packages = c('sampling'),
                  .combine = rbind) %dopar% scen.fun(i)

table = matrix(data = mapply(mean, result), nrow = nrow(scenario))
colnames(result) = c('n.NP', 'n.P', 'est', 'sim.var', 'CI')
colnames(table) = c('n.NP', 'n.P', 'est', 'sim.var', 'CI')

capture.output(round(table, 4), file = "//cbsp.nl/Productie/secundair/MPOnderzoek/Werk/Combineren/Projecten/nonprobability/variance/var_ind.txt")

total.time = Sys.time() - start_time

save.image("//cbsp.nl/Productie/secundair/MPOnderzoek/Werk/Combineren/Projecten/nonprobability/variance/var_ind.RData")



library("xtable")
table = as.data.frame(table)
table$RB = 100*(table$est-table$sim.var)/table$sim.var
table$CI = 100*table$CI
xtable(table[,c(1,2,6,5)], digits = c(0,2,2,2,2), caption = 'Variance Estimates of Ind')

