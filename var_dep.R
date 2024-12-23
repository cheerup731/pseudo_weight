
### please run the 'data.R' file for importing population ### 

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
M = 1000; boot = 500; D = 1

i = m = 11

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
  x_bar = sum(sample.NP$y*f11)/sum(f11)
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
    
    bsample.p$Sstar = 0
    bsample.p$Sstar[indP %in% A] = 1
    glmL = glm(Sstar ~  x1+x2+x4, data = bsample.p, family = "binomial")
    L = 1 - predict(glmL, newdata = bsample.np, type = "response")
    
    W.NP = 1/bsample.np$includP
    bf11 = (W.NP-1)/(O*L) 
    p11 = sum(y*bf11)/sum(bf11) 
    return(c(p11))
  }
  
  # choose one of the two following lines as a basis for bootstrapping
  f.norm <- N * f11 / sum(f11)
  
  v = 0
  
  for (d in 1:D) {

  ## finite population bootstrap based on NP
  f.norm.rounded <- trunc(f.norm) + UPrandomsystematic(f.norm - trunc(f.norm)) # round weights using the method of Fellegi (1975)
  inflate.NP <- rep(1:nrow(sample.NP), times = f.norm.rounded)
  fakepop.NP <- sample.NP[inflate.NP, ]
  
  fakepop.NP$W.NP <- W.NP[inflate.NP]
  fakepop.NP$O <- O[inflate.NP]
  fakepop.NP$L <- L[inflate.NP]
  fakepop.NP$includP <- inclusionprobabilities(1 / fakepop.NP$W.NP, n = nrow(sample.P))
  
  fakepop.NP$includNP <- inclusionprobabilities((1 / f.norm)[inflate.NP], n = nrow(sample.NP))
  fakepop.NP$includNP.condS1 <- 1 - fakepop.NP$L
  fakepop.NP$includNP.condS0 <- (fakepop.NP$includNP - (1 - fakepop.NP$L)/fakepop.NP$W.NP) / (1 - 1/fakepop.NP$W.NP)
  eps = 1e-6
  fakepop.NP$includNP.condS0 <- ifelse(fakepop.NP$includNP.condS0 > eps, fakepop.NP$includNP.condS0, eps)
  fakepop.NP$includNP.condS1 <- ifelse(fakepop.NP$includNP.condS1 > eps, fakepop.NP$includNP.condS1, eps)
  
  
  bsample.p <- replicate(boot, fakepop.NP[as.logical(UPrandomsystematic(fakepop.NP$includP, eps = eps/10)), ], simplify = F)
  
  bsample.np <- lapply(1:boot, function(b) {
    S1 <- fakepop.NP[row.names(fakepop.NP) %in% row.names(bsample.p[[b]]), ]
    S0 <- fakepop.NP[!(row.names(fakepop.NP) %in% row.names(bsample.p[[b]])), ]
    rbind(S1[as.logical(UPrandomsystematic(S1$includNP.condS1, eps = eps/10)), ],
          S0[as.logical(UPrandomsystematic(S0$includNP.condS0, eps = eps/10)), ])
  })
  
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

capture.output(round(table, 4), file = "//cbsp.nl/Productie/secundair/MPOnderzoek/Werk/Combineren/Projecten/nonprobability/variance/var_dep.txt")

total.time = Sys.time() - start_time

save.image("//cbsp.nl/Productie/secundair/MPOnderzoek/Werk/Combineren/Projecten/nonprobability/variance/var_dep.RData")


library("xtable")
table = as.data.frame(table)
table$RB = 100*(table$est-table$sim.var)/table$sim.var
table$CI = 100*table$CI
xtable(table[,c(6,5)], digits = c(0,2,2), caption = 'Variance Estimates of Dep')


