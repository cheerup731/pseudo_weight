
### please run the 'data.R' file for importing population ### 

# probability sample on x1 
# nonprobability sample is selective at random (on x2 and x4)
# P(S*=0|S=1, x) = P(S*=0|x)


library(sampling)

set.seed(14)

colnames(pop) = c('x1', 'x2', 'x3', 'x4', 'y') 
# x1 = register year, x2 = engine type, x3 = fuel type, x4 = owner age, y = distance

N = nrow(pop)
# pop = pop[order(pop$x2),] # for stratified sampling

scenario = expand.grid(P = c(0.01, 0.1),
                       NP= c(0.05, 0.3, 0.5))

result = matrix(NA, nrow(scenario), 20)

M = 1000
i = m = 1

for (i in 1:nrow(scenario)){
  cat(nrow(scenario)-i)
  icf.P = scenario[i, 1] ; icf.NP = scenario[i, 2] # inclusion fraction 
  bias = mse = 0
  
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
    
    # estimate L
    sample.P$Sstar = Sstar[S == 1]
    glmL = glm(Sstar ~ x1+x2+x4, data = sample.P, family = "binomial")
    L = 1 - predict(glmL, newdata = sample.NP, type = "response")
    
    # formula (11) dep
    fdep = (W.NP-1)/(O*L) 
    
    # formula (18) ind
    find = 1+(W.NP-1)/O
    
    # direct estimate in probability sample (Kim)
    glmKim = glm(Sstar ~ x1+x2+x4, data = sample.P, family = "binomial", weights = W.P)
    fKim = 1/predict(glmKim, newdata = sample.NP, type = "response")
    
    # CLW
    xsa = cbind(rep(1,nrow(sample.NP)), as.matrix(sample.NP[,c('x1','x2','x4')])) # working matrix
    
    xsb = cbind(rep(1,nrow(sample.P)), as.matrix(sample.P[,c('x1','x2','x4')]))
    
    col_xsa = colSums(xsa)
    theta_0 = rep(0, length(col_xsa))
    theta_1 = solve(t(xsb)%*%(W.P*xsb), col_xsa-t(xsb)%*%W.P)
    
    while (abs(max(theta_1-theta_0))>10^(-8)){
      theta_0 = theta_1
      ps = 1/(1+exp(-xsb%*%theta_0))
      dev = tryCatch(solve(t(xsb)%*%(c(W.P*ps*(1-ps))*xsb), col_xsa - t(xsb)%*%c(W.P*ps)), error = function(e) dev <- 0)
      theta_1 = theta_0 + dev
    }
    
    fCLW = 1+exp(-xsa%*%theta_1)
    
    
    # estimates with A+B set (new)
    
    AB = rbind(sample.NP[,-5], sample.P[,-5])
    AB$Z = c(rep(1,nrow(sample.NP)), rep(0,nrow(sample.P)))
    
    # estimate O
    glmO = glm(Z ~ x1+x2+x4, data = AB, family = "binomial")
    O = exp(predict(glmO, newdata = sample.NP))
    
    # formula (18)
    ABfind = 1+(W.NP-1)/O
    
    # formula (20) EV
    ABfEV = W.NP/O
    
    # Valliant and Dever
    AB$W = 1
    AB$W[AB$Z == 0] = W.P*( round(sum(W.P))-nrow(sample.P) )/round(sum(W.P))
    glmV = glm(Z ~  x1+x2+x4, data = AB, family = "binomial", weights = W)
    fV = 1/predict(glmV, newdata = sample.NP, type = "response")
    
    # Wang 
    AB$W[AB$Z == 0] = length(W.P)*W.P/sum(W.P) # normalized weights
    glmW = glm(Z ~  x1+x2+x4, data = AB, family = "binomial", weights = W)
    oddW = exp(predict(glmW, newdata = sample.NP))
    fW = 1/oddW
    
    ### results
    est = c(mean(sample.NP$y), 
            sum(sample.NP$y*find)/sum(find), 
            sum(sample.NP$y*fdep)/sum(fdep), 
            sum(sample.NP$y*ABfind)/sum(ABfind), 
            sum(sample.NP$y*ABfEV)/sum(ABfEV),
            sum(sample.NP$y*fCLW)/sum(fCLW),
            sum(sample.NP$y*fV)/sum(fV),
            sum(sample.NP$y*fW)/sum(fW),
            sum(sample.NP$y*fKim)/sum(fKim))
    
    bias = (est - mean(pop$y))/mean(pop$y)/M + bias
    
    mse = (est - mean(pop$y))^2/M + mse
    
  }
  
  result[i,] = c(icf.NP, icf.P, bias*100, sqrt(mse))
  cat(result[i,])
  
}

colnames(result) = c('n.NP', 'n.P', 
                     'Naive', 'w_ind', 'w_dep', 'AB_w_ind', 'AB_EV', 'CLW', 'VD', 'Wang','KW',
                     'Naive', 'w_ind', 'w_dep', 'AB_w_ind', 'AB_EV', 'CLW', 'VD', 'Wang','KW')
View(result)

capture.output(round(result, 4), file = '//cbsp.nl/Productie/secundair/MPOnderzoek/Werk/Combineren/Projecten/nonprobability/est/ind.txt')

save.image("//cbsp.nl/Productie/secundair/MPOnderzoek/Werk/Combineren/Projecten/nonprobability/est/ind.RData")

library("xtable")
xtable(result[,c(1:11)], digits = c(0,rep(2,11)), caption = 'Relative Bias')
xtable(result[,c(1,2,12:20)], digits = c(0,2,2,rep(0,9)), caption = 'RMSE')
