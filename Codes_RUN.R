rm(list = ls())  
gc()
dloc = 'c:/'
setwd(dloc)

if(!require(R2OpenBUGS)) install.packages("R2OpenBUGS")
if(!require(evd)) install.packages("evd")
library(splines)
library(evd)

n = 60   
nsim = 100
x = 1:n

cloc = rep( c(0, 40), each = 2 )
slp = rep( c( 0, 0.00732 * 120 / n ), 2) 

for ( dd in c(1, 2, 3, 4) )
{  
set.seed(666)
dat.sim = matrix(0, nsim, n)
for(isim in 1:nsim)
{ 
    {  
        dat.sim[isim , 1:(n/2)] = rgev( n / 2, shape = 0.1, 
        loc = 150 , scale = exp( 3 + slp[dd] * ( (1:(n/2)) ) ) )
        dat.sim[isim , ( 1 + n/2 ):n] = rgev(n / 2, shape = 0.1, 
        loc = 150 + cloc[dd], scale = exp( 3 + slp[dd] * ( ((1 + n / 2):n) ) ) )
    }
}
}

## Codes for the first replicaiton ## 

dd = 1
y = dat.sim[dd, ]
cx = (1:length(y)) / 100  
dx1 = dx2 = cbind(1, cx, cx^2, cx^3)
fit = fgev(y, std.err = FALSE)
    
for ( m in c(1,2,4,5) ) 
{  
    if ( m == 1 ) model.file='changepoint1_th.txt'
    if ( m == 2 ) model.file='changepoint2_th.txt'
    if ( m == 4 ) model.file='changepoint0_th.txt'
    if ( m == 5 ) model.file='changepoint0v_th.txt'

    fit = fgev(y)
    n = length(y) 
    
    if ( m == 1 ) 
    { 
      dat <- list(y = y, T = n, ymed = median(y), m.mle = fit$estimate[1], 
                  x1 = dx2[,1], x2 = dx2[,2], x3 = dx2[,3] ) 
      
      inits <- list(
        list(  mu.1 = rnorm(1, fit$estimate[1], 0.0000), sigma = 30, xi = 0.0, tau1 = 0.01, 
               dof1 = 10 ), 
        list(  mu.1 = rnorm(1, fit$estimate[1], 0.0000), sigma = 30, xi = 0.0, tau1 = 0.01,
               dof1 = 10 )
      )
      
      parameters <- c('mu','sigma','nu','delta','eta','xi')
    }
    
    if ( m == 2 ) 
    { 
      dat <- list(y = y, T = n, ymed = median(y), m.mle = fit$estimate[1],
                  x1 = dx2[,1], x2 = dx2[,2], x3= dx2[,3] ) 
      
      inits <- list(
        list(  mu.1 = rnorm(1, fit$estimate[1], 0.0000), a = 30, b = 0, 
               tau1 = 0.01, tau2 = 0.01, tau3 = 0.01, tau4 = 0.01, dof1 = 10,  
               eta = rnorm(1, 0.1, 0.001), xi = 0.0 ), 
        
        list(  mu.1 = rnorm(1, fit$estimate[1], 0.0000), a = 30, b = 0, 
               tau1 = 0.01, tau2 = 0.01, tau3 = 0.01, tau4 = 0.01, dof1 = 10,   
               eta = rnorm(1, 0.1, 0.001), xi = 0.0 )
      )
      parameters <- c('mu','sigma','nu','delta','eta','xi','a','b')
    }
    
    if ( m == 4 ) 
    { 
      dat <- list(y = y, T = n, x1 = dx1[,1], x2 = dx1[,2], 
                  ymed = median(y), m.mle = fit$estimate[1] ) 
      
      inits <- list(
        list( mu.1 = rnorm(1, fit$estimate[1], 0.00), a = 30, b = 0, p1 = 0,   
              c = 0 , p2 = 0, tau3 = 0.01, tau4 = 0.01, dof1 = 10, dof2 = 10,  
              xi = 0.0, tau1 = 0.01, tau2 = 0.01, sigma = 30 ),
        list( mu.1 = rnorm(1, fit$estimate[1], 0.00), a = 30, b = 0, p1 = 0,   
              c =0, p2 = 0, tau3 = 0.01, tau4 = 0.01, dof1 = 10, dof2 = 10,   
              xi = 0.0, tau1 = 0.01, tau2 = 0.01, sigma = 30 )
                    )
      parameters <- c('mu','sigma','xi')
    }
    
    if ( m == 5 )
    {
      dat <- list(y = y, T = n, x1 = dx1[,1], x2 = dx1[,2], 
                  ymed = median(y), lsigma0 = fit$estimate[2], m.mle = fit$estimate[1] )
      
      inits <- list(
        
        list( mu.1 = rnorm(1, fit$estimate[1], 0.00), a = 30, b = 0, p1 = 0,   
              c = 0 , p2 = 0, tau3 = 0.01, tau4 = 0.01, dof1 = 10,  
              xi = 0.0, tau1 = 0.01, tau2 = 0.01 ),
        list( mu.1 = rnorm(1, fit$estimate[1], 0.00), a = 30, b = 0, p1 = 0,   
              c =0, p2 = 0, tau3 = 0.01, tau4 = 0.01, dof1 = 10,  
              xi = 0.0, tau1 = 0.01, tau2 = 0.01 )
      )
      parameters <- c('mu','sigma','xi', 'a', 'b')
      
    }
    
    n.iter = 5000 * 1  
    if (m <= 3) n.iter = n.iter / 1.0  
    
    psfit <- bugs(data = dat,inits = inits, parameters.to.save = parameters,
                  model.file = model.file, debug = F, DIC = T, codaPkg = F,
                  n.chains = 2, n.iter = n.iter, n.thin = 1 * 30 )  
                                                             
    psfit$sims.matrix -> vv
    vv[, ncol(vv)] -> vv
    DVV = mean(vv) + var(vv)/2 
    
    if (m == 1) list(psfit, DVV) -> res_1
    if (m == 2) list(psfit, DVV) -> res_2
    if (m == 4) list(psfit, DVV) -> res_0
    if (m == 5) list(psfit, DVV) -> res_v
    
    print(m) 
}

save(file = 'Res1', res_1, res_2, res_0, res_v2, y) 



