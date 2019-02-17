
### load libraries (particularly before reading in data)-----
library(highcharter)
library(xts)
library(dygraphs)
library(nanotime)
## to use milliseconds
opts <- options(digits.secs = 15)
## to display milliseconds ---
options("digits.secs" = 15)
### read in book data------

### maybe compare bid/offer vol between futures and options--we have an option book built----

bookNov9ES <- readRDS("../../../../data/built_orderbooks/book/ES/1109/book_sending_ESZ6_nov9.rds")

trade_size_dollars <- bookNov9ES$tradeSize * (bookNov9ES$price / 100) * 50

## get minute bid/ask data--------

trade_size_dollars_minute <- xts::to.minutes(trade_size_dollars)[, 4]
change_trade_price_minute <- diff(xts::to.minutes(bookNov9ES$price)[, 4])
## trade_size_minute <- xts::to.minutes(bookNov9ES$tradeSize)
# bid_minute <- xts::to.minutes(bookNov9ES$bidP1)
# offer_minute <- xts::to.minutes(bookNov9ES$offerP1)

### kyles lambda calc-----

### To estiamte Kyle's Lambda we use: lambda = |delt(price)| / volume in $

price <- change_trade_price_minute[-1]

volume <- trade_size_dollars_minute[-1]

## remove unnecessary data------

rm(bookNov9ES)
rm(trade_size_dollars)
rm(trade_size_dollars_minute)
rm(change_trade_price_minute)

## Should we kalman filter Kyle's lambda???

#{{{

lik <- function(theta, volume, price){


    ## R and Q transformed below (squared) -- need to sqrt in filter
    R <- theta[1]
    Q <- theta[2]
    F <- theta[3]
    mu <- theta[4]
    
    beta_tt <- rep(0, length(volume))
    beta_tt_1 <- rep(0, length(volume))

    eta <- rep(0, length(volume))
    f <- rep(0, length(volume))

    Ptt <- rep(0, length(volume))

    Ptt_1 <- rep(0, length(volume))
        
    beta_tt[1] <- lm(price~volume)$coef[2]
    Ptt[1] <- (summary(lm(price~volume))$coef[4])^2

    for(i in 2:length(volume)){
    ## Prediction
        beta_tt_1[i] <- mu + F*beta_tt[i-1]
        Ptt_1[i] <- F*Ptt[i-1]*F+Q^2
        eta[i] <- price[i]-volume[i]*beta_tt_1[i]
        f[i] <- volume[i]*Ptt_1[i]*volume[i]+R^2
    ## Updating
        beta_tt[i] <- beta_tt_1[i]+Ptt_1[i]*volume[i]*(1/f[i])*eta[i]
        Ptt[i] <- Ptt_1[i]-Ptt_1[i]*volume[i]*(1/f[i])*volume[i]*Ptt_1[i]
    }

    ## This is returning Inf because length(volume) > 700
    ## logl <- -0.5 * sum(log((((2 * pi)^{length(volume)}) * abs(f))[-1])) - 0.5 * sum(eta * eta * (1 / f),na.rm=T)

    ## this gives exact same result as above on subsample -- use
    logl <- -0.5 * sum(log((abs(f))[-1])) - 0.5 * sum(eta * eta * (1 / f),na.rm=T)

    return(-logl)
}

theta.start <- c(0.01, 0.01, 0.1, 0.1)
max.lik.optim <- optim(theta.start, lik, volume=volume, price=price, hessian = TRUE)

## Run though filter to get betas

R.hat <- (max.lik.optim$par[1])^2
Q.hat <- (max.lik.optim$par[2])^2
F.hat <- max.lik.optim$par[3]
mu.hat <- max.lik.optim$par[4]

beta_tt <- rep(0, length(volume))
    beta_tt_1 <- rep(0, length(volume))

    eta <- rep(0, length(volume))
    f <- rep(0, length(volume))

    Ptt <- rep(0, length(volume))

    Ptt_1 <- rep(0, length(volume))
        
    beta_tt[1] <- lm(price~volume)$coef[2]
    Ptt[1] <- (summary(lm(price~volume))$coef[4])^2

    for(i in 2:length(volume)){
    ## Prediction
        beta_tt_1[i] <- mu.hat + F.hat*beta_tt[i-1]
        Ptt_1[i] <- F.hat*Ptt[i-1]*F.hat+Q.hat
        eta[i] <- price[i]-volume[i]*beta_tt_1[i]
        f[i] <- volume[i]*Ptt_1[i]*volume[i]+R.hat
    ## Updating
        beta_tt[i] <- beta_tt_1[i]+Ptt_1[i]*volume[i]*(1/f[i])*eta[i]
        Ptt[i] <- Ptt_1[i]-Ptt_1[i]*volume[i]*(1/f[i])*volume[i]*Ptt_1[i]
    }
    logl <- -0.5*sum(log((((2*pi)^length(volume))*abs(f))[-1]))-.5*sum(eta*eta*(1/f),na.rm=T)

dygraph(as.xts(0.01 / beta_tt, order.by = index(twtr.s.d)[-1]))




par(mfrow=c(2,1))
plot(ts(0.01 / beta_tt, start =c(2000,2), frequency=12))
plot(ts(Ptt, start =c(2000,2), frequency=12))


par(mfrow=c(1,1))
pdf("lambda.pdf")
plot.xts(as.xts((0.01 / beta_tt), order.by = index(twtr.s.d)[-1]), ylab = "$", main = "Volume Required to Move Price by $0.01")
dev.off()
#}}}
