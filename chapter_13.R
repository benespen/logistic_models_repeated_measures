# Examples from Chapter 13 of Statistical Rethinking

library(rethinking)
library(ggplot2)

## R code 13.1
data(reedfrogs)
d <- reedfrogs
str(d)

# view in RStudio
View(d)

# use Stan to estimate a single level model of the survival rates

## R code 13.2
# make the tank cluster variable
d$tank <- 1:nrow(d)

# map a list of columns from the data frame to pass to Stan
dat <- list(
  S = d$surv,
  N = d$density,
  tank = d$tank )

# approximate posterior
m13.1 <- ulam(
  alist(
    S ~ dbinom( N , p ) ,
    logit(p) <- a[tank] ,
    a[tank] ~ dnorm( 0 , 1.5 )
  ), data=dat , chains=4 , log_lik=TRUE )

# inspect the posterior
precis(m13.1 , depth=2)

# extend the model to be multi-level
## R code 13.3
m13.2 <- ulam(
  alist(
    S ~ dbinom( N , p ) ,
    logit(p) <- a[tank] ,
    a[tank] ~ dnorm( a_bar , sigma ) ,
    a_bar ~ dnorm( 0 , 1.5 ) ,
    sigma ~ dexp( 1 )
  ), data=dat , chains=4 , log_lik=TRUE )

# inspect the posterior
precis(m13.2 , depth=2)

## R code 13.4
compare( m13.1 , m13.2 )

## R code 13.5
# extract Stan samples
post <- extract.samples(m13.2)

# compute mean intercept for each tank
# also transform to probability with logistic
d$propsurv.est <- logistic( apply( post$a , 2 , mean ) )

ggplot( d , aes(x=tank , y=propsurv)) +
  geom_point() +
  labs(
    title = "Reed frog tadpole mortality",
    subtitle = "Proportion surviving",
    x = "Tank", y = "Proportion"
  )

# display raw proportions surviving in each tank
plot( d$propsurv , ylim=c(0,1) , pch=16 , xaxt="n" ,
      xlab="tank" , ylab="proportion survival" , col=rangi2 )
axis( 1 , at=c(1,16,32,48) , labels=c(1,16,32,48) )

# overlay posterior means
points( d$propsurv.est )

# mark posterior mean probability across tanks
abline( h=mean(inv_logit(post$a_bar)) , lty=2 )

# mark empirical mean probability across tanks
abline( h=mean(d$propsurv) , lty=3 )

# draw vertical dividers between tank densities
abline( v=16.5 , lwd=0.5 )
abline( v=32.5 , lwd=0.5 )
text( 8 , 0 , "small tanks" )
text( 16+8 , 0 , "medium tanks" )
text( 32+8 , 0 , "large tanks" )

## R code 13.6
# show first 100 populations in the posterior
plot( NULL , xlim=c(-3,4) , ylim=c(0,0.35) ,
      xlab="log-odds survive" , ylab="Density" )
for ( i in 1:100 )
  curve( dnorm(x,post$a_bar[i],post$sigma[i]) , add=TRUE ,
         col=col.alpha("black",0.2) )

# sample 8000 imaginary tanks from the posterior distribution
sim_tanks <- rnorm( 8000 , post$a_bar , post$sigma )

# transform to probability and visualize
dens( inv_logit(sim_tanks) , lwd=2 , adj=0.1 )
