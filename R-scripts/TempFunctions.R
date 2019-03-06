require(minpack.lm)

briere<-function(t, c, T0, Tm){
  b=c()
  for (i in 1:length(t)){
    if (t[i]>T0 && t[i]<Tm)  {b[i]<-(c*t[i]*(t[i]-T0)*sqrt(Tm-t[i]))} else {b[i]<-(0)}
  }
  b
}

quad<-function(t, c, T0, Tm){
  b=c()
  for (i in 1:length(t)){
    if (t[i]>T0 && t[i]<Tm)  {b[i]<-(-c*(t[i]-T0)*(t[i]-Tm))} else {b[i]<-(0)}
  }
  b
}

eppley<-function(t, a, b2){
  b=c()
  for (i in 1:length(t)){
    b[i]<-a*exp(b2*t[i])
  }
  b
}

eppleynorberg<-function(t, z, w, a, b2){
  b=c()
  for (i in 1:length(t)){
    b[i]<- (1 - ((t[i] - z)/w)^2)*a*exp(b2*t[i])
  }
  b
}

temps <- seq(1,40,0.5)

plot(briere(temps, 0.0001, 5, 33) ~ temps, type = "l")
plot(quad(temps, 0.001, 5, 33) ~ temps, type = "l")
plot(eppley(temps, 0.59, 0.0633) ~ temps, type = "l")

plot(norberg(temps, 25, 10, 0.59, 0.0633) ~ temps, type = "l", ylim = c(0,5), xlim = c(0,40))
lines(norberg(temps, 25, 15, 0.59, 0.0633) ~ temps, col = "red")
lines(norberg(temps, 20, 10, 0.59, 0.0633) ~ temps, col = "black", lty = 2)
lines(norberg(temps, 20, 15, 0.59, 0.0633) ~ temps, col = "red", lty = 2)