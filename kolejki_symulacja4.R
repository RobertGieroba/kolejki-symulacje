N <- 30000
#warunek początkowy
n1 <- 25
#alpha <- c(1,3/5)
alpha <- c(0.127,0.254)
arrival_times0 <- rbind(cumsum(c(rep(0,n1),rexp(N-n1, alpha[1]))),cumsum(c(rep(0,n1),rexp(N-n1,alpha[2]))))

#v <- rbind(runif(N),runif(N,min=2/3))
#v <- rbind(sample(c(0.5, 1, 1.5), N+n1, T),sample(c(0.2, 0.3, 0.4, 0.6, 1.5), N+n1, T))
v <- rbind(rgeom(N, 1/2)+1,rgeom(N, 1/2)+2)
w <- v
t_0 <- min(arrival_times0)
arrival_times <- arrival_times0 - t_0

min2 <- function(x){
  min(x[which(x>0)])
}
queued_jobs <- function(){
  which(arrival_times==0 & w>0,arr.ind = T)
}
serviced_job <- function(){
  srpt <- which(w==min2(w[queued_jobs()]),arr.ind = T)
  srpt <- matrix(srpt[arrival_times0[srpt]==min(arrival_times0[srpt])],ncol=2)
  srpt[1,]
}
unit_vector <- function(n,k){
  x <- rep(0,n)
  x[k] <- 1
  x
}
pos <- function(x){
  ifelse(x>0,x,0)
}

t <- 0
i <- 1
timestamps <- rep(0, N)

Q <- c(n1,n1)
if(Q[1]+Q[2]==0){
  Q <- unit_vector(2,serviced_job()[1])
}
W <- c(sum(w[1,1:n1]),sum(w[2,1:n1]))
if(W[1]+W[2]==0){
  W <- w[serviced_job()]*Q
}
Q_history <- matrix(0, ncol = 2, nrow = N)
W_history <- matrix(0, ncol = 2, nrow = N)
while(t<=10000 & i <= N){
  if(sum(Q)>=1){
    t2 <- ifelse(min2(w[queued_jobs()])<min2(arrival_times), min2(w[queued_jobs()]), min2(arrival_times))
    w[serviced_job()[1],serviced_job()[2]] <- pos(w[serviced_job()[1],serviced_job()[2]]-t2)
  } else {
    t2 <- min2(arrival_times)
  }
  t <- t+t2
  Q_history[i,] <- Q
  W_history[i,] <- W
  timestamps[i] <- t
  i <- i+1
  arrival_times <- pos(arrival_times-t2)
  Q <- c(-sum(queued_jobs()[,1]-2),sum(queued_jobs()[,1]-1))
  W <- c(-sum(w[queued_jobs()]*(queued_jobs()[,1]-2)),sum(w[queued_jobs()]*(queued_jobs()[,1]-1)))
}
# plot(timestamps, W_history[,2], type='l', col='red',xlab='t',ylab='W_k(t)')
# lines(timestamps, W_history[,1], type='l', col='blue')
# legend(35, legend=c("W_1(t)", "W_2(t)"), col=c("blue", "red"),lty=1)
# plot(timestamps, Q_history[,2], type='l', col='red',xlab='t',ylab='Q_k(t)')
# lines(timestamps, Q_history[,1], type='l', col='blue')
# legend(35, legend=c("Q_1(t)", "Q_2(t)"), col=c("blue", "red"),lty=1)
# hist(as.vector(w[queued_jobs()]), main = 'Residual transmission times of jobs present in the queue',xlab='Residual transmission time')
plot(timestamps[1:7500], W_history[1:7500,1]/(W_history[1:7500,1]+W_history[1:7500,2]), type='l', xlab='t',ylab='W_1(t)/W_Sigma(t)', col = 'red',
     cex.lab=1.5, cex.axis=2, ylim = c(0, 1))
abline(h=1/5, lty=4, lwd=3)
plot(timestamps[1:7500], W_history[1:7500,2], type='l', xlab='t',ylab='', col='red',
     cex.lab=1.5, cex.axis=2)
lines(timestamps[1:7500], 4/5*(W_history[1:7500,1]+W_history[1:7500,2]), type='l', col='blue')
legend(200, legend=c("W_1(t)", "p_1*Delta_1*W_Sigma(t)"), col=c("red", "blue"),lty=1,
       cex=1.75, lwd = 4)

mean_w <- c()
W_history2 <- matrix(0,nrow=10,ncol=2)
for(i in 1:20){
  N <- 6000
  alpha <- c(1,3/5)
  arrival_times0 <- rbind(cumsum(rexp(N, alpha[1])),cumsum(rexp(N,alpha[2])))
  v <- rbind(runif(N),runif(N,min=2/3))
  w <- v
  t_0 <- min(arrival_times0)
  arrival_times <- arrival_times0 - t_0
  t <- 0
  timestamps <- c(0)
  Q <- unit_vector(2,serviced_job()[1])
  W <- w[serviced_job()]*unit_vector(2,serviced_job()[1])
  while(t<=3000){
    if(sum(Q)>=1){
      t2 <- ifelse(min2(w[queued_jobs()])<min2(arrival_times), min2(w[queued_jobs()]), min2(arrival_times))
      w[serviced_job()] <- pos(w[serviced_job()]-t2)
    } else {
      t2 <- min2(arrival_times)
    }
    t <- t+t2
    timestamps <- c(timestamps, t)
    arrival_times <- pos(arrival_times-t2)
    Q <- c(-sum(queued_jobs()[,1]-2),sum(queued_jobs()[,1]-1))
    W <- c(-sum(w[queued_jobs()]*(queued_jobs()[,1]-2)),sum(w[queued_jobs()]*(queued_jobs()[,1]-1)))
  }
  mean_w <- c(mean_w, mean(w[queued_jobs()]))
  W_history2[i,] <- W
}
hist(mean_w, main = 'Means of residual transmission times of jobs present in the queue',xlab='Residual transmission time')
a <- (W_history2[,1]/apply(W_history2, 1, sum))
a <- a[which(is.nan(a)==F)]
summary(a)
hist(a, main='Histogram of W_1/W_Sigma (Mean: 0.355, expected value: 0.357)', xlab = 'W_1/W_Sigma')
