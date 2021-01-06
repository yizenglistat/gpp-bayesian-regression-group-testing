
ise_all = c()

for(rep in 1:200){
	source('Equal_Homo.R')
	ise_all = c(ise_all, ise)
	med = median(ise_all)
	iqr = IQR(ise_all)
	if(rep%%10==0){
		print(rep)
		print(cbind(med,iqr))
	}
}


rm(list=ls())
n = 100
x = rnorm(n,5,1)
y = rpois(n,5)
d = x-y

xbar = mean(x)
ybar = mean(y)
dbar = mean(d)

sx = sd(x)
sy = sd(y)
sp = ((n-1)*sx+(n-1)*sy)/(2*n-2)

alpha = 0.05
c = seq(-1,1,0.01)


# one-sided for mu_x - mu_y = 0 vs mu_x - mu_y > c

T_two = qt(1-alpha,n+n-2) - (c-0)/(sp*sqrt(1/n+1/n))
T_pair = qt(1-alpha,n-1) - (c-0)/(sd(d)*sqrt(1/n))

pi_two = 1- pt(T_two,n+n-2)
pi_pair = 1- pt(T_pair,n-1)

plot(c, pi_two,type='l',col='red',ylim=c(0,1),lty=1,ylab='Power',xlab='Alternative')
lines(c, pi_pair,type='l',col='blue',lty=2)
legend('left',c('two-sample','paired'),lty=1:2,col=c('red','blue'))



# one-sided for mu_x - mu_y = 0 vs mu_x - mu_y < c

T_two = qt(alpha,n+n-2) - (c-0)/(sp*sqrt(1/n+1/n))
T_pair = qt(alpha,n-1) - (c-0)/(sd(d)*sqrt(1/n))

pi_two = pt(T_two,n+n-2)
pi_pair = pt(T_pair,n-1)

plot(c, pi_two,type='l',col='red',ylim=c(0,1))
lines(c, pi_pair,type='l',col='blue')

N = 1000
out = c()
for(i in 1:N){
	x = c()
	while(T){
		x = c(x, runif(1))
		condition1 = length(x) > 1
		if(condition1){
			condition2 = tail(x,2)[2] - tail(x,2)[1] > 0
			if(condition2){
				out = c(out, length(x))
				break
			}
		}
		
	}
}
mean(out)
