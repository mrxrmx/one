call_option <- function(x, K) {
  answer <- max(x-K, 0)
  return(answer)
}

put_option <- function(x, K) {
  answer <- max(K-x, 0)
  return(answer)
}

getExternalCall <- function (case, vstrike) {
  if (case == 1) {
    if (vstrike == 70) {
      result <- c (30.0000, 31.1991, 31.8758, 32.4202, 32.9421, 33.4518, 33.9508, 34.4396, 34.9188, 
                   35.3887, 35.8498, 36.3024, 36.7469, 37.1836, 37.6128, 38.0347)
    } else
      if (vstrike == 100) {
        result <- c (0, 4.4034, 5.6275, 6.7343, 7.7729, 8.7569, 9.6945, 10.5918, 11.4535, 
                     12.2835, 13.0847, 13.8598, 14.6111, 15.3404, 16.0493, 16.7394)
      } else
        if (vstrike == 140) {
          result <- c (0.0000, 0.0022, 0.0089, 0.0166, 0.0275, 0.0434, 0.0663, 0.0990, 0.1450, 
                       0.2088, 0.2958, 0.4120, 0.5638, 0.7564, 0.9937, 1.2769)
        }
  } else
    if (case == 2) {
      if (vstrike == 70) {
        result <- c (30.0000, 31.0144, 31.7409, 32.2622, 32.7232, 33.1607, 33.5864, 34.0046, 
                     34.4169, 34.8241, 35.2265, 35.6242, 36.0172, 36.4058, 36.7899, 37.1697)
      } else
        if (vstrike == 100) {
          result <- c (0, 5.0998, 6.4175, 7.4775, 8.4337, 9.3287, 10.1796, 10.9950, 11.7801,
                       12.5385, 13.2728, 13.9854, 14.6778, 15.3518, 16.0085, 16.6492)
        } else
          if (vstrike == 140) {
            result <- c (0.0000, 0.3192, 0.7359, 1.0613, 1.3588, 1.6533, 1.9542, 2.2650, 2.5869, 
                         2.9202, 3.2648, 3.6201, 3.9858, 4.3611, 4.7455, 5.1382)
          }
      
    } else
      
      if (case == 3) {
        if (vstrike == 70) {
          result <- c (30.0000, 32.0027, 33.9058, 35.6236, 37.2407, 38.7720, 40.2244, 41.6041, 42.9171, 
                       44.1691, 45.3653, 46.5102, 47.6081, 48.6626, 49.6770, 50.6542)
        } else
          if (vstrike == 100) {
            result <- c (0, 9.7738, 13.5681, 16.6689, 19.3726, 21.7953, 24.0021, 26.0355, 27.9255, 29.6939, 
                         31.3579, 32.9307, 34.4229, 35.8433, 37.1992, 38.4965)
          } else
            if (vstrike == 140) {
              result <- c (0, 1.5514, 3.5857, 5.6995, 7.8513, 9.9831, 12.0600, 14.0654, 15.9935, 
                           17.8438, 19.6187, 21.3220, 22.9577, 24.5301, 26.0434, 27.5012)
            }
        
      }
  
}

#Euler_S1<-function(S,r,Vt,h,w1){   ## <---f--------kai beta ir lambda = 1
#  next_S=S*exp((r-abs(Vt/2))*h+sqrt(abs(Vt))*w1)
#  
#  return(next_S)
#} 

Euler_S2 <- function (S, r, Vt, h, w3, lambda, beta) {
  
  next_S = S * exp((r-Vt*(S^(2*beta-2))*0.5*lambda^2)*h + lambda*sqrt(Vt)*S^(beta-1)*w3)
  
  return(next_S)
} ## bendras atvejis

Euler_V<-function(V0, teta, k, h, sigma, alfa, w1, skaicius){
  if (skaicius==1) {next_v = Euler_V_absorption (V0, teta, k, h, sigma, alfa, w1) }
  if (skaicius==2) {next_v = Euler_V_reflection (V0, teta, k, h, sigma, alfa, w1) }
  if (skaicius==3) {next_v = Euler_V_higham (V0, teta, k, h, sigma, alfa, w1) }
  if (skaicius==4) {next_v = Euler_V_partial (V0, teta, k, h, sigma, alfa, w1) }
  if (skaicius==5) {next_v = Euler_V_full (V0, teta, k, h, sigma, alfa, w1) }
  return(next_v)
}

Euler_V_absorption <- function(V0, teta, k, h, sigma, alfa, w1){
  
  next_V = max(V0,0) + k*(teta-max(V0, 0))*h + sigma*(max(V0,0))^(alfa)*w1
  
  return(next_V)
}

Euler_V_reflection <- function (V0,teta,k,h,sigma,alfa,w1) {

  next_V = abs(V0) + k*(teta-abs(V0))*h + sigma*(abs(V0))^(alfa)*w1
  
  return(next_V)
}


Euler_V_higham <- function (V0,teta,k,h,sigma,alfa,w1) {
  
  next_V = V0 + k*(teta-V0)*h + sigma*(abs(V0))^(alfa)*w1
  
  return(next_V)
}

Euler_V_partial <- function (V0,teta,k,h,sigma,alfa,w1) {
  
  next_V = V0 + k*(teta-V0)*h + sigma*(max(V0,0))^(alfa)*w1
  
  return(next_V)
}


Euler_V_full <- function (V0,teta,k,h,sigma,alfa,w1) {
  
  next_V = V0 + k*(teta-max(V0,0))*h + sigma*(max(V0,0))^(alfa)*w1
  
  return(next_V)
}

Naujas_Euler <- function (Pr,T,h,S0,V0,r,k,teta,sigma,strike,alfa,vvid,rho,beta,lambda,skaicius) {
  n = (T-Pr)/h + 1
  E <- rep (0, n)
  V = vector()
  
  for (j in 1:vvid) {
    Vt=V0
    St = S0
    E[1] = E[1] + call_option(St,strike)

    # w1 ir W3 priklausomi 
    w1 = rnorm(n,0,sqrt(h))
    w2 = rnorm(n,0,sqrt(h))
    w3 = rho*w1 + sqrt(1-rho*rho) * w2
    
    for (i in 1:(n-1)){
      V[i] = Euler_V (Vt, teta, k, h, sigma, alfa, w1[i], skaicius)
      
      if (skaicius==1) { Vt_ = max (Vt, 0) }
      if (skaicius==2) { Vt_ = abs (Vt) }
      if (skaicius==3) { Vt_ = abs (Vt) }
      if (skaicius==4) { Vt_ = max (Vt, 0) }
      if (skaicius==5) { Vt_ = max (Vt, 0) }
      
      St = Euler_S2 (St, r, Vt_, h, w3[i], lambda, beta)
      
      Vt = V[i]
      #St = exp (Xt)

# Errr modulius reikia taikyt V_t      
# S_t logaritmuotas, tad nebus neigiamas (imsim eksponentÄ™)...      
      
      #Next_Vt=Euler_V(Vt,teta,k,h,sigma,alfa,w3,skaicius)

#      if (skaicius==1 | skaicius==4 | skaicius==5){
#        St=Euler_S2(St,r,max(V[i],0),h,w1[i],lambda,beta)
#      }
#      
#      if (skaicius==2 | skaicius==3){
#        St=Euler_S2(St,r,abs(V[i]),h,w1[i],lambda,beta)
#      }
      
      
      #Vt=V[i]
      
      E[i+1]=E[i+1]+call_option(St,strike)
    }
  }
  
  E = E / vvid
  return(E)
}

#skaicius = 1 ---> Euler_V_absorption
#skaicius = 2 ---> Euler_V_reflection
#skaicius = 3 ---> Euler_V_higham
#skaicius = 4 ---> Euler_V_partial
#skaicius = 5 ---> Euler_V_full
#Pr,T,h,S0,V0,r,k,teta,sigma,strike,alfa,vvid
#plot(Naujas_Euler(0,3,0.1,100,0.09,0.02,1,0.12,0.2,100,0.5,1000,-0.5),type="l")


# HESTON ----------------------------------
Case = 3

if (Case == 1) {
  # Case 1
  vPr = 0
  vs0 = 100
  vteta = 0.04
  vy0 = vteta
  vrho = -0.9
  vk = 0.5
  vsigma = 1
  vr = 0
} else if (Case==2) {
  # Case 2
  vPr = 0 
  vs0 = 100
  vteta = 0.04
  vy0 = vteta
  vrho = -0.5
  vk = 0.3
  vsigma = 0.9
  vr = 0
} else if (Case == 3) {
  # Case 3
  vPr = 0
  vs0 = 100
  vteta = 0.09
  vy0 = vteta
  vrho = -0.3
  vk = 1
  vsigma = 1
  vr = 0
} 

# -----------------------------------------------------------
  
  
S0 = vs0
V0 = vy0
# r = miu
r = vr 
k = vk
teta = vteta
sigma = vsigma


# Äia klaida - rho nelygus 0 !!!
rho = vrho

lambda = 1
alfa = 0.5
beta = 1

Pr = 0
T = 10
h = 0.01 #0.05
strike = 70
vvid = 100000 #50000

E1 <- Naujas_Euler(Pr,T,h,S0,V0,r,k,teta,sigma,strike,alfa,vvid,rho,beta,lambda,1)
print ("E1 baigtas")
E2 <- Naujas_Euler(Pr,T,h,S0,V0,r,k,teta,sigma,strike,alfa,vvid,rho,beta,lambda,2)
print ("E2 baigtas")
E3 <- Naujas_Euler(Pr,T,h,S0,V0,r,k,teta,sigma,strike,alfa,vvid,rho,beta,lambda,3)
print ("E3 baigtas")
E4 <- Naujas_Euler(Pr,T,h,S0,V0,r,k,teta,sigma,strike,alfa,vvid,rho,beta,lambda,4)
print ("E4 baigtas")
E5 <- Naujas_Euler(Pr,T,h,S0,V0,r,k,teta,sigma,strike,alfa,vvid,rho,beta,lambda,5)
print ("E5 baigtas")

result <- getExternalCall (Case, strike)
result <- head (result, T+1)

max = max (max(E1), max(E2), max(E3), max(E4), max(E5), max(result))
min = min (min(E1), min(E2), min(E3), min(E4), min(E5), min(result))
plot (c(0,(T/h)), c(min, max), type="n", xlab="Time", ylab="Call", xlim = range(0:T), main = "Eulers")
x <- seq(Pr, T, h)
lines(x, E1, col="black", lwd=1)
lines(x, E2, col="red", lwd=1) 
lines(x, E3, col="green", lwd=1)
lines(x, E4, col="blue", lwd=1) 
lines(x, E5, col="orange", lwd=1) 
y<-seq(Pr, T, 1)
points(y, result, lwd=1) 
print ("Baigta")
print ("-----------------")

#lines(Naujas_Euler(0,10,0.05,100,0.04,0,0.5,0.04,1,100,0.5,3000,-0.9,1,1),col="red",lwd=2.5) 
#legend('bottomright',c("pirmas","antras"), lty=c(1,1), lwd=c(2.5,2.5),col=c("green","blue"), box.lty=0)
