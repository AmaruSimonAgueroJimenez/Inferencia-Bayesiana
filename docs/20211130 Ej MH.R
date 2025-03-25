## Ejemplo de Metropolis-Hastings

#yi ~ n(mu,1), i=1,...,n
#mu ~ t(0,1,1)

#logaritmo de la verosimilitud
lg = function(mu, n, ybar){
     mu2 = mu^2
     n*(ybar*mu - mu2/2) - log(1 + mu2)
}

mh = function(n, ybar, n_iter, mu_init, cand_sd){
     mu_out = numeric(n_iter)
     accpt = 0
     mu_now = mu_init
     lg_now = lg(mu=mu_now, n=n, ybar=ybar)

     for(i in 1:n_iter){
         mu_cand = rnorm(1, mean=mu_now, sd=cand_sd)
         lg_cand = lg(mu=mu_cand, n=n, ybar=ybar)
         lalpha = lg_cand - lg_now
         alpha = exp(lalpha)

         u = runif(1)
         if (u < alpha){
         mu_now = mu_cand
         accpt = accpt + 1
         lg_now = lg_cand
         }

     mu_out[i] = mu_now
     }
     
     list(mu=mu_out, accpt = accpt/n_iter)
}

##Datos

y = c(1.2, 1.4, -.5, .3, .9, 2.3, 1.0, .1, 1.3, 1.9)

ybar = mean(y)
n = length(y)

hist(y, freq=F, xlim=c(-1.0, 3.0))
points(y, rep(0.0, n))
points(ybar, 0.0, pch=19)
curve(dt(x, df=1), lty=2, add=T, col="red")

##Muestreo a posteriori
set.seed(43)
post = mh(n=n, ybar=ybar, n_iter=1e3, mu_init=0.0, cand_sd=3.0)
str(post)

library(coda)
traceplot(as.mcmc(post$mu))
#si acepta muy poco el modelo no es tan ?til porque se 
#queda muy pegado
#?c?mo puedo hacer que acepte siempre, 
#aumentar la tasa de aceptaci?n?
#disminuyendo la sd

post = mh(n=n, ybar=ybar, n_iter=1e3, mu_init=0.0, cand_sd=.05)
str(post)
traceplot(as.mcmc(post$mu))
#la tasa de aceptaci?n es muy alta, el modelo no aprende nada

#una cadena bonita se ve as?:
traceplot(as.mcmc(rnorm(1000,0,1)))
#recorre bien el espacio de par?metros

post = mh(n=n, ybar=ybar, n_iter=1e3, mu_init=0.0, cand_sd=.9)
str(post)
traceplot(as.mcmc(post$mu))
#esta cadena se ve mucho mejor

#ya sabemos que el sd=.9 funciona bien, as? que jugaremos con el valor inicial de mu
post = mh(n=n, ybar=ybar, n_iter=1e3, mu_init=30.0, cand_sd=.9)
str(post)
traceplot(as.mcmc(post$mu))
#la cadena estuvo s?per perdida y convergi? como en el 100
#entonces tengo que bajar el valor m?s cerca de donde convergi? para que
#sea m?s r?pido


##An?lisis a posteriori

post$mu_keep = post$mu[-c(1:100)]
str(post)
traceplot(as.mcmc(post$mu_keep))
#la misma cadena que ten?a antes, borr?ndole los primeros 100, queda mucho mejor

#cuando un proceso recorre bien un espacio de par?metros
#se dice que tiene un buen mixing

plot(density(post$mu_keep), xlim=c(-1,3)) #a posteriori
points(ybar, 0.0, pch=19)
curve(dt(x, df=1), lty=2, add=T, col="red") #a priori

autocorr.plot(as.mcmc(post$mu))
#ac? veo la memoria del proceso
#una muestra dependiente tiene el problema de qu? tan ?til es la aproximaci?n ...
#quiz?s no tengo realmente 1000 valores distintos sino menos

autocorr.diag(as.mcmc(post$mu))
#esto quiere decir que son muy dependientes las muestras que tengo
autocorr.plot(as.mcmc(post$mu), lag.max=500)

effectiveSize(as.mcmc(post$mu))
#tengo un tama?o muestral efectivo de ~12
#en el gr?fico vemos que la dependecia de 1 individuo afecta m?s o menos
#hasta el individuo 100, entonces 1/100 es independiente

thin_interval = 100
thin_indx = seq(from=thin_interval, 1e3, by=5)
length(thin_indx)

par(mfrow = c(2,1))
traceplot(as.mcmc(post$mu_keep))
traceplot(as.mcmc(post$mu[thin_indx]))
#se aprecia que en la segunda opci?n las muestras son m?s independientes 
#que en las muestras iniciales

autocorr.plot(as.mcmc(post$mu[thin_indx]))
#este es un comportamiento m?s adecuado
#esto se hace para quedarme con un tama?o de muestras m?s manejable
#pero no implica una diferencia en la inferencia
#se usa cuando tienes muuuchos datos y as? hay la misma info con menos datos
effectiveSize(as.mcmc(post$mu[thin_indx]))

#otra opci?n es, al simular, saltarse el guardado de datos cada cierto rato

############################################################

##Procedimiento de

post1 = mh(n=n, ybar=ybar, n_iter=1e3, mu_init=30.0, cand_sd=.9)
post2 = mh(n=n, ybar=ybar, n_iter=1e3, mu_init=10.0, cand_sd=1.2)
post3 = mh(n=n, ybar=ybar, n_iter=1e3, mu_init=0.0, cand_sd=.9)
post4 = mh(n=n, ybar=ybar, n_iter=1e3, mu_init=-10.0, cand_sd=1.2)
post5 = mh(n=n, ybar=ybar, n_iter=1e3, mu_init=-30.0, cand_sd=.9)

pmc = mcmc.list(as.mcmc(post1$mu),as.mcmc(post2$mu),as.mcmc(post3$mu),as.mcmc(post4$mu),as.mcmc(post5$mu))

par(mfrow = c(1,1))
traceplot(pmc)
#n la pr?ctica uno debiera correr varias cadenas con distintos par?metros

gelman.diag(pmc)
#indica si no hay evidencia de convergencia entre las cadenas
#cuando obtengo valores muy cercanos a 1 estamos bien
#estamos convergiendo en las cadenas
#lo m?s peligroso ser?a que no convergiera
#otra estrategia es generar una cadena muy muy larga pero con eso se consumen m?s recursos

#entonces, es importante ver si la distribuci?n tiene buen mixing
#y ver cu?nta dependencia hay para ver si se puede obtener una muestra m?s eficiente

#estimador de bayes
mean(post3$mu)
ybar #verosimilitud

summary(as.mcmc(post2$mu))
#me da todo lo que necesito para el modelo

