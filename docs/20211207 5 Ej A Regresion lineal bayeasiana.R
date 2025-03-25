# EJEMPLO A DE LA SESIÓN 13 ####
# Simulación de los datos
# Tamaño de la muestra
n <- 10
beta_0 <- 100
beta_1 <- 20
sigma2 <- 25
# Variable predictora: x
set.seed(1984); x <- rnorm(n = n)
# Variable respuesta: y
set.seed(1313); y <- rnorm(n = n,
                           mean = beta_0 + beta_1*x,
                           sd = sqrt(sigma2))
# Crear una tabla de datos y graficar
print(datos <- data.frame(y, x))
plot(y ~ x, data = datos)
# Análisis frecuentista
# Ajuste del modelo lineal
ajuste <- lm(y ~ x, data = datos)
# Tabla de ANOVA
anova(ajuste)
# Resumen del ajuste
summary(ajuste)
#
# ESPECIFICACIÓN DEL MODELO EN JAGS ####
modelo_texto <- "model{

  # Modelo probabilístico
  for (i in 1:length(y)) {
    y[i] ~ dnorm(mu[i], tau)
    mu[i] = b0 + b1*x[i]
  }

  # Distribuciones a priori 'poco informativas'
    b0 ~ dnorm(0.0, 1.0/1.0e6)
    b1 ~ dnorm(0.0, 1.0/1.0e4)
    tau ~ dgamma(1.0/1.0e3, 1.0/1.0e3)

  # Funciones parametrales de interés
    sigma2 <- 1/tau
}"
#
# EJECUCIÓN DEL MÉTODO MCMC ####
# Preparativos
datos_jags <- as.list(datos)
parametros <- c("b0", "b1", "sigma2")
modelo_jags <- jags.model(
  file = textConnection(modelo_texto),
  data = datos_jags,
  n.chains = 3)
# Burn-in
update(object = modelo_jags, n.iter = 1e3)
# Generar las cadenas
modelo_coda <- coda.samples(
  model = modelo_jags,
  variable.names = parametros,
  n.iter = 5e3)
modelo_mcmc <- as.mcmc(do.call(rbind, modelo_coda))
#
# DIAGNÓSTICOS DE CONVERGENCIA ####
# a) Haga gráficos de trazas para evaluar convergencia y mezcla
# b) Haga gráficos de autocorrelación para evaluar si la mezcla es lenta
# c) Haga gráficos del estadístico de Gelman-Rubin para evaluar convergencia
# d) Haga un gráfico de las correlaciones cruzadas
# e) Haga gráficos de dispersión de los parámetros
#
# ANÁLISIS DE LA DISTRIBUCIÓN A POSTERIORI ####
# a) Haga gráficos de densidad de las distribuciones a posteriori
# b) Calcule el tamaño efectivo de las muestras de cada parámetro
# c) Calcule los estadísticos de Gelman-Rubin para cada parámetro
# d) Calcule intervalos de credibilidad para cada parámetro