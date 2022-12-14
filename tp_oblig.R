# Genero una muestra aleatoria de Bernoulli de prob p y de tamaño n 
# El codigo mas formateado se encuentra en el siguiente colab
# https://colab.research.google.com/drive/1V3YPi8S0r9cAOL3UgJiguCJ5Kz94keli?usp=sharing


generar_muestra <- function(n,p) {
  return (rbinom(n,1,p))
  }


  # Estimacion Metodo 1
estimacion_metodo_1 <- function(muestra, p_real) {

  promedio <- mean(muestra)
  n <- length(muestra)
  sum_Xi <- sum(muestra)
  quantil <- qnorm(0.975) # z_{1-\alpha/2}

  k <- quantil*sqrt(n*promedio*(1-promedio))

  extremo_izq <- (sum_Xi - k)/n
  extremo_der <- (sum_Xi + k)/n

  longitud <- extremo_der - extremo_izq

  contiene_p_real <- ( (p_real<=extremo_der) & (p_real >= extremo_izq) )

  return (c(extremo_izq, extremo_der, contiene_p_real, longitud))

  }

  # Estimacion Metodo 2
 estimacion_metodo_2 <- function(muestra, p_real) {

  promedio <- mean(muestra)
  n <- length(muestra)
  quantil <- qnorm(0.975) # z_{1-\alpha/2}

  k <- quantil*sqrt(4*n*promedio + quantil^2 - 4*n*promedio^2)

  extremo_izq <- (2*n*promedio + quantil^2 - k)/(2*(n+quantil^2))
  extremo_der <- (2*n*promedio + quantil^2 + k)/(2*(n+quantil^2))

  longitud <- extremo_der - extremo_izq

  contiene_p_real <- ( (p_real<=extremo_der) & (p_real >= extremo_izq) )

  return (c(extremo_izq, extremo_der, contiene_p_real, longitud))

}

experimento <- function(n,p,k) {

  muestras <- matrix(nrow=k,ncol=n)

  estimacion_1 <- matrix(nrow=k,ncol=4)
  estimacion_2 <- matrix(nrow=k,ncol=4)

  colnames(estimacion_1) <- c("ext_izq", "ext_der", "contiene_p", "long")
  colnames(estimacion_2) <- c("ext_izq", "ext_der", "contiene_p", "long")

  for (i in 1:k){
    muestras[i,] <- generar_muestra(n,p)

    estimacion_1[i,] <- estimacion_metodo_1(muestras[i,], p_real=p)

    estimacion_2[i,] <- estimacion_metodo_2(muestras[i,], p_real=p)
  }

  # Devuelo el promedio del cubrimiento empirico y 
  # el promedio de la estimacion de la longitud

  datos <- data <- matrix(nrow=2, ncol=2)
  datos[1,] <- c(sum(estimacion_1[,3])/k, sum(estimacion_1[,4])/k)  # Metodo 1 (cubrimiento , long)
  datos[2,] <- c(sum(estimacion_2[,3])/k, sum(estimacion_2[,4])/k)  # Metodo 2 (cubrimiento , long)

  return( datos ) 
}

# Grafico del histograma para la proporcion de cubrimiento empirico

graf_exp_prop <- function(exp,metodo,n,p) {
  hist(exp,
  main=paste("Proporcion de cubrimiento empirico", metodo),
  xlab="Proporcion de aciertos",
  ylab="Frecuencia",
  col="chocolate",
  border="brown"
  )

  title(sub = paste(paste("n=",n),paste("y p=",p)))

  abline(v = mean(exp),col = "red",lwd = 3)                       # Add line for mean
}

# Grafico del histograma para la longitudes de IC

graf_exp_long <- function(exp,metodo,n,p) {
  hist(exp,
  main=paste("Longitud de IC", metodo),
  xlab="Longitudes",
  ylab="Frecuencia",
  col="blue",
  border="yellow"
  )

  title(sub = paste(paste("n=",n),paste("y p=",p)))

  abline(v = mean(exp),col = "red",lwd = 3)                       # Add line for mean
}

# Repito el experimento k veces (para el histograma)

k_experimentos <- function(n,p,k){

  exp_metodo_1 <- matrix(nrow=k,ncol=2)
  exp_metodo_2 <- matrix(nrow=k,ncol=2)

  for (i in 1:k){
  exp <- experimento(n[1],p[1],k)

  exp_metodo_1[i,] <- exp[1,]
  exp_metodo_2[i,] <- exp[2,]
}

  graf_exp_prop(exp_metodo_1[,1], "metodo 1", n, p)
  graf_exp_long(exp_metodo_1[,2], "metodo 1", n, p)

  graf_exp_prop(exp_metodo_2[,1], "metodo 2", n, p)
  graf_exp_long(exp_metodo_2[,2], "metodo 2", n, p)

  data <- matrix(nrow=2, ncol=5)

  data[1, ] <- c(n,p,1,mean(exp_metodo_1[,1]), mean(exp_metodo_1[,2]))
  data[2, ] <- c(n,p,2,mean(exp_metodo_2[,1]), mean(exp_metodo_2[,2]))

  return(data)
}


k<-100
n<-c(20, 50, 100)
p<-c(0.10,0.5)


data <- matrix(nrow=12, ncol=7)
colnames(data) <- c("n", "p", "metodo", "ext_izq", "ext_der", "cubrimiento_emp","long_estimada")

# n=20, p=0.1
n_aux <- 20
p_aux <- 0.1
muestra <- generar_muestra(n_aux,p_aux)
data[1,]<-c(n_aux, p_aux, 1, estimacion_metodo_1(muestra,p_aux))
data[2,]<-c(n_aux, p_aux, 2, estimacion_metodo_2(muestra,p_aux))

# n=20, p=0.5
n_aux <- 20
p_aux <- 0.5
muestra <- generar_muestra(n_aux,p_aux)
data[3,]<-c(n_aux, p_aux, 1, estimacion_metodo_1(muestra,p_aux))
data[4,]<-c(n_aux, p_aux, 2, estimacion_metodo_2(muestra,p_aux))

# n=50, p=0.1
n_aux <- 50
p_aux <- 0.1
muestra <- generar_muestra(n_aux,p_aux)
data[5,]<-c(n_aux, p_aux, 1, estimacion_metodo_1(muestra,p_aux))
data[6,]<-c(n_aux, p_aux, 2, estimacion_metodo_2(muestra,p_aux))

# n=50, p=0.5
n_aux <- 50
p_aux <- 0.5
muestra <- generar_muestra(n_aux,p_aux)
data[7,]<-c(n_aux, p_aux, 1, estimacion_metodo_1(muestra,p_aux))
data[8,]<-c(n_aux, p_aux, 2, estimacion_metodo_2(muestra,p_aux))

# n=100, p=0.1
n_aux <- 100
p_aux <- 0.1
muestra <- generar_muestra(n_aux,p_aux)
data[9,]<-c(n_aux, p_aux, 1, estimacion_metodo_1(muestra,p_aux))
data[10,]<-c(n_aux, p_aux, 2, estimacion_metodo_2(muestra,p_aux))

# n=100, p=0.5
n_aux <- 100
p_aux <- 0.5
muestra <- generar_muestra(n_aux,p_aux)
data[11,]<-c(n_aux, p_aux, 1, estimacion_metodo_1(muestra,p_aux))
data[12,]<-c(n_aux, p_aux, 2, estimacion_metodo_2(muestra,p_aux))


print(data)

# Aca en adelante se realiza k veces cada experimento
# Aca guardo la data que va generando para presentarlo como una tabla al final

data <- matrix(nrow=12, ncol=5)
colnames(data) <- c("n", "p", "metodo", "cubrimiento_emp","long_estimada")


# Analizo el cubrimiento empirico con ambos metodos siendo n = 20, p = 0.1
data_aux <- k_experimentos(n[1],p[1], k)

data[1,] <- data_aux[1,]
data[2,] <- data_aux[2,]

# Analizo el cubrimiento empirico con ambos metodos siendo n = 20, p = 0.5
print("Siendo n = 20, p = 0.5")
data_aux <- k_experimentos(n[1],p[2], k)

data[3,] <- data_aux[1,]
data[4,] <- data_aux[2,]

# Analizo el cubrimiento empirico con ambos metodos siendo n = 50, p = 0.1
print("Siendo n = 50, p = 0.1")
data_aux <- k_experimentos(n[2],p[1], k)

data[5,] <- data_aux[1,]
data[6,] <- data_aux[2,]

# Analizo el cubrimiento empirico con ambos metodos siendo n = 50, p = 0.5
print("Siendo n = 50, p = 0.5")
data_aux <- k_experimentos(n[2],p[2], k)

data[7,] <- data_aux[1,]
data[8,] <- data_aux[2,]

# Analizo el cubrimiento empirico con ambos metodos siendo n = 100, p = 0.1
print("Siendo n = 100, p = 0.1" )
data_aux <- k_experimentos(n[3],p[1],k)

data[9,] <- data_aux[1,]
data[10,] <- data_aux[2,]


# Analizo el cubrimiento empirico con ambos metodos siendo n = 100, p = 0.5
print("Siendo n = 100, p = 0.5" )
data_aux <- k_experimentos(n[3],p[2],k)

data[11,] <- data_aux[1,]
data[12,] <- data_aux[2,]

print(data)