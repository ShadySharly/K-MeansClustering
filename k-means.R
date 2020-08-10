# ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////// #
# ATRIBUTOS # ///////////////////////////////////////////////////////////////////////////////////////////////////////// #
# ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////// #

library(ggpubr)
library(cowplot)
library(corrplot)
library(mice)
library(missForest)
library(VIM)
library(Hmisc)
library(factoextra)
library(gridExtra)
library(RColorBrewer)
library(tidyr)
library(scatterplot3d)
library(cluster)
library(GGally)
library(plotly)
library(dplyr)

# Para Test de Normalidad
library(normtest)
library(nortest)
library(moments)

# - code: Numero del c?digo de la muestra
# - clumpThickness: Grosor del grupo (1 - 10)
# - unifCellSize: Tama?o de c?lula uniforme (1 - 10)
# - unifCellShape: Forma de c?lula uniforme (1 - 10)
# - marginalAdhesion: Adhesi?n marginal (1 - 10)
# - epithCellSize: Tama?o de c?lula epitelial (1 - 10)
# - bareNuclei: N?cleos desnudos (1 - 10)
# - blandChromatin: Cromatina suave (1 - 10)
# - normalNucleoli: Nucleolos normales (1 - 10) 
# - mitoses: Mitosis (1 - 10)
# - class: Clase (2 para BENIGNO, 4 para MALIGNO)

# Se crean los nombres que representaran a cada columna, relativos a los par?metros que son de relevancia
# en cada observaci?n.
columns = c("code",
            "clumpThickness",
            "unifCellSize",
            "unifCellShape",
            "marginalAdhesion",
            "epithCellSize",
            "bareNuclei",
            "blandChromatin",
            "normalNucleoli",
            "mitoses",
            "class"
            )

# Se procede a almacenar los datos desde el repositorio web "Breast Cancer Wisconsin" (Original), esto en
# un data frame llamado "df"
url = "https://archive.ics.uci.edu/ml/machine-learning-databases/breast-cancer-wisconsin/breast-cancer-wisconsin.data"
#url = "breast-cancer-wisconsin.data"
df = read.csv(url, header = F, sep=",", col.names = columns)

set.seed(20)
# ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////// #
# PRE-PROCESAMIENTO # ///////////////////////////////////////////////////////////////////////////////////////////////// #
# ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////// #
# 
# Antes de iniciar con el tratamiento de los valores omitidos o missing values, es necesario
# dejar de lado un par de variables que no aportan informaci?n, o que no tiene sentido 
# conservarlas para la implementaci?n y an?lisis del m?todo de las k-medias.

# En primera instancia, el c?digo de cada observaci?n es mas un identificador que no tiene ning?n
# tipo de relaci?n e influencia con las dem?s variables ni el problema de estudio, es solo un
# numero, por lo que se puede descartar
df <- subset(df, select = -c(code))
df.original <- df

# En segunda instancia, para efectos de lo que se busca lograr con el m?todo de las k-medias,
# tampoco tiene sentido conservar la variable "class", pues este m?todo de clustering es del
# tipo "No Supervisado", osea no existe una variable explicativa que categoriza las observaciones
# y las agrupe, de hecho no tendr?a sentido realizar el m?todo si de ante mano ya se conocen los 
# grupos. Posterior a la implementaci?n del m?todo y su respectivo an?lisis, se puede volver a
# esta variable para efectos de comparaci?n.
df <- subset(df, select = -c(class))

#______________________________________________________________________________________________________________________ #
#______________________________________________________________________________________________________________________ #
# I. MISSING VALUES

# Se sabe que el conjunto de datos cuenta con 16 observaciones que presentan missing values para 
# la variable "bareNuclei", denotados por un caracter "?", sin embargo el lenguaje R normalmente
# asocia este tipo de valores con el s?mbolo "NA" al igual que todos los paquetes relativos a los
# missing values, por lo que para trabajar de mejor manera se procede a cambiar los "?" por "NA".
df.n = nrow(df)
df.m = ncol(df)

for (row in 1:df.n) {
  for (col in 1:df.m) {
    if (df[row, col] == "?") {
      df[row, col] <- NA
    }
  }
}

# Debido a que la variable bareNuclei conten?a valores "?" la variable esta clasificada como de tipo
# "character". por lo que es necesario modificarla para que sea del tipo "integer".
df$bareNuclei <- as.integer(df$bareNuclei)

# De acuerdo a la descripci?n de los datos del archivo "names", que viene junto con el conjunto
# datos disponible en el repositorio, existe 1 missing value para la variable "bareNuclei" para 
# 16 observaciones distintas, ahora bien para corroborar que esto es correcto se utiliza el 
# package "mice" que contiene una funci?n conocida como "md.pattern()", que retorna una tabla
# con las variables que presentan missing values y la cantidad.
md.pattern(df)

# Interpretando esta tabla, se entiende que existen 683 observaciones que no presentan
# missing values para ninguna de sus variables, y por otro lado existen 16 observaciones que
# para alguna de sus variables se presenta un missing value, ahora caracterizado por "NA", y
# que ademas estos missing values son para sola una variable, la cual es "bareNuclei".

# Otra manera de verificar esto, es a trav?s del paquete "VIM", para dilucidar esto a trav?s
# de gr?ficos mas agradables a la vista.
  
missing_data <- aggr(
  df, 
  col = mdc(1: 2),
  numbers = TRUE,
  sortVars = TRUE,
  labels = names(df),
  cex.axis = .55, 
  gap = 3,
  ylab = c("Proporci?n de Datos Omitidos", "Patr?n de Datos Omitidos"))

# Para el gr?fico anterior, las barras de color amarillo representan porcentajes de 
# missing values, mientras que las azules representan aquellos que no lo son, por ende se
# puede interpretar que del conjunto de datos total, existe un 97.7% de valores que no son
# missing values, mientras que un 2.3% corresponde a missing values pertenecientes 
# ?nicamente a la variable "bareNuclei".

###########################
# a. Listwise Deletion    #
###########################

# Una de las formas de manejar los valores omitidos, consiste en simplemente eliminar cada
# observaci?n que en sus variables presente uno o mas missing values, ahora bien la simplicidad
# de este m?todo viene perjudicado por el hecho de que el modelo pierde fuerza, ya que se
# pierde informaci?n relevante, ahora bien dependiendo de la raz?n entre el numero de observaciones
# que presentan missing values y el total de observaciones, puede afectar en menor o mayor
# medida la precisi?n del modelo para predecir una variable de estudio. En este caso, raz?n
# de observaciones que se perder?an al aplicar este m?todo corresponde a un 2.3% aproximadamente,
# un numero bastante bajo para considerar este m?todo.
df.listwise <- na.omit(df)
summary(df.listwise)

####################
# b. Imputation    #
####################

# Si se busca minimizar el impacto que tiene la perdida de informaci?n derivada de los m?todos de
# eliminaci?n, resulta mas pertinente emplear los denominados m?todos de imputaci?n, los cuales en
# vez de eliminar observaciones que presentan valores omitidos, tratan de predecir que valor podr?a 
# tomar la variable omitida en cuesti?n a partir de las variables no omitidas, esto a trav?s de la
# generaci?n de modelos de regresi?n lineal utilizando diversas medidas, como la "media" o la
# "mediana" por ejemplo. Existe gran variedad de paquetes provistos para R para realizar esta 
# tarea, dependiendo de la distribuci?n que sigue el conjunto de datos algunos paquetes podr?an ser
# efectivos y otros no, as? tambi?n como la naturaleza de los datos.

### Usando "mice Package" ###
# Este paquete permite generar m?ltiples imputaciones asumiendo que los valores omitidos siguen una
# MAR (Missing at Random), lo que quiere decir que no existe una relaci?n de causalidad que relacione
# un valor omitido, con alguna variable u otra observaci?n, en otras palabras, el valor omitido
# depende ?nicamente de la observaci?n a la que pertenece. A trav?s de la funci?n "mice" se pueden 
# generar distintos conjuntos de datos a partir de un conjunto de datos inicial que contiene
# missing values, estos conjuntos generados difieren precisamente en los valores que son imputados
# para cada missing value. 
miceImp <- mice(df, m = 5, maxit = 50, method = 'pmm', seed = 500)
summary(miceImp)

# Donde "pmm" especifica el m?todo empleado para la imputaci?n, y que corresponde a "Predictive Mean
# Matching", es decir se utiliza la media para predecir los valores imputados.

# Es posible visualizar los valores imputados para cada uno de los conjuntos de datos, donde se 
# puede observar que para cada conjunto se generan distintas combinaciones de valores imputados para
# la variable en cuesti?n.
miceImp$imp$bareNuclei

# As? tambi?n es posible manipular cada uno de los conjuntos por separado.
df.miceImp1 <- complete(miceImp, 1)
df.miceImp2 <- complete(miceImp, 2)
df.miceImp3 <- complete(miceImp, 3)
df.miceImp4 <- complete(miceImp, 4)
df.miceImp5 <- complete(miceImp, 5)

### Usando "missForest Package" ###
# Otra de las formas de realizar imputaciones, consiste en implementar un m?todo basado en el 
# algoritmo de "Random Forest", aplicable para datos no parametricos, osea aquellos que no siguen una
# distribuci?n normal. Para esto se crea un modelo de Random Forest para cada variable, para 
# predecir los missing values de la variable basados en los valores observados.
forestImp <- missForest(df)
summary(forestImp)

# Valores imputados
df.forestImp <- forestImp$ximp
df.forestImp$bareNuclei = round(df.forestImp$bareNuclei)

# Error de Imputaci?n
forestImp$OOBerror

# - Usando "Hmisc Package"
# Hmisc es otro de los paquetes que permite manejar missing values, a trav?s de la imputaci?n de sus
# valores, a trav?s de dos funciones principales impute() y aregImpute().

# La funci?n impute() simplemente imputa el missing value utilizando alg?n m?todo estad?stico a 
# elecci?n, como media, m?ximo, etc. Por defecto esta funci?n emplea la mediana como m?todo de
# imputaci?n
df.miscMeanImp <- df
df.miscMeanImp$bareNuclei <- with(df, impute(bareNuclei, mean))
df.miscMeanImp$bareNuclei = round(df.miscMeanImp$bareNuclei)

df.miscRandImp <- df
df.miscRandImp$bareNuclei <- with(df, impute(bareNuclei, 'random'))

df.miscMedianImp <- df
df.miscMedianImp$bareNuclei <- with(df, impute(bareNuclei))

df.miscMinImp <- df
df.miscMinImp$bareNuclei <- with(df, impute(bareNuclei, min))

df.miscMaxImp <- df
df.miscMaxImp$bareNuclei <- with(df, impute(bareNuclei, max))


#______________________________________________________________________________________________________________________ #
#______________________________________________________________________________________________________________________ #
# II. REDUCCION DE DIMENSIONALIDAD





# ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////// #
# CLUSTERING # //////////////////////////////////////////////////////////////////////////////////////////////////////// #
# ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////// #

# Como bien sabemos, este modelo o algoritmo de agrupamiento consta de tomar un conjunto de datos 
# inicial, y para cada una de las observaciones pertenecientes a este, distribuirlas en un numero K
# de grupos distintos, cada uno de los cuales consta de un "centroide" cuya distancia hacia cada una
# de los puntos u observaciones pertenecientes a ese grupo, es la menor posible. En resumen, el 
# algoritmo consta de los siguientes pasos:

# 1. Definir un total de k centroides al azar.
# 
# 2. Calcular las distancias de cada uno de los puntos de entrada a los k centroides, y asignar cada
#    punto al centroide cuya distancia sea menor.
#   
# 3. Actualizar la posici?n de los k centroides, calculando la posici?n promedio de todos los puntos
#    que pertenecen a cada clase.
#    
# 4. Repetir los pasos 2 y 3 hasta que los centroides no cambien de posici?n y, por lo tanto, las 
#    asignaciones de puntos entre clases no cambie.

#______________________________________________________________________________________________________________________ #
#______________________________________________________________________________________________________________________ #
# I. TEST DE NORMALIDAD

df.current <- df.listwise

###################
# a. Hip?tesis    #
###################

# H0: La muestra proviene de una distribuci?n normal.
# H1: La muestra no proviene de una distribuci?n normal.

################################
# b. Nivel de Significancia    #
################################

# El nivel de significancia que se trabajar? es de 0.05. Alfa=0.05

# Criterio de Decisi?n

# Si p < Alfa Se rechaza Ho
# Si p >= Alfa No se rechaza Ho

# c. Histogramas    #
#####################

hist.var1 <- ggplot(data = df.current, aes(clumpThickness)) + 
  geom_histogram(breaks = seq(1, 10, by = 1), col = "red", aes(fill = ..count.., y = ..density..)) +
  scale_fill_gradient("Cantidad", low = "yellow", high = "red") +
  geom_density(col = 1)

hist.var1

hist.var2 <- ggplot(data = df.current, aes(unifCellSize)) + 
  geom_histogram(breaks = seq(1, 10, by = 1), col = "red", aes(fill = ..count.., y = ..density..)) +
  scale_fill_gradient("Count", low = "yellow", high = "red") +
  geom_density(col = 1)

hist.var2

hist.var3 <- ggplot(data = df.current, aes(unifCellShape)) + 
  geom_histogram(breaks = seq(1, 10, by = 1), col = "red", aes(fill = ..count.., y = ..density..)) +
  scale_fill_gradient("Count", low = "yellow", high = "red") +
  geom_density(col = 1)

hist.var3

hist.var4 <- ggplot(data = df.current, aes(marginalAdhesion)) + 
  geom_histogram(breaks = seq(1, 10, by = 1), col = "red", aes(fill = ..count.., y = ..density..)) +
  scale_fill_gradient("Count", low = "yellow", high = "red") +
  geom_density(col = 1)

hist.var4

hist.var5 <- ggplot(data = df.current, aes(epithCellSize)) + 
  geom_histogram(breaks = seq(1, 10, by = 1), col = "red", aes(fill = ..count.., y = ..density..)) +
  scale_fill_gradient("Count", low = "yellow", high = "red") +
  geom_density(col = 1)

hist.var5

hist.var6 <- ggplot(data = df.current, aes(bareNuclei)) + 
  geom_histogram(breaks = seq(1, 10, by = 1), col = "red", aes(fill = ..count.., y = ..density..)) +
  scale_fill_gradient("Count", low = "yellow", high = "red") +
  geom_density(col = 1)

hist.var6

hist.var7 <- ggplot(data = df.current, aes(blandChromatin)) + 
  geom_histogram(breaks = seq(1, 10, by = 1), col = "red", aes(fill = ..count.., y = ..density..)) +
  scale_fill_gradient("Count", low = "yellow", high = "red") +
  geom_density(col = 1)

hist.var7

hist.var8 <- ggplot(data = df.current, aes(normalNucleoli)) + 
  geom_histogram(breaks = seq(1, 10, by = 1), col = "red", aes(fill = ..count.., y = ..density..)) +
  scale_fill_gradient("Count", low = "yellow", high = "red") +
  geom_density(col = 1)

hist.var8

hist.var9 <- ggplot(data = df.current, aes(mitoses)) + 
  geom_histogram(breaks = seq(1, 10, by = 1), col = "red", aes(fill = ..count.., y = ..density..)) +
  scale_fill_gradient("Count", low = "yellow", high = "red") +
  geom_density(col = 1)

hist.var9

grid.arrange(hist.var1,
             hist.var2,
             hist.var3,
             hist.var4,
             hist.var5,
             hist.var6,
             hist.var7,
             hist.var8,
             hist.var9,
             nrow = 3)

#########################################################
# d. Pruebas de Normalidad con el paquete "normtest"    #
#########################################################

### Prueba de Pearson chi-square ###
# basada en una distribuci?n Ji cuadrado y que corresponde a una prueba de bondad de ajuste.

pearson.test(df.current$clumpThickness)
pearson.test(df.current$unifCellSize)
pearson.test(df.current$unifCellShape)
pearson.test(df.current$marginalAdhesion)
pearson.test(df.current$epithCellSize)
pearson.test(df.current$bareNuclei)
pearson.test(df.current$blandChromatin)
pearson.test(df.current$normalNucleoli)
pearson.test(df.current$mitoses)

### Prueba de Shapiro-Francia ###

sf.test(df.current$clumpThickness)
sf.test(df.current$unifCellSize)
sf.test(df.current$unifCellShape)
sf.test(df.current$marginalAdhesion)
sf.test(df.current$epithCellSize)
sf.test(df.current$bareNuclei)
sf.test(df.current$blandChromatin)
sf.test(df.current$normalNucleoli)
sf.test(df.current$mitoses)

########################################################
# e. Pruebas de Normalidad con el paquete "nortest"    #
########################################################

### Prueba de Jarque Bera ###

jb.norm.test(df.current$clumpThickness)
jb.norm.test(df.current$unifCellSize)
jb.norm.test(df.current$unifCellShape)
jb.norm.test(df.current$marginalAdhesion)
jb.norm.test(df.current$epithCellSize)
jb.norm.test(df.current$bareNuclei)
jb.norm.test(df.current$blandChromatin)
jb.norm.test(df.current$normalNucleoli)
jb.norm.test(df.current$mitoses)

### Prueba de Geary ###
# Usa los valores acumulados muestrales, sus medias y desviaciones est?ndar.

geary.norm.test(df.current$clumpThickness)
geary.norm.test(df.current$unifCellSize)
geary.norm.test(df.current$unifCellShape)
geary.norm.test(df.current$marginalAdhesion)
geary.norm.test(df.current$epithCellSize)
geary.norm.test(df.current$bareNuclei)
geary.norm.test(df.current$blandChromatin)
geary.norm.test(df.current$normalNucleoli)
geary.norm.test(df.current$mitoses)

########################################################
# f. Pruebas de Normalidad con el paquete "moments"    #
########################################################

### Prueba de Agostino ###

agostino.test(df.current$clumpThickness)
agostino.test(df.current$unifCellSize)
agostino.test(df.current$unifCellShape)
agostino.test(df.current$marginalAdhesion)
agostino.test(df.current$epithCellSize)
agostino.test(df.current$bareNuclei)
agostino.test(df.current$blandChromatin)
agostino.test(df.current$normalNucleoli)
agostino.test(df.current$mitoses)

###########################################################
# g. Funciones incluidas en los paquetes b?sicos de R.    #
###########################################################

### Prueba de Shapiro-Wilk ###
# Es m?s poderosa cuando se compara con otras pruebas de normalidad cuando la muestra es peque?a.

shapiro.test(df.current$clumpThickness)
shapiro.test(df.current$unifCellSize)
shapiro.test(df.current$unifCellShape)
shapiro.test(df.current$marginalAdhesion)
shapiro.test(df.current$epithCellSize)
shapiro.test(df.current$bareNuclei)
shapiro.test(df.current$blandChromatin)
shapiro.test(df.current$normalNucleoli)
shapiro.test(df.current$mitoses)

#######################
# h. Conclusiones.    #
#######################

# De acuerdo a los resultados obtenidos en un principio para el m?todo gr?fico, ilustrado en los 
# histogramas obtenidos, ademas de los resultados arrojados por cada uno de los test realizados, se
# puede concluir que existe suficiente evidencia estad?stica para rechazar Ho, por lo que se puede 
# decir que las variables no siguen una distribuci?n normal.

#______________________________________________________________________________________________________________________ #
#______________________________________________________________________________________________________________________ #
# II. DETERMINAR K ?PTIMO

# Ahora bien existe una relaci?n entre el numero de conjuntos a formar, y la WCSS (Within Clusters 
# Summed Squares), es decir, la suma de los cuadrados dentro del cluster, y es que generalmente
# un mayor numero de grupos implica un WCSS menor, lo cual se explica en que la agrupaci?n es mucho
# mas especializada, por otro lado el fin mismo del clustering consta en agrupar un
# conjunto de datos inicial, pero manteniendo un numero considerable de grupos, ya que emplear el
# algoritmo con un k muy elevado hace que este pierda el sentido, por lo que es necesario mantener un
# cierto equilibrio entre el numero de grupos formados, y la WCSS.

# El primer paso consta de determinar el k mas apropiado, que logre optimizar al agrupamiento de
# forma de equilibrar el numero de centroides que corresponde al numero de clusters creados, y la
# WCSS. Para lograr un consenso respecto al k, se emplean 3 m?todos que emplean distintos 
# par?metros en la formaci?n de cada cluster. Para cada uno de los m?todos se define un rango de
# 10 valores de prueba para k.

# Ahora bien, para la utilizaci?n de cada uno de estos m?todo es necesario tener en consideraci?n 
# la distribuci?n bajo la cual se comporta este conjunto de datos, y tal como se pudo concluir en el
# apartado del Test de Normalidad, los datos no siguen una distribuci?n normal, raz?n por la cual no
# resulta recomendable utilizar 

euclidean_dist <- daisy(df.current, metric = "euclidean", stand = FALSE)
gower_dist <- daisy(df.current, metric = "gower", stand = FALSE)
manhattan_dist <- daisy(df.current, metric = "manhattan", stand = FALSE)
distance <- gower_dist
distance_matrix <- as.matrix(distance)

#########################
# a. M?todo del Codo    #
#########################

# El primero de los m?todos se basa principalmente en la suma de los cuadrados dentro de los clusters
# para un determinado rango de valores para k, dibujando una curva continua desde k con valor 0 hasta 
# el m?ximo escogido. Cada transici?n de un k al siguiente significa una variaci?n en la WCSS, al 
# aumentar en uno el valor de k, la idea b?sica es escoger un k cuya transici?n signifique una ca?da
# significante de la WCSS, lo suficiente para compensar el ingreso de un nuevo cluster.

fviz_nbclust(distance_matrix, kmeans, nstart = 25, method = "wss", k.max = 10) + 
  ggtitle("M?todo del Codo") + 
  xlab("k") +
  ylab("WCSS") +
  geom_vline(xintercept = 2, linetype = "dashed", color = "steelblue")

# De acuerdo al gr?fico anterior, en el eje "Y" se tiene la WCSS, mientras que para el eje "X" se
# tiene la cantidad de centroides k, que tambi?n coincide con el numero de grupos resultantes del
# clustering. Como es posible observar, en la medida que el numero de centroides se incrementa, la 
# WCSS cae de manera brusca entre los k 1 y 2, mientras que para los dem?s lo hace de una forma 
# mucho mas pausada, d?ndole a la curva la forma de un "codo". Para escoger el valor optimo de k, la
# idea consiste en encontrar el punto para el cual la WCSS ya no sufre variaciones significantes al
# aumentar el k, lo que significa que la WCSS ganada no es suficientemente significante en 
# comparaci?n con aumentar el k una unidad. Finalmente, para este m?todo el valor para k que 
# resulta ser optimo es un numero de 2 clusters.

##############################
# b. M?todo de la Silueta    #
##############################

# En resumen, la aproximaci?n del promedio de la silueta mide la calidad de un clustering.
# Esto significa que determina que tan bien un objeto cae dentro de un determinado 
# cluster, por consiguiente, una gran promedio en el ancho de la silueta indica un buen
# clustering.

fviz_nbclust(distance_matrix, kmeans, nstart = 25, method = "silhouette", k.max = 10) + 
  ggtitle("Silhouette Method") + 
  xlab("k") +
  ylab("Ancho Promedio de Silueta")

# Observando los resultados del gr?fico anterior, los valores de k para los cuales el ancho
# de silueta promedio se maximiza en orden descendente corresponden a 2, 3 y 9.

#########################################
# c. M?todo de la Brecha Estad?stica    #
#########################################

# Este m?todo compara la WCSS para diferentes valores de k con sus valores esperados bajo una
# distribuci?n con nula referencia de los datos, es decir una distribuci?n sin un clustering tan
# obvio. Este conjunto es generado utilizando las simulaciones de Monte Carlo en el proceso de 
# muestreo. 

# fviz_nbclust(as.matrix(distance), kmeans, nstart = 25, method = "gap_stat", k.max = 10) +
#  ggtitle("M?todo de la Brecha Estadistica") + 
#  xlab("k") +
#  ylab("Brecha")
  
# Finalmente para este m?todo, el k optimo recae en el valor 6.

#______________________________________________________________________________________________________________________ #
#______________________________________________________________________________________________________________________ #
# III. ALGORITMO K-MEDIAS

kmeans2 <- kmeans(distance_matrix, 2, nstart = 25)
kmeans2
kmeans3 <- kmeans(distance_matrix, 3, nstart = 25)
kmeans3
kmeans4 <- kmeans(distance_matrix, 4, nstart = 25)
kmeans4
kmeans5 <- kmeans(distance_matrix, 5, nstart = 25)
kmeans5
kmeans6 <- kmeans(distance_matrix, 6, nstart = 25)
kmeans6

kmeans2.p <- fviz_cluster(
  kmeans2,
  data = distance_matrix,
  geom = "point",
  ellipse.type = "norm",
  main = "Clusters k = 2",
  xlab = "X",
  ylab = "Y",
  palette = "jco",
)
kmeans2.p

kmeans3.p <- fviz_cluster(
  kmeans3,
  data = distance_matrix,
  geom = "point",
  ellipse.type = "norm",
  main = "Clusters k = 3",
  xlab = "X",
  ylab = "Y",
  palette = "jco"
)
kmeans3.p

kmeans4.p <- fviz_cluster(
  kmeans4,
  data = distance_matrix,
  geom = "point",
  ellipse.type = "norm",
  main = "Clusters k = 4",
  xlab = "X",
  ylab = "Y",
  palette = "jco"
)
kmeans4.p

kmeans5.p <- fviz_cluster(
  kmeans5,
  data = distance_matrix,
  geom = "point",
  ellipse.type = "norm",
  main = "Clusters k = 5",
  xlab = "X",
  ylab = "Y",
  palette = "jco"
)
kmeans5.p

kmeans6.p <- fviz_cluster(
  kmeans6,
  data = distance_matrix,
  geom = "point",
  ellipse.type = "norm",
  main = "Clusters k = 6",
  xlab = "X",
  ylab = "Y",
  palette = "jco"
)
kmeans6.p

grid.arrange(kmeans2.p, kmeans3.p, kmeans4.p, kmeans5.p, nrow = 2)


#______________________________________________________________________________________________________________________ #
#______________________________________________________________________________________________________________________ #
# IV. MAPAS DE CALOR

kmeans2.c <- kmeans2$centers
kmeans3.c <- kmeans3$centers
kmeans4.c <- kmeans4$centers
kmeans5.c <- kmeans5$centers
kmeans6.c <- kmeans6$centers

cluster2 <- c(1: 2)
cluster3 <- c(1: 3)
cluster4 <- c(1: 4)
cluster5 <- c(1: 5)
cluster6 <- c(1: 6)

center2_df <- data.frame(cluster2, kmeans2.c)
center3_df <- data.frame(cluster3, kmeans3.c)
center4_df <- data.frame(cluster4, kmeans4.c)
center5_df <- data.frame(cluster5, kmeans5.c)
center6_df <- data.frame(cluster6, kmeans6.c)

hm.palette <- colorRampPalette(rev(brewer.pal(10, 'RdYlGn')), space='Lab')

center2Reshape <- gather(center2_df, caracteristica, valor, 1: 9)
center3Reshape <- gather(center3_df, caracteristica, valor, 1: 9)
center4Reshape <- gather(center4_df, caracteristica, valor, 1: 9)
center5Reshape <- gather(center5_df, caracteristica, valor, 1: 9)
center6Reshape <- gather(center6_df, caracteristica, valor, 1: 9)

ggplot(data = center2Reshape, aes(x = caracteristica, y = cluster2, fill = valor)) +
  scale_y_continuous(breaks = seq(1, 2, by = 1)) +
  geom_tile() +
  ggtitle("Mapa de Calor k = 2") +
  xlab("caracteristica") +
  ylab("cluster") +
  coord_equal() +
  scale_fill_gradientn(colours = hm.palette(90)) +
  theme_classic()

ggplot(data = center3Reshape, aes(x = caracteristica, y = cluster3, fill = valor)) +
  scale_y_continuous(breaks = seq(1, 3, by = 1)) +
  geom_tile() +
  ggtitle("Mapa de Calor k = 3") +
  xlab("caracteristica") +
  ylab("cluster") +
  coord_equal() +
  scale_fill_gradientn(colours = hm.palette(90)) +
  theme_classic()

ggplot(data = center4Reshape, aes(x = caracteristica, y = cluster4, fill = valor)) +
  scale_y_continuous(breaks = seq(1, 4, by = 1)) +
  geom_tile() +
  ggtitle("Mapa de Calor k = 4") +
  xlab("caracteristica") +
  ylab("cluster") +
  coord_equal() +
  scale_fill_gradientn(colours = hm.palette(90)) +
  theme_classic()

ggplot(data = center5Reshape, aes(x = caracteristica, y = cluster5, fill = valor)) +
  scale_y_continuous(breaks = seq(1, 5, by = 1)) +
  geom_tile() +
  ggtitle("Mapa de Calor k = 5") +
  xlab("caracteristica") +
  ylab("cluster") +
  coord_equal() +
  scale_fill_gradientn(colours = hm.palette(90)) +
  theme_classic()

ggplot(data = center6Reshape, aes(x = caracteristica, y = cluster6, fill = valor)) +
  scale_y_continuous(breaks = seq(1, 6, by = 1)) +
  geom_tile() +
  ggtitle("Mapa de Calor k = 6") +
  xlab("caracteristica") +
  ylab("cluster") +
  coord_equal() +
  scale_fill_gradientn(colours = hm.palette(90)) +
  theme_classic()

# ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////// #
# ANALISIS DE RESULTADOS # //////////////////////////////////////////////////////////////////////////////////////////// #
# ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////// #

#______________________________________________________________________________________________________________________ #
#______________________________________________________________________________________________________________________ #
# I. GRAFICOS DE SILUETA

sil2 <- silhouette(kmeans2$cluster, distance)
sil3 <- silhouette(kmeans3$cluster, distance)
sil4 <- silhouette(kmeans4$cluster, distance)
sil5 <- silhouette(kmeans5$cluster, distance)

silohuette.p1 <- fviz_silhouette(sil2)
silohuette.p1
silohuette.p2 <- fviz_silhouette(sil3)
silohuette.p2
silohuette.p3 <- fviz_silhouette(sil4)
silohuette.p3
silohuette.p4 <- fviz_silhouette(sil5)
silohuette.p4

####

# Comparacion de resultados respecto a si la muestra de un tumor efectivamente
# correspondia a uno canceroso o no

# La comparacion se realiza con el dataframe original, que incluye la variable class
df.original.n = nrow(df.original)
df.original.m = ncol(df.original)

for (row in 1:df.original.n) {
  for (col in 1:df.original.m) {
    if (df.original[row, col] == "?") {
      df.original[row, col] <- NA
    }
  }
}

# Se le aplica el mismo preprocesamietno
df.original$bareNuclei <- as.integer(df.original$bareNuclei)
df.original <- na.omit(df.original)
classId <- df.original$class
classId <- sapply(classId, function(x) if(x == 2) "Benigno" else "Canceroso")

#Con k=2
clusterId2 <- kmeans2$cluster
by_clusterId <- data.frame(clusterId2, classId) %>% group_by(clusterId2,classId)
results2 <- by_clusterId %>% summarise(
   cases = n()
)

geomcol2 <- ggplot(results2, aes(clusterId2,cases)) + 
  geom_col(aes(fill=classId), position = "dodge") +
  scale_x_continuous(breaks = seq(1, 2, by = 1)) +
  ggtitle("Gráfico de barras de Clustering de 2 grupos") +
  xlab("Cluster") +
  ylab("Casos") +
  scale_fill_discrete(name = "Clase") +
  annotate(geom = "text", x = 0.79 , y = 22, label = "3%") +
  annotate(geom = "text", x = 1.23 , y = 230, label = "97%") +
  annotate(geom = "text", x = 1.79 , y = 455, label = "95%") +
  annotate(geom = "text", x = 2.23 , y = 40, label = "5%")
geomcol2

#Con k=3
clusterId3 <- kmeans3$cluster
by_clusterId <- data.frame(clusterId3, classId) %>% group_by(clusterId3,classId)
results3 <- by_clusterId %>% summarise(
  cases = n()
)

geomcol3 <- ggplot(results3, aes(clusterId3,cases)) + 
  geom_col(aes(fill=classId), position = "dodge") +
  scale_x_continuous(breaks = seq(1, 3, by = 1)) +
  ggtitle("Gráfico de barras de Clustering de 3 grupos") +
  xlab("Cluster") +
  ylab("Casos") +
  scale_fill_discrete(name = "Clase") +
  annotate(geom = "text", x = 0.79 , y = 28, label = "10%") +
  annotate(geom = "text", x = 1.23 , y = 140, label = "90%") +
  annotate(geom = "text", x = 1.79 , y = 444, label = "98%") +
  annotate(geom = "text", x = 2.23 , y = 25, label = "2%") +
  annotate(geom = "text", x = 3 , y = 118, label = "100%")
geomcol3

#Con k=4
clusterId4 <- kmeans4$cluster
by_clusterId <- data.frame(clusterId4, classId) %>% group_by(clusterId4,classId)
results4 <- by_clusterId %>% summarise(
  cases = n()
)

geomcol4 <- ggplot(results4, aes(clusterId4,cases)) + 
  geom_col(aes(fill=classId), position = "dodge") +
  scale_x_continuous(breaks = seq(1, 4, by = 1)) +
  ggtitle("Gráfico de barras de Clustering de 4 grupos") +
  xlab("Cluster") +
  ylab("Casos") +
  scale_fill_discrete(name = "Clase") +
  annotate(geom = "text", x = 0.79 , y = 445, label = "99%") +
  annotate(geom = "text", x = 1.23 , y = 18, label = "1%") +
  annotate(geom = "text", x = 2 , y = 78, label = "100%") +
  annotate(geom = "text", x = 2.79 , y = 26, label = "14%") +
  annotate(geom = "text", x = 3.23 , y = 99, label = "86%") +
  annotate(geom = "text", x = 3.79 , y = 18, label = "2%") +
  annotate(geom = "text", x = 4.23 , y = 105, label = "98%")
geomcol4

#Con k=5
clusterId5 <- kmeans5$cluster
by_clusterId <- data.frame(clusterId5, classId) %>% group_by(clusterId5,classId)
results5 <- by_clusterId %>% summarise(
  cases = n()
)

geomcol5 <- ggplot(results5, aes(clusterId5,cases)) + 
  geom_col(aes(fill=classId), position = "dodge") +
  scale_x_continuous(breaks = seq(1, 5, by = 1)) +
  ggtitle("Gráfico de barras de Clustering de 5 grupos") +
  xlab("Cluster") +
  ylab("Casos") +
  scale_fill_discrete(name = "Clase") +
  annotate(geom = "text", x = 0.79 , y = 23, label = "14%") +
  annotate(geom = "text", x = 1.23 , y = 88, label = "86%") +
  annotate(geom = "text", x = 2 , y = 76, label = "100%") +
  annotate(geom = "text", x = 3 , y = 383, label = "100%") +
  annotate(geom = "text", x = 3.79 , y = 17, label = "2%") +
  annotate(geom = "text", x = 4.23 , y = 102, label = "98%") +
  annotate(geom = "text", x = 4.79 , y = 73, label = "83%") +
  annotate(geom = "text", x = 5.27 , y = 23, label = "17%")
geomcol5

#Con k=6
clusterId6 <- kmeans6$cluster
by_clusterId <- data.frame(clusterId6, classId) %>% group_by(clusterId6,classId)
results6 <- by_clusterId %>% summarise(
  cases = n()
)

geomcol6 <- ggplot(results6, aes(clusterId6,cases)) + 
  geom_col(aes(fill=classId), position = "dodge") +
  scale_x_continuous(breaks = seq(1, 6, by = 1)) +
  ggtitle("Gráfico de barras de Clustering de 6 grupos") +
  xlab("Cluster") +
  ylab("Casos") +
  scale_fill_discrete(name = "Clase")+
  annotate(geom = "text", x = 1 , y = 385, label = "100%") +
  annotate(geom = "text", x = 2 , y = 75, label = "100%") +
  annotate(geom = "text", x = 2.79 , y = 30, label = "14%") +
  annotate(geom = "text", x = 3.23 , y = 87, label = "86%") +
  annotate(geom = "text", x = 3.79 , y = 18, label = "3%") +
  annotate(geom = "text", x = 4.23 , y = 93, label = "97%") +
  annotate(geom = "text", x = 5 , y = 34, label = "100%") +
  annotate(geom = "text", x = 5.79 , y = 74, label = "83%") +
  annotate(geom = "text", x = 6.32 , y = 28, label = "17%")
geomcol6

grid.arrange(geomcol2, geomcol3, geomcol4, geomcol5, geomcol6, nrow = 3)

########

# Se calcula la media y mediana para los elementos de cada cluster

#Para k=2
clusterInfo2 <- data.frame(clusterId2, df.original) %>% group_by(clusterId2)
clusterInfo2.mean <- clusterInfo2 %>% summarise(
  clumpThickness = mean(clumpThickness),
  unifCellSize = mean(unifCellSize),
  unifCellShape = mean(unifCellShape),
  marginalAdhesion = mean(marginalAdhesion),
  epithCellSize = mean(epithCellSize),
  bareNuclei = mean(bareNuclei),
  blandChromatin = mean(blandChromatin),
  normalNucleoli = mean(normalNucleoli),
  mitoses = mean(mitoses)
)
clusterInfo2.mean

clusterInfo2.median <- clusterInfo2 %>% summarise(
  clumpThickness = median(clumpThickness),
  unifCellSize = median(unifCellSize),
  unifCellShape = median(unifCellShape),
  marginalAdhesion = median(marginalAdhesion),
  epithCellSize = median(epithCellSize),
  bareNuclei = median(bareNuclei),
  blandChromatin = median(blandChromatin),
  normalNucleoli = median(normalNucleoli),
  mitoses = median(mitoses)
)
clusterInfo2.median

#Para k=3
clusterInfo3 <- data.frame(clusterId3, df.original) %>% group_by(clusterId3)
clusterInfo3.mean <- clusterInfo3 %>% summarise(
  clumpThickness = mean(clumpThickness),
  unifCellSize = mean(unifCellSize),
  unifCellShape = mean(unifCellShape),
  marginalAdhesion = mean(marginalAdhesion),
  epithCellSize = mean(epithCellSize),
  bareNuclei = mean(bareNuclei),
  blandChromatin = mean(blandChromatin),
  normalNucleoli = mean(normalNucleoli),
  mitoses = mean(mitoses)
)
clusterInfo3.mean

clusterInfo3.median <- clusterInfo3 %>% summarise(
  clumpThickness = median(clumpThickness),
  unifCellSize = median(unifCellSize),
  unifCellShape = median(unifCellShape),
  marginalAdhesion = median(marginalAdhesion),
  epithCellSize = median(epithCellSize),
  bareNuclei = median(bareNuclei),
  blandChromatin = median(blandChromatin),
  normalNucleoli = median(normalNucleoli),
  mitoses = median(mitoses)
)
clusterInfo3.median

#Para k=4
clusterInfo4 <- data.frame(clusterId4, df.original) %>% group_by(clusterId4)
clusterInfo4.mean <- clusterInfo4 %>% summarise(
  clumpThickness = mean(clumpThickness),
  unifCellSize = mean(unifCellSize),
  unifCellShape = mean(unifCellShape),
  marginalAdhesion = mean(marginalAdhesion),
  epithCellSize = mean(epithCellSize),
  bareNuclei = mean(bareNuclei),
  blandChromatin = mean(blandChromatin),
  normalNucleoli = mean(normalNucleoli),
  mitoses = mean(mitoses)
)
clusterInfo4.mean

clusterInfo4.median <- clusterInfo4 %>% summarise(
  clumpThickness = median(clumpThickness),
  unifCellSize = median(unifCellSize),
  unifCellShape = median(unifCellShape),
  marginalAdhesion = median(marginalAdhesion),
  epithCellSize = median(epithCellSize),
  bareNuclei = median(bareNuclei),
  blandChromatin = median(blandChromatin),
  normalNucleoli = median(normalNucleoli),
  mitoses = median(mitoses)
)
clusterInfo4.median

#Para k=5
clusterInfo5 <- data.frame(clusterId5, df.original) %>% group_by(clusterId5)
clusterInfo5.mean <- clusterInfo5 %>% summarise(
  clumpThickness = mean(clumpThickness),
  unifCellSize = mean(unifCellSize),
  unifCellShape = mean(unifCellShape),
  marginalAdhesion = mean(marginalAdhesion),
  epithCellSize = mean(epithCellSize),
  bareNuclei = mean(bareNuclei),
  blandChromatin = mean(blandChromatin),
  normalNucleoli = mean(normalNucleoli),
  mitoses = mean(mitoses)
)
clusterInfo5.mean

clusterInfo5.median <- clusterInfo5 %>% summarise(
  clumpThickness = median(clumpThickness),
  unifCellSize = median(unifCellSize),
  unifCellShape = median(unifCellShape),
  marginalAdhesion = median(marginalAdhesion),
  epithCellSize = median(epithCellSize),
  bareNuclei = median(bareNuclei),
  blandChromatin = median(blandChromatin),
  normalNucleoli = median(normalNucleoli),
  mitoses = median(mitoses)
)
clusterInfo5.median

#Para k=6
clusterInfo6 <- data.frame(clusterId6, df.original) %>% group_by(clusterId6)
clusterInfo6.mean <- clusterInfo6 %>% summarise(
  clumpThickness = mean(clumpThickness),
  unifCellSize = mean(unifCellSize),
  unifCellShape = mean(unifCellShape),
  marginalAdhesion = mean(marginalAdhesion),
  epithCellSize = mean(epithCellSize),
  bareNuclei = mean(bareNuclei),
  blandChromatin = mean(blandChromatin),
  normalNucleoli = mean(normalNucleoli),
  mitoses = mean(mitoses)
)
clusterInfo6.mean

clusterInfo6.median <- clusterInfo6 %>% summarise(
  clumpThickness = median(clumpThickness),
  unifCellSize = median(unifCellSize),
  unifCellShape = median(unifCellShape),
  marginalAdhesion = median(marginalAdhesion),
  epithCellSize = median(epithCellSize),
  bareNuclei = median(bareNuclei),
  blandChromatin = median(blandChromatin),
  normalNucleoli = median(normalNucleoli),
  mitoses = median(mitoses)
)
clusterInfo6.median


#______________________________________________________________________________________________________________________ #
#______________________________________________________________________________________________________________________ #
# II. INTERPRETACION DEL CLUSTERING

df.current$cluster2 <- as.factor(kmeans2$cluster)
df.current$cluster3 <- as.factor(kmeans3$cluster)
df.current$cluster4 <- as.factor(kmeans4$cluster)
df.current$cluster5 <- as.factor(kmeans5$cluster)

p1 <- ggparcoord(data = df.current, columns = c(1:9), groupColumn = "cluster2", scale = "std") + 
  labs(x = "Caracteristicas", y = "Valor (en unidades de SD)", title = "Clustering k = 2")

ggplotly(p1)

p2 <- ggparcoord(data = df.current, columns = c(1:9), groupColumn = "cluster3", scale = "std") + 
  labs(x = "Caracteristicas", y = "Valor (en unidades de SD)", title = "Clustering k = 3")

ggplotly(p2)

p3 <- ggparcoord(data = df.current, columns = c(1:9), groupColumn = "cluster4", scale = "std") + 
  labs(x = "Caracteristicas", y = "Valor (en unidades de SD)", title = "Clustering k = 4")

ggplotly(p3)

p4 <- ggparcoord(data = df.current, columns = c(1:9), groupColumn = "cluster5", scale = "std") + 
  labs(x = "Caracteristicas", y = "Valor (en unidades de SD)", title = "Clustering k = 5")

ggplotly(p4)

