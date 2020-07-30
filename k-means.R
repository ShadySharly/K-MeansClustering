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
library(mi)
library(factoextra)
library(gridExtra)
library(RColorBrewer)
library(tidyr)
library(scatterplot3d)

# - code: Numero del código de la muestra
# - clumpThickness: Grosor del grupo (1 - 10)
# - unifCellSize: Tamaño de célula uniforme (1 - 10)
# - unifCellShape: Forma de célula uniforme (1 - 10)
# - marginalAdhesion: Adhesión marginal (1 - 10)
# - epithCellSize: Tamaño de célula epitelial (1 - 10)
# - bareNuclei: Núcleos desnudos (1 - 10)
# - blandChromatin: Cromatina suave (1 - 10)
# - normalNucleoli: Nucleolos normales (1 - 10) 
# - mitoses: Mitosis (1 - 10)
# - class: Clase (2 para BENIGNO, 4 para MALIGNO)

# Se crean los nombres que representaran a cada columna, relativos a los parámetros que son de relevancia
# en cada observación.
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
df = read.csv(url, header = F, sep=",", col.names = columns)

# ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////// #
# PRE-PROCESAMIENTO # ///////////////////////////////////////////////////////////////////////////////////////////////// #
# ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////// #
# 
# Antes de iniciar con el tratamiento de los valores omitidos o missing values, es necesario
# dejar de lado un par de variables que no aportan información, o que no tiene sentido 
# conservarlas para la implementación y análisis del método de las k-medias.

# En primera instancia, el código de cada observación es mas un identificador que no tiene ningún
# tipo de relación e influencia con las demás variables ni el problema de estudio, es solo un
# numero, por lo que se puede descartar
df <- subset(df, select = -c(code))

# En segunda instancia, para efectos de lo que se busca lograr con el método de las k-medias,
# tampoco tiene sentido conservar la variable "class", pues este método de clustering es del
# tipo "No Supervisado", osea no existe una variable explicativa que categoriza las observaciones
# y las agrupe, de hecho no tendría sentido realizar el método si de ante mano ya se conocen los 
# grupos. Posterior a la implementación del método y su respectivo análisis, se puede volver a
# esta variable para efectos de comparación.
df <- subset(df, select = -c(class))

#______________________________________________________________________________________________________________________ #
# I. MISSING VALUES

# Se sabe que el conjunto de datos cuenta con 16 observaciones que presentan missing values para 
# la variable "bareNuclei", denotados por un caracter "?", sin embargo el lenguaje R normalmente
# asocia este tipo de valores con el símbolo "NA" al igual que todos los paquetes relativos a los
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

# Debido a que la variable bareNuclei contenía valores "?" la variable esta clasificada como de tipo
# "character". por lo que es necesario modificarla para que sea del tipo "integer".
df$bareNuclei <- as.integer(df$bareNuclei)

# De acuerdo a la descripción de los datos del archivo "names", que viene junto con el conjunto
# datos disponible en el repositorio, existe 1 missing value para la variable "bareNuclei" para 
# 16 observaciones distintas, ahora bien para corroborar que esto es correcto se utiliza el 
# package "mice" que contiene una función conocida como "md.pattern()", que retorna una tabla
# con las variables que presentan missing values y la cantidad.
md.pattern(df)

# Interpretando esta tabla, se entiende que existen 683 observaciones que no presentan
# missing values para ninguna de sus variables, y por otro lado existen 16 observaciones que
# para alguna de sus variables se presenta un missing value, ahora caracterizado por "NA", y
# que ademas estos missing values son para sola una variable, la cual es "bareNuclei".

# Otra manera de verificar esto, es a través del paquete "VIM", para dilucidar esto a través
# de gráficos mas agradables a la vista
mice_plot <- aggr(df, 
                  col = c('navyblue','yellow'),
                  numbers = TRUE, 
                  sortVars = TRUE,
                  labels = names(df), 
                  cex.axis = .55,
                  gap = 3, 
                  ylab = c("Missing data","Pattern")
                  )

# Para el gráfico anterior, las barras de color amarillo representan porcentajes de 
# missing values, mientras que las azules representan aquellos que no lo son, por ende se
# puede interpretar que del conjunto de datos total, existe un 97.7% de valores que no son
# missing values, mientras que un 2.3% corresponde a missing values pertenecientes 
# únicamente a la variable "bareNuclei".

# a. Listwise Deletion
# 
# Una de las formas de manejar los valores omitidos, consiste en simplemente eliminar cada
# observación que en sus variables presente uno o mas missing values, ahora bien la simplicidad
# de este método viene perjudicado por el hecho de que el modelo pierde fuerza, ya que se
# pierde información relevante, ahora bien dependiendo de la razón entre el numero de observaciones
# que presentan missing values y el total de observaciones, puede afectar en menor o mayor
# medida la precisión del modelo para predecir una variable de estudio. En este caso, razón
# de observaciones que se perderían al aplicar este método corresponde a un 2.3% aproximadamente,
# un numero bastante bajo para considerar este método.
df.listwise <- na.omit(df)

# b. Pairwise Deletion

# c. Imputation
#
# Si se busca minimizar el impacto que tiene la perdida de información derivada de los métodos de
# eliminación, resulta mas pertinente emplear los denominados métodos de imputación, los cuales en
# vez de eliminar observaciones que presentan valores omitidos, tratan de predecir que valor podría 
# tomar la variable omitida en cuestión a partir de las variables no omitidas, esto a través de la
# generación de modelos de regresión lineal utilizando diversas medidas, como la "media" o la
# "mediana" por ejemplo. Existe gran variedad de paquetes provistos para R para realizar esta 
# tarea, dependiendo de la distribución que sigue el conjunto de datos algunos paquetes podrían ser
# efectivos y otros no, así también como la naturaleza de los datos.

# - Usando "mice Package"
# Este paquete permite generar múltiples imputaciones asumiendo que los valores omitidos siguen una
# MAR (Missing at Random), lo que quiere decir que no existe una relación de causalidad que relacione
# un valor omitido, con alguna variable u otra observación, en otras palabras, el valor omitido
# depende únicamente de la observación a la que pertenece. A través de la función "mice" se pueden 
# generar distintos conjuntos de datos a partir de un conjunto de datos inicial que contiene
# missing values, estos conjuntos generados difieren precisamente en los valores que son imputados
# para cada missing value. 
miceImp <- mice(df, m = 5, maxit = 50, method = 'pmm', seed = 500)
summary(miceImp)

# Donde "pmm" especifica el método empleado para la imputación, y que corresponde a "Predictive Mean
# Matching", es decir se utiliza la media para predecir los valores imputados.

# Es posible visualizar los valores imputados para cada uno de los conjuntos de datos, donde se 
# puede observar que para cada conjunto se generan distintas combinaciones de valores imputados para
# la variable en cuestión.
miceImp$imp$bareNuclei

# Así también es posible manipular cada uno de los conjuntos por separado.
df.miceImp1 <- complete(miceImp)


# Ahora bien, como se generan 5 conjuntos de datos distintos en este caso, una buena practica
# consiste en combinar los resultados para obtener una predicción mas confiable.


# - Usando "missForest Package"
# Otra de las formas de realizar imputaciones, consiste en implementar un método basado en el 
# algoritmo de "Random Forest", aplicable para datos no parametricos, osea aquellos que no siguen una
# distribución normal. Para esto se crea un modelo de Random Forest para cada variable, para 
# predecir los missing values de la variable basados en los valores observados.
forestImp <- missForest(df)

# Valores imputados
df.forestImp <- forestImp$ximp

# Error de Imputación
forestImp$OOBerror

# - Usando "Hmisc Package"
# Hmisc es otro de los paquetes que permite manejar missing values, a través de la imputación de sus
# valores, a través de dos funciones principales impute() y aregImpute().

# La función impute() simplemente imputa el missing value utilizando algún método estadístico a 
# elección, como media, máximo, etc. Por defecto esta función emplea la mediana como método de
# imputación
df.miscMeanImp <- df
df.miscMeanImp$bareNuclei <- with(df, impute(bareNuclei, mean))

df.miscRandImp <- df
df.miscRandImp$bareNuclei <- with(df, impute(bareNuclei, 'random'))

df.miscMedianImp <- df
df.miscMedianImp$bareNuclei <- with(df, impute(bareNuclei))

df.miscMinImp <- df
df.miscMinImp$bareNuclei <- with(df, impute(bareNuclei, min))

df.miscMaxImp <- df
df.miscMaxImp$bareNuclei <- with(df, impute(bareNuclei, max))

# - Usando "mi Package"
mi_data <- mi(df, seed = 335)
summary(mi_data)
#______________________________________________________________________________________________________________________ #
# II. REDUCCION DE DIMENSIONALIDAD

apply(df.listwise, 2, var)
apply(df.forestImp, 2, var)

# ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////// #
# CLUSTERING # //////////////////////////////////////////////////////////////////////////////////////////////// #
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
# 3. Actualizar la posición de los k centroides, calculando la posición promedio de todos los puntos
#    que pertenecen a cada clase.
#    
# 4. Repetir los pasos 2 y 3 hasta que los centroides no cambien de posición y, por lo tanto, las 
#    asignaciones de puntos entre clases no cambie.

#______________________________________________________________________________________________________________________ #
# I. DETERMINAR K PTIMO

# Ahora bien existe una relación entre el numero de conjuntos a formar, y la WCSS (Within Clusters 
# Summed Squares), es decir, la suma de los cuadrados dentro del cluster, y es que generalmente
# un mayor numero de grupos implica un WCSS menor, lo cual se explica en que la agrupación es mucho
# mas especializada, por otro lado el fin mismo del clustering consta en agrupar un
# conjunto de datos inicial, pero manteniendo un numero considerable de grupos, ya que emplear el
# algoritmo con un k muy elevado hace que este pierda el sentido, por lo que es necesario mantener un
# cierto equilibrio entre el numero de grupos formados, y la WCSS.

# El primer paso consta de determinar el k mas apropiado, que logre optimizar al agrupamiento de
# forma de equilibrar el numero de centroides que corresponde al numero de clusters creados, y la
# WCSS. Para esto se aplica una técnica, que permite visualmente tener una mejor idea del valor
# o conjunto de valores k para los cuales resulta optimo aplicar el algoritmo, correspondiente al 
# llamado "Elbow Method" o "Método del Codo", el cual en pocas palabras busca seleccionar numero de 
# grupos ideal optimizando el WCSS.

# Entonces se procede a crear un vector que almacenara las WCSS para un rango determinado de valores
# de k de prueba, e identificar cual de estos es mas conveniente utilizar. En este caso particular,
# se crea un vector de 20 valores, correspondientes a la WCSS considerando que el clustering es
# realizado considerando desde 1 grupo (donde la WCSS debería ser mayor), hasta 20 grupos (donde la
# WCSS debería ser menor).

wcss <- vector()
for(i in 1:20){
  wcss[i] <- sum(kmeans(df.listwise, i)$withinss)
}

ggplot() + geom_point(aes(x = 1:20, y = wcss), color = 'blue') + 
  geom_line(aes(x = 1:20, y = wcss), color = 'blue') + 
  ggtitle("Método del Codo") + 
  xlab('Cantidad de Centroides k') + 
  ylab('WCSS')

# De acuerdo al gráfico anterior, en el eje "Y" se tiene la WCSS, mientras que para el eje "X" se
# tiene la cantidad de centroides k, que también coincide con el numero de grupos resultantes del
# clustering. Como es posible observar, en la medida que el numero de centroides se incrementa, la 
# WCSS cae de manera brusca entre los k 1 y 2, mientras que para los demás lo hace de una forma 
# mucho mas pausada, dándole a la curva la forma de un "codo". Para escoger el valor optimo de k, la
# idea consiste en encontrar el punto para el cual la WCSS ya no sufre variaciones significantes al
# aumentar el k, lo que significa que la WCSS ganada no es suficientemente significante en 
# comparación con aumentar el k una unidad. Esto resulta mas o menos evidente entre los puntos
# correspondientes los valores de k 2 y 5, y para los cuales se aplica el modelo a continuación.

kmeans2 <- kmeans(df.listwise, 2, iter.max = 1000, nstart = 10)
kmeans3 <- kmeans(df.listwise, 3, iter.max = 1000, nstart = 10)
kmeans4 <- kmeans(df.listwise, 4, iter.max = 1000, nstart = 10)
kmeans5 <- kmeans(df.listwise, 5, iter.max = 1000, nstart = 10)

kmeans2.p <- fviz_cluster(
  kmeans2,
  data = df.listwise,
  geom = "point",
  ellipse.type = "norm",
  main = "Clusters k = 2",
  xlab = "X",
  ylab = "Y",
  palette = "Set2"
)
kmeans2.p

kmeans3.p <- fviz_cluster(
  kmeans3,
  data = df.listwise,
  geom = "point",
  ellipse.type = "norm",
  main = "Clusters k = 3",
  xlab = "X",
  ylab = "Y",
  palette = "Set2"
)
kmeans3.p

kmeans4.p <- fviz_cluster(
  kmeans4,
  data = df.listwise,
  geom = "point",
  ellipse.type = "norm",
  main = "Clusters k = 4",
  xlab = "X",
  ylab = "Y",
  palette = "Set2"
)
kmeans4.p

kmeans5.p <- fviz_cluster(
  kmeans5,
  data = df.listwise,
  geom = "point",
  ellipse.type = "norm",
  main = "Clusters k = 5",
  xlab = "X",
  ylab = "Y",
  palette = "Set2"
)
kmeans5.p

grid.arrange(kmeans2.p, kmeans3.p, kmeans4.p, kmeans5.p, nrow = 2)

# DENDOGRAMAS 

dend2 <- hcut(df.listwise, k = 2, hc_method = "complete")
dend3 <- hcut(df.listwise, k = 3, hc_method = "complete")
dend4 <- hcut(df.listwise, k = 4, hc_method = "complete")
dend5 <- hcut(df.listwise, k = 5, hc_method = "complete")

dend2.p <- fviz_dend(
  dend2,
  show_labels = FALSE,
  rect = TRUE,
  main = "Dendograma k = 2",
  ylab = "Altura",
  ggtheme = theme_grey()
)
dend2.p

dend3.p <- fviz_dend(
  dend3, 
  show_labels = FALSE,
  rect = TRUE,
  main = "Dendograma k = 3",
  ylab = "Altura",
  ggtheme = theme_grey()
)
dend3.p

dend4.p <- fviz_dend(
  dend4, 
  show_labels = FALSE, 
  rect = TRUE,
  main = "Dendograma k = 4",
  ylab = "Altura",
  ggtheme = theme_grey()
)
dend4.p

dend5.p <- fviz_dend(
  dend5, show_labels = FALSE,
  rect = TRUE,
  main = "Dendograma k = 5",
  ylab = "Altura",
  ggtheme = theme_grey()
)
dend5.p

# MAPAS DE CALOR

kmeans2$size
kmeans3$size
kmeans4$size
kmeans5$size

kmeans2.c <- kmeans2$centers
kmeans3.c <- kmeans3$centers
kmeans4.c <- kmeans4$centers
kmeans5.c <- kmeans5$centers

cluster2 <- c(1: 2)
cluster3 <- c(1: 3)
cluster4 <- c(1: 4)
cluster5 <- c(1: 5)

center2_df <- data.frame(cluster2, kmeans2.c)
center3_df <- data.frame(cluster3, kmeans3.c)
center4_df <- data.frame(cluster4, kmeans4.c)
center5_df <- data.frame(cluster5, kmeans5.c)

hm.palette <- colorRampPalette(rev(brewer.pal(10, 'RdYlGn')), space='Lab')

center2Reshape <- gather(center2_df, caracteristica, valor, clumpThickness: mitoses)
center3Reshape <- gather(center3_df, caracteristica, valor, clumpThickness: mitoses)
center4Reshape <- gather(center4_df, caracteristica, valor, clumpThickness: mitoses)
center5Reshape <- gather(center5_df, caracteristica, valor, clumpThickness: mitoses)

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

ggplot(data = kmeans2 , aes(x = dim_1, y = dim_2)) +
  geom_point(aes(color = cluster2)) + 
  theme_bw()


scatterplot3d(x = resultados$dim_1,
              y = resultados$dim_2,
              z = resultados$dim_3,
              pch = 20, color = colores, cex.lab = 0.8,
              grid = TRUE, box = FALSE)
legend("bottom", legend = levels(resultados$numero),
       col = colores, pch = 16, 
       inset = -0.23, xpd = TRUE, horiz = TRUE)

