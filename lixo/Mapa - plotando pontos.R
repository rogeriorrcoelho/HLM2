################################################################################
#               INSTALAÇÃO E CARREGAMENTO DE PACOTES NECESSÁRIOS               #
################################################################################
#Pacotes utilizados
pacotes <- c("plotly","tidyverse","reshape2","knitr","kableExtra","rgl",
             "gghalves","ggdist","tidyquant","car","nlme","lmtest",
             "fastDummies","msm","lmeInfo","jtools","gganimate","ggridges",
             "viridis","hrbrthemes","leaflet","geosphere")

options(rgl.debug = TRUE)

if(sum(as.numeric(!pacotes %in% installed.packages())) != 0){
  instalador <- pacotes[!pacotes %in% installed.packages()]
  for(i in 1:length(instalador)) {
    install.packages(instalador, dependencies = T)
    break()}
  sapply(pacotes, require, character = T) 
} else {
  sapply(pacotes, require, character = T) 
}

# Crie um dataframe de exemplo com as coordenadas de latitude e longitude dos crimes
df_latlon <- read.csv("cpvlatlon.csv", header = TRUE)
df_latlon <- df_latlon[,c("lat","lon")]
df_latlon <- na.omit(df_latlon)
df <- df_latlon
# Criando o data frame resultante
resultante <- df %>%
  group_by_all() %>%
  mutate(ocorrencias = n())
df <- df_latlon
# Crie o mapa interativo com base no OpenStreetMap
mapa <- leaflet() %>%
  addTiles() %>%
  setView(lng =-45.7028 , lat = -23.09412, zoom = 12)

size <- 1
mapa <- mapa %>%
  addCircleMarkers(data = df_latlon, lng = ~lon, lat = ~lat, 
                   radius = ~sqrt(size), color = "red", fill = TRUE, fillOpacity = 0.6,
                   popup = ~paste("Ocorrências:", 1))
mapa

#==================================================
# exemplo 2

# Instale as bibliotecas necessárias, se ainda não estiverem instaladas
install.packages("leaflet")
install.packages("geosphere")

# Carregue as bibliotecas necessárias
library(leaflet)
library(geosphere)

# Crie um dataframe de exemplo com as coordenadas de latitude e longitude e número de ocorrências dos crimes
dados_crimes <- data.frame(
  latitude = c(40.7128, 40.7306, 40.7128, 40.7219, 40.7589),
  longitude = c(-74.0060, -74.0650, -74.0060, -74.0460, -73.9850),
  ocorrencias = c(10, 5, 7, 3, 8)
)

# Defina o raio de 1 quilômetro em metros
raio <- 1000

# Crie uma matriz com as coordenadas em radianos 
# ERRO!!!!
coord_rad <- matrix(geosphere::rad(dados_crimes[c("longitude", "latitude")]), ncol = 2)

# Realize a clusterização utilizando o algoritmo DBSCAN
clusters <- dbscan(coord_rad, eps = raio, minPts = 2)

# Adicione o rótulo de cluster e o número de ocorrências ao dataframe original
dados_crimes$cluster <- clusters$cluster
dados_crimes$size <- dados_crimes$ocorrencias

# Crie o mapa interativo com base no OpenStreetMap
mapa <- leaflet() %>%
  addTiles() %>%
  setView(lng = -74.0060, lat = 40.7128, zoom = 12)

# Adicione os pontos de crimes ao mapa, com tamanho proporcional ao número de ocorrências
mapa <- mapa %>%
  addCircleMarkers(data = dados_crimes, lng = ~longitude, lat = ~latitude, 
                   radius = ~sqrt(size), color = "red", fill = TRUE, fillOpacity = 0.6,
                   popup = ~paste("Ocorrências:", ocorrencias))

# Exiba o mapa
mapa

# código corrigido
# Instale as bibliotecas necessárias, se ainda não estiverem instaladas
install.packages("leaflet")
install.packages("geosphere")

# Carregue as bibliotecas necessárias
library(leaflet)
library(geosphere)

# Crie um dataframe de exemplo com as coordenadas de latitude e longitude e número de ocorrências dos crimes
dados_crimes <- data.frame(
  latitude = c(40.7128, 40.7306, 40.7128, 40.7219, 40.7589),
  longitude = c(-74.0060, -74.0650, -74.0060, -74.0460, -73.9850),
  ocorrencias = c(10, 5, 7, 3, 8)
)

# Defina o raio de 1 quilômetro em metros
raio <- 1000

# Crie uma matriz com as coordenadas em radianos
coord_rad <- sapply(dados_crimes[c("longitude", "latitude")], function(coord) {
  destPoint(coord, raio/1000, 90)$lon
})

# Realize a clusterização utilizando o algoritmo DBSCAN
clusters <- dbscan(coord_rad, eps = raio, minPts = 2)

# Adicione o rótulo de cluster e o número de ocorrências ao dataframe original
dados_crimes$cluster <- clusters$cluster
dados_crimes$size <- dados_crimes$ocorrencias

# Crie o mapa interativo com base no OpenStreetMap
mapa <- leaflet() %>%
  addTiles() %>%
  setView(lng = -74.0060, lat = 40.7128, zoom = 12)

# Adicione os pontos de crimes ao mapa, com tamanho proporcional ao número de ocorrências
mapa <- mapa %>%
  addCircleMarkers(data = dados_crimes, lng = ~longitude, lat = ~latitude, 
                   radius = ~sqrt(size), color = "red", fill = TRUE, fillOpacity = 0.6,
                   popup = ~paste("Ocorrências:", ocorrencias))

# Exiba o mapa
mapa

# tentativa 3=====================================================
# Instale as bibliotecas necessárias, se ainda não estiverem instaladas
install.packages("leaflet")
install.packages("geosphere")

# Carregue as bibliotecas necessárias
library(leaflet)
library(geosphere)

# Crie um dataframe de exemplo com as coordenadas de latitude e longitude e número de ocorrências dos crimes
dados_crimes <- data.frame(
  latitude = c(40.7128, 40.7306, 40.7128, 40.7219, 40.7589),
  longitude = c(-74.0060, -74.0650, -74.0060, -74.0460, -73.9850),
  ocorrencias = c(10, 5, 7, 3, 8)
)

# Defina o raio de 1 quilômetro em metros
raio <- 1000

# Converta as coordenadas em radianos usando a função destPoint()
coord_rad <- geosphere::destPoint(dados_crimes[, c("longitude", "latitude")], 
                                  b = raio/1000, 
                                  a = seq(0, 360, length.out = nrow(dados_crimes)))

# Realize a clusterização utilizando o algoritmo DBSCAN
clusters <- dbscan(coord_rad, eps = raio, minPts = 2)

# Adicione o rótulo de cluster e o número de ocorrências ao dataframe original
dados_crimes$cluster <- clusters$cluster
dados_crimes$size <- dados_crimes$ocorrencias

# Crie o mapa interativo com base no OpenStreetMap
mapa <- leaflet() %>%
  addTiles() %>%
  setView(lng = -74.0060, lat = 40.7128, zoom = 12)

# Adicione os pontos de crimes ao mapa, com tamanho proporcional ao número de ocorrências
mapa <- mapa %>%
  addCircleMarkers(data = dados_crimes, lng = ~longitude, lat = ~latitude, 
                   radius = ~sqrt(size), color = "red", fill = TRUE, fillOpacity = 0.6,
                   popup = ~paste("Ocorrências:", ocorrencias))

# Exiba o mapa
mapa


#=========================================================
# mais uma tentativa kkk
# Instale as bibliotecas necessárias, se ainda não estiverem instaladas
install.packages("leaflet")
install.packages("geosphere")

# Carregue as bibliotecas necessárias
library(leaflet)
library(geosphere)

# Crie um dataframe de exemplo com as coordenadas de latitude e longitude e número de ocorrências dos crimes
dados_crimes <- data.frame(
  latitude = c(40.7128, 40.7306, 40.7128, 40.7219, 40.7589),
  longitude = c(-74.0060, -74.0650, -74.0060, -74.0460, -73.9850),
  ocorrencias = c(10, 5, 7, 3, 8)
)

# Defina o raio de 1 quilômetro em metros
raio <- 1000

# Converta as coordenadas em radianos usando a função destPoint() e apply()
coord_rad <- apply(dados_crimes[, c("longitude", "latitude")], 1, function(coord) {
  geosphere::destPoint(p = coord, b = raio/1000, a = 0)$lon
})

# Realize a clusterização utilizando o algoritmo DBSCAN
clusters <- dbscan(coord_rad, eps = raio, minPts = 2)

# Adicione o rótulo de cluster e o número de ocorrências ao dataframe original
dados_crimes$cluster <- clusters$cluster
dados_crimes$size <- dados_crimes$ocorrencias

# Crie o mapa interativo com base no OpenStreetMap
mapa <- leaflet() %>%
  addTiles() %>%
  setView(lng = -74.0060, lat = 40.7128, zoom = 12)

# Adicione os pontos de crimes ao mapa, com tamanho proporcional ao número de ocorrências
mapa <- mapa %>%
  addCircleMarkers(data = dados_crimes, lng = ~longitude, lat = ~latitude, 
                   radius = ~sqrt(size), color = "red", fill = TRUE, fillOpacity = 0.6,
                   popup = ~paste("Ocorrências:", ocorrencias))

# Exiba o mapa
mapa








# Instale as bibliotecas necessárias, se ainda não estiverem instaladas
install.packages("leaflet")
install.packages("geosphere")

# Carregue as bibliotecas necessárias
library(leaflet)
library(geosphere)

# Crie um dataframe de exemplo com as coordenadas de latitude e longitude e número de ocorrências dos crimes
dados_crimes <- data.frame(
  latitude = c(40.7128, 40.7306, 40.7128, 40.7219, 40.7589),
  longitude = c(-74.0060, -74.0650, -74.0060, -74.0460, -73.9850),
  ocorrencias = c(100, 5, 7, 38, 8)
)

# Defina o raio de 1 quilômetro em metros
raio <- 50000

# Converta as coordenadas em radianos usando a função distGeo() e sapply()
coord_rad <- sapply(1:nrow(dados_crimes), function(i) {
  distGeo(c(dados_crimes[i, "longitude"], dados_crimes[i, "latitude"]), c(dados_crimes[i, "longitude"] + 1, dados_crimes[i, "latitude"]))$lon
})

# Realize a clusterização utilizando o algoritmo DBSCAN
clusters <- dbscan(coord_rad, eps = raio, minPts = 2)

# Adicione o rótulo de cluster e o número de ocorrências ao dataframe original
dados_crimes$cluster <- clusters$cluster
dados_crimes$size <- dados_crimes$ocorrencias

# Crie o mapa interativo com base no OpenStreetMap
mapa <- leaflet() %>%
  addTiles() %>%
  setView(lng = -74.0060, lat = 40.7128, zoom = 12)

# Adicione os pontos de crimes ao mapa, com tamanho proporcional ao número de ocorrências
mapa <- mapa %>%
  addCircleMarkers(data = dados_crimes, lng = ~longitude, lat = ~latitude, 
                   radius = ~sqrt(size), color = "red", fill = TRUE, fillOpacity = 0.6,
                   popup = ~paste("Ocorrências:", ocorrencias))

# Exiba o mapa
mapa


###################### mais uma tentativa, agora com versão nova do R
# Instale as bibliotecas necessárias, se ainda não estiverem instaladas

# Carregue as bibliotecas necessárias
library(leaflet)
library(geosphere)

# Crie um dataframe de exemplo com as coordenadas de latitude e longitude e número de ocorrências dos crimes
dados_crimes <- data.frame(
  latitude = c(40.7128, 40.7306, 40.7128, 40.7219, 40.7589),
  longitude = c(-74.0060, -74.0650, -74.0060, -74.0460, -73.9850),
  ocorrencias = c(10, 5, 7, 3, 8)
)

# Defina o raio de 1 quilômetro em metros
raio <- 1000

# Converta as coordenadas em radianos usando a função destPoint() e lapply()
coord_rad <- lapply(1:nrow(dados_crimes), function(i) {
  coord <- dados_crimes[i, c("longitude", "latitude")]
  dest <- destPoint(p = geosphere::as.geosphere(coord), dist = raio / 1000, b = seq(0, 360, length.out = 10))
  rad <- geosphere::as.data.frame(geosphere::as.geosphere(dest))
  rad$lon <- geosphere::rad(rad$lon)
  rad$lat <- geosphere::rad(rad$lat)
  rad
})

# Realize a clusterização utilizando o algoritmo DBSCAN
clusters <- dbscan(do.call(rbind, coord_rad), eps = raio, minPts = 2)

# Adicione o rótulo de cluster e o número de ocorrências ao dataframe original
dados_crimes$cluster <- clusters$cluster
dados_crimes$size <- dados_crimes$ocorrencias

# Crie o mapa interativo com base no OpenStreetMap
mapa <- leaflet() %>%
  addTiles() %>%
  setView(lng = -74.0060, lat = 40.7128, zoom = 12)

# Adicione os pontos de crimes ao mapa, com tamanho proporcional ao número de ocorrências
mapa <- mapa %>%
  addCircleMarkers(data = dados_crimes, lng = ~longitude, lat = ~latitude, 
                   radius = ~sqrt(size), color = "red", fill = TRUE, fillOpacity = 0.6,
                   popup = ~paste("Ocorrências:", ocorrencias))

# Exiba o mapa
mapa

