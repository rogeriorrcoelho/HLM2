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

# Função para agrupar coordenadas próximas
agrupar_coordenadas <- function(lat, lon, distancia_maxima) {
  coords <- cbind(lon, lat) # Coordenadas em formato [lon, lat]
  n <- nrow(coords) # Número de coordenadas
  
  # Matriz para armazenar as distâncias entre as coordenadas
  distancias <- matrix(0, nrow = n, ncol = n)
  
  # Preencher a matriz de distâncias
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      distancias[i, j] <- geosphere::distHaversine(coords[i, ], coords[j, ]) / 1000 # Distância em km
      distancias[j, i] <- distancias[i, j] # A matriz é simétrica
    }
  }
  
  # Vetor para armazenar os grupos de coordenadas próximas
  grupos <- rep(NA, n)
  
  # Identificar os grupos
  grupo_atual <- 1
  for (i in 1:n) {
    if (is.na(grupos[i])) {
      grupo <- which(distancias[i, ] <= distancia_maxima)
      grupos[grupo] <- grupo_atual
      grupo_atual <- grupo_atual + 1
    }
  }
  
  return(grupos)
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
df <- resultante
df_latlon<-resultante
# Crie o mapa interativo com base no OpenStreetMap
mapa <- leaflet() %>%
  addTiles() %>%
  setView(lng =-45.7028 , lat = -23.09412, zoom = 12)

size <- resultante$ocorrencias
mapa <- mapa %>%
  addCircleMarkers(data = df_latlon, lng = ~lon, lat = ~lat, 
                   radius = ~sqrt(size), color = "red", fill = TRUE, fillOpacity = 0.6,
                   popup = ~paste("Ocorrências:", size))
#mapadados_crimes$size <- dados_crimes$ocorrencias
mapa


#agrupando coordenadas para refazer o mapa
df<- df[,c("lat","lon")]

# Agrupando as coordenadas próximas
grupos <- agrupar_coordenadas(df$lat, df$lon, distancia_maxima = 5) # Distância máxima em quilômetros

# Adicionando a informação do grupo ao data frame original
df$grupo <- grupos

# Exibindo o data frame com a informação do grupo
print(df)

