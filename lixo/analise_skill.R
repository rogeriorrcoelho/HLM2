install.packages("stopwords")
install.packages("tidytext")
install.packages("wordcloud2")

library(wordcloud2)
library(tidytext)
library(dplyr)
library(stopwords)
library(stringr)

# Ler o arquivo Excel usando a biblioteca 'readxl'
dados <- readxl::read_excel("OPORTUNIDADES.xlsx")

# Concatenar todas as células em uma única string
texto_completo <- paste(dados$skill, collapse = " ")

# Identificar palavras compostas
palavras_compostas <- str_extract_all(texto_completo, "\\b\\w+\\s\\w+\\b")

# Calcular a frequência das palavras-chave
contagem <- table(palavras_compostas)
contagem_ordenada <- sort(contagem, decreasing = TRUE)
df_frequencia <- as.data.frame(contagem_ordenada)
head(df_frequencia)
wordcloud2(data = df_frequencia, size = 0.85)
