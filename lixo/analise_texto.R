install.packages("stopwords")
install.packages("tidytext")
install.packages("wordcloud2")

library(wordcloud2)
library(tidytext)
library(dplyr)
library(stopwords)

#abordagem básica para a extração de palavras-chave de textos em português

# Carregar os dados do arquivo CSV
data <- NULL
data <- read.csv("LinkedIn.csv", header = TRUE)

# Pré-processamento dos textos
stpwrd <- stopwords::stopwords("pt")
stpwrd_en <- stopwords::stopwords("en")
stpwrd <- c(stpwrd, stpwrd_en, "and", "or", "ter", "with", "in", "off", "to")

preprocessed_data <- data %>%
  mutate(descritivo = tolower(descritivo)) %>%
  unnest_tokens(palavra, descritivo, token = "words") %>%
  filter(!(palavra %in% stpwrd))

preprocessed_data_compostas <- data %>%
  mutate(descritivo = tolower(descritivo)) %>%
  unnest_tokens(word, input = descritivo) %>%
  filter(!word %in% stpwrd) %>%
  mutate(bigram = paste(word, lead(word), sep = " ")) %>%
  filter(!is.na(bigram))

preprocessed_data_skill <- data %>%
  mutate(skill = tolower(skill)) %>%
  unnest_tokens(palavra, skill, token = "words") %>%
  filter(!(palavra %in% stpwrd))

# Identificação das palavras-chave
keywords <- preprocessed_data %>%
  group_by(empresa, cargo) %>%
  count(palavra, sort = TRUE) %>%
  ungroup()

# Identificação das palavras-chave compostas
keywords_compostas <- preprocessed_data %>%
  group_by(empresa, cargo) %>%
  mutate(bigram = paste(palavra, lead(palavra), sep = " ")) %>%
  filter(!is.na(bigram)) %>%
  count(bigram, sort = TRUE)

keywords_skill <- preprocessed_data_skill %>%
  group_by(empresa, cargo) %>%
  count(palavra, sort = TRUE) %>%
  ungroup()

# Exibir as palavras-chave e palavras-chave compostas
head(keywords)
head(keywords_compostas)
head(keywords_skill)

# Calcular a frequência das palavras-chave
freq_table <- preprocessed_data %>%
  count(palavra, sort = TRUE)

freq_table2 <- preprocessed_data_compostas %>%
  count(bigram, sort = TRUE)

freq_table3 <- preprocessed_data_skill %>%
  count(palavra, sort = TRUE)

# Exibir as principais palavras-chave
head(freq_table)
head(freq_table2)
head(freq_table3)

freq_table2a <- freq_table2 %>%
  rename(palavra = bigram)

# Unir as tabelas de frequências
combined_freq <- bind_rows(freq_table, freq_table2a, freq_table3)

# Ordenar as frequências em ordem decrescente
combined_freq <- combined_freq %>% arrange(desc(n))

wordcloud2(data = combined_freq, size = 0.7)

wordcloud2(data = freq_table, size = 0.7)
wordcloud2(data = freq_table2, size = 0.7)
wordcloud2(data = freq_table3, size = 0.7)

