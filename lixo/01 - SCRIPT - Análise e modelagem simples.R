################################################################################
#               INSTALAÇÃO E CARREGAMENTO DE PACOTES NECESSÁRIOS               #
################################################################################
#Pacotes utilizados
pacotes <- c("plotly","tidyverse","reshape2","knitr","kableExtra","rgl",
             "gghalves","ggdist","tidyquant","car","nlme","lmtest",
             "fastDummies","msm","lmeInfo","jtools","gganimate","ggridges",
             "viridis","hrbrthemes")

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

#Algoritmo para determinação dos erros-padrão das variâncias no componente de
#efeitos aleatórios
#ATENÇÃO: A função abaixo é plenamente funcional para modelos do tipo HLM2
#e HLM3, desde que estimados pelo pacote nlme

stderr_nlme <- function(model){
  if(base::class(model) != "lme"){
    base::message("Use a lme object model from nlme package")
    stop()}
  resume <- base::summary(model)
  if(base::length(base::names(model$groups))==1){
    m.type <- "HLM2"
  } else if(base::length(base::names(model$groups))==2){
    m.type <- "HLM3"
  }
  if(m.type == "HLM2"){
    vcov_matrix <- model$apVar
    logs_sd_re <- base::attr(vcov_matrix,"Pars")
    if(base::length(logs_sd_re)==2){
      stderr_tau00 <- msm::deltamethod(~exp(x1)^2,logs_sd_re,vcov_matrix)
      stderr_sigma <- msm::deltamethod(~exp(x2)^2,logs_sd_re,vcov_matrix)
      results <- base::data.frame(`RE Components`=base::c("Var(v0j)","Var(e)"),
                                  `Variance Estimatives`= base::c(base::exp(logs_sd_re)[[1]]^2,
                                                                  base::exp(logs_sd_re[[2]])^2),
                                  `Std Err.`=base::c(stderr_tau00,
                                                     stderr_sigma),
                                  z=base::c(base::exp(logs_sd_re)[[1]]^2/stderr_tau00,
                                            base::exp(logs_sd_re[[2]])^2/stderr_sigma),
                                  `p-value`=base::round(stats::pnorm(q=base::c(base::exp(logs_sd_re)[[1]]^2/stderr_tau00,
                                                                               base::exp(logs_sd_re[[2]])^2/stderr_sigma),
                                                                     lower.tail=F)*2,3))
      return(results)
    }
    else{
      stderr_tau00 <- msm::deltamethod(~exp(x1)^2,logs_sd_re,vcov_matrix)
      stderr_tau01 <- msm::deltamethod(~exp(x2)^2,logs_sd_re,vcov_matrix)
      stderr_sigma <- msm::deltamethod(~exp(x4)^2,logs_sd_re,vcov_matrix)
      results <- base::data.frame(Components=base::c("Var(v0j)","Var(v1j)","Var(e)"),
                                  Estimatives= base::c(base::exp(logs_sd_re)[[1]]^2,
                                                       base::exp(logs_sd_re[[2]])^2,
                                                       base::exp(logs_sd_re[[4]])^2),
                                  Std_Err=base::c(stderr_tau00,
                                                  stderr_tau01,
                                                  stderr_sigma),
                                  z=base::c(base::exp(logs_sd_re)[[1]]^2/stderr_tau00,
                                            base::exp(logs_sd_re[[2]])^2/stderr_tau01,
                                            base::exp(logs_sd_re[[4]])^2/stderr_sigma),
                                  `p-value`=base::round(stats::pnorm(q=base::c(base::exp(logs_sd_re)[[1]]^2/stderr_tau00,
                                                                               base::exp(logs_sd_re[[2]])^2/stderr_tau01,
                                                                               base::exp(logs_sd_re[[4]])^2/stderr_sigma),
                                                                     lower.tail=F)*2,3))
      return(results)
    }
  }
  if(m.type == "HLM3"){
    vcov_matrix <- model$apVar
    logs_sd_re <-  base::attr(vcov_matrix,"Pars")
    if(base::length(logs_sd_re) == 3){
      stderr_tau_r000 <- msm::deltamethod(~exp(x1)^2,logs_sd_re,vcov_matrix)
      stderr_tau_u000 <- msm::deltamethod(~exp(x2)^2,logs_sd_re,vcov_matrix)
      stderr_sigma <- msm::deltamethod(~exp(x3)^2,logs_sd_re,vcov_matrix)
      results <- base::data.frame(Components=base::c("Var(t00k)","Var(v0jk)","Var(e)"),
                                  Estimatives=base::c(base::exp(logs_sd_re)[[2]]^2,
                                                      base::exp(logs_sd_re)[[1]]^2,
                                                      base::exp(logs_sd_re)[[3]]^2),
                                  Std_Err=base::c(stderr_tau_u000,
                                                  stderr_tau_r000,
                                                  stderr_sigma),
                                  z=base::c(base::exp(logs_sd_re)[[2]]^2/stderr_tau_u000,
                                            base::exp(logs_sd_re)[[1]]^2/stderr_tau_r000,
                                            base::exp(logs_sd_re)[[3]]^2/stderr_sigma),
                                  `p-value`=base::round(stats::pnorm(q=base::c(base::exp(logs_sd_re)[[2]]^2/stderr_tau_u000,
                                                                               base::exp(logs_sd_re)[[1]]^2/stderr_tau_r000,
                                                                               base::exp(logs_sd_re)[[3]]^2/stderr_sigma),
                                                                     lower.tail=F)*2,3))
      return(results)
    } 
    else{
      stderr_tau_r000 <- msm::deltamethod(~exp(x1)^2,logs_sd_re,vcov_matrix)
      stderr_tau_r100 <- msm::deltamethod(~exp(x2)^2,logs_sd_re,vcov_matrix)
      stderr_tau_u000 <- msm::deltamethod(~exp(x4)^2,logs_sd_re,vcov_matrix)
      stderr_tau_u100 <- msm::deltamethod(~exp(x5)^2,logs_sd_re,vcov_matrix)
      stderr_sigma <- msm::deltamethod(~exp(x7)^2,logs_sd_re,vcov_matrix)
      results <- base::data.frame(`RE_Components`=base::c("Var(t00k)","Var(t10k)",
                                                          "Var(v0jk)","Var(v1jk)",
                                                          "Var(e)"),
                                  `Variance Estimatives`=base::c(base::exp(logs_sd_re)[[4]]^2,
                                                                 base::exp(logs_sd_re)[[5]]^2,
                                                                 base::exp(logs_sd_re)[[1]]^2,
                                                                 base::exp(logs_sd_re)[[2]]^2,
                                                                 base::exp(logs_sd_re)[[7]]^2),
                                  `Std Err.`=base::c(stderr_tau_u000,
                                                     stderr_tau_u100,
                                                     stderr_tau_r000,
                                                     stderr_tau_r100,
                                                     stderr_sigma),
                                  z=base::c(base::exp(logs_sd_re)[[4]]^2/stderr_tau_u000,
                                            base::exp(logs_sd_re)[[5]]^2/stderr_tau_u100,
                                            base::exp(logs_sd_re)[[1]]^2/stderr_tau_r000,
                                            base::exp(logs_sd_re)[[2]]^2/stderr_tau_r100,
                                            base::exp(logs_sd_re)[[7]]^2/stderr_sigma),
                                  `p-value`=base::round(stats::pnorm(q=base::c(base::exp(logs_sd_re)[[4]]^2/stderr_tau_u000,
                                                                               base::exp(logs_sd_re)[[5]]^2/stderr_tau_u100,
                                                                               base::exp(logs_sd_re)[[1]]^2/stderr_tau_r000,
                                                                               base::exp(logs_sd_re)[[2]]^2/stderr_tau_r100,
                                                                               base::exp(logs_sd_re)[[7]]^2/stderr_sigma),
                                                                     lower.tail=F)*2,3))
      return(results)
    }
  }
}


################################################################################
#                      DESCRIÇÃO E EXPLORAÇÃO DO DATASET                       #
################################################################################

#Carregamento da base de dados
dados_f_pedidos <- read.csv("Fato_Pedido.csv", header = TRUE)
dados_f_pedidos$data_ref <- as.Date(dados_f_pedidos$Dt_Pedido, format = "%d/%m/%Y")
dados_f_pedidos$data_ref <- format(dados_f_pedidos$data_ref, "%m-%Y")

dados_f_pedidos <- dados_f_pedidos[, c("id_pedido", "id_cliente", "id_cidade", "id_produto", "Dt_Pedido", "Preço_Venda", "Preço_Compra")]


dados_d_cliente <- read.csv("Dominio_Cliente.csv", header = TRUE)
dados_d_cliente <- dados_d_cliente[, c("id_cliente", "id_sexo")]

#juntando os dados
fato_pedidos <- merge(dados_f_pedidos, dados_d_cliente, by = "id_cliente")

#Visualização da base de dados
fato_pedidos %>%
  kable() %>%
  kable_styling(bootstrap_options = "striped", 
                full_width = F, 
                font_size = 12)

#ajustando as colunas com seus tipos corretos
fato_pedidos$id_cliente <-as.factor(fato_pedidos$id_cliente)
fato_pedidos$id_pedido <-as.factor(fato_pedidos$id_pedido)
fato_pedidos$id_produto <-as.factor(fato_pedidos$id_produto)
fato_pedidos$id_cidade <-as.factor(fato_pedidos$id_cidade)

fato_pedidos$Preço_Venda <- gsub("[^0-9\\.\\,]", "", fato_pedidos$Preço_Venda) # remove caracteres não numéricos
fato_pedidos$Preço_Venda <- gsub(",", ".", fato_pedidos$Preço_Venda) # substitui vírgulas por pontos
fato_pedidos$Preço_Venda <- as.numeric(fato_pedidos$Preço_Venda) # converte para numérico

fato_pedidos$Preço_Compra <- gsub("[^0-9\\.\\,]", "", fato_pedidos$Preço_Compra) # remove caracteres não numéricos
fato_pedidos$Preço_Compra <- gsub(",", ".", fato_pedidos$Preço_Compra) # substitui vírgulas por pontos
fato_pedidos$Preço_Compra <- as.numeric(fato_pedidos$Preço_Compra) # converte para numérico

fato_pedidos$id_sexo <-as.factor(fato_pedidos$id_sexo)

#Estatísticas descritivas
summary(fato_pedidos)


#valor médio dos pedidos por cidade 
valor_por_pedido <- aggregate(Preço_Venda ~ id_pedido, data = fato_pedidos, FUN = sum)
#juntando os dados data frames
fato_pedidos <- merge(fato_pedidos, valor_por_pedido, by = "id_pedido")
valor_por_pedido_por_cidade <- aggregate(Preço_Venda.y ~ id_cidade, data = fato_pedidos, FUN = mean)

fato_pedidos <- fato_pedidos %>%
  group_by(id_cidade) %>%
  mutate(preco_medio_pedido_por_cidade = mean(Preço_Venda.y))

fato_pedidos %>%
  kable() %>%
  kable_styling(bootstrap_options = "striped",
                full_width = F,
                font_size = 12)

#Estudo sobre o desbalanceamento dos dados
fato_pedidos %>% 
  group_by(id_cidade) %>% 
  summarise(quantidade = n()) %>% 
  kable() %>%
  kable_styling(bootstrap_options = "striped", 
                full_width = F, 
                font_size = 12)

fato_pedidos <- fato_pedidos %>%
  rename(Preco_Pedido = Preço_Venda.y) %>%
  rename(Preco_Venda = Preço_Venda.x) 
fato_pedidos <- fato_pedidos %>%  
  rename(Preco_Compra = Preço_Compra)


#Valor do ticket de compra médio por cliente
#===========================================
fato_pedidos <- fato_pedidos %>%
  group_by(id_cliente) %>%
  mutate(preco_medio_pedido_por_cliente = mean(Preco_Pedido))

#gráfico de distribuição clientes x valor médio de compra
#verificar por CEP ... bairro onde moram!?!?!

#grava o data frame em arquivo
save(fato_pedidos, file = "TCC_dados.RData")


###########################################################################
# GLM - Regressão Simples
# Qual % de influência da data (data_ref = mm/aaaa) no preço do pedido
# Yhat = intercepto + B.X(i)
###########################################################################
df <- fato_pedidos[, c("Dt_Pedido", "Preco_Pedido")]
df$data_ref <- as.Date(fato_pedidos$Dt_Pedido, format = "%m/%d/%Y")
df$mes_ano <- format(as.Date(df$data_ref, format = "%m/%d/%y"), "%Y-%m")
df <- distinct(df)
df <- df[,c("data_ref","Preco_Pedido")]
df$num_dt <- as.numeric(df$data_ref)
ggplotly(
  ggplot(df, aes(x = num_dt, y = Preco_Pedido)) +
  geom_point() +
  labs(x = "Mês e Ano", y = "Preço") +
  ggtitle("Gráfico de Preços")
)

glimpse(df)
summary(df)

ggplotly(
  ggplot(df, aes(x = num_dt, y = Preco_Pedido)) +
    geom_point(color ="blue" , size = 1) +
    geom_smooth(aes(color = "fitted values"),
                method = "lm", formula = y ~ x, se = F, size = 1 ) +
    labs(x = "Mês e Ano", 
         y = "Preço",
         title = paste("R2:",
                       round(((cor(df$Preco_Pedido, df$num_dt))^2),4))) +
    scale_color_manual("legenda:",
                       values = "grey50") +
    theme_classic()
)
modelo_precodata <- lm(formula= Preco_Pedido ~ num_dt ,
                         data = df)
summary(modelo_precodata)
##########################################################################

###########################################################################
# GLM - Regressão Simples (OLS linear regression)
# Qual % de influência do preço de compra no preço de venda
# Yhat = intercepto + B.X(i)
###########################################################################
df <- fato_pedidos[, c("Preco_Compra", "Preco_Venda","id_cidade")]
df <- distinct(df)

ggplotly(
  ggplot(df, aes(x = Preco_Compra, y = Preco_Venda, color = id_cidade)) +
    geom_point() +
    labs(x = "Compra", y = "Venda", color="id_cidade") +
    theme_minimal() +
    ggtitle("Gráfico de Preços")
)

glimpse(df)
summary(df)

ggplotly(
  ggplot(df, aes(x = Preco_Compra, y = Preco_Venda)) +
    geom_point(color ="blue" , size = 0.7) +
    geom_smooth(aes(color = "fitted values"),
                method = "lm", formula = y ~ x, se = F, size = 0.5 ) +
    labs(x = "Compra", 
         y = "Venda",
         title = paste("R2:",
                       round(((cor(df$Preco_Venda, df$Preco_Compra))^2),4))) +
    scale_color_manual("legenda:",
                       values = "gray70") +
    theme_classic()
)
##########################################################################
#Modelagem Regressão Linear Simples => OLS linear regression
#Estimando o modelo x=compra y=venda

modelo_compravenda <- lm(formula= Preco_Venda ~ Preco_Compra ,
                         data = df)
summary(modelo_compravenda)

summ(modelo_compravenda, confint = T, digits = 4, ci.width = .95)
export_summs(modelo_compravenda, scale = F, digits = 4)

loglik <- logLik(modelo_compravenda)
loglik
###########################################################################
# GLM - Regressão Simples (OLS linear regression)
# Qual % de influência do IPCA no preço de venda
# Yhat = intercepto + B.X(i)
###########################################################################
dados_ipca <- read.csv("ipca.csv", header = TRUE)
dados_ipca$ipca <- gsub("[^0-9\\.\\,]", "", dados_ipca$ipca) # remove caracteres não numéricos
dados_ipca$ipca <- gsub(",", ".", dados_ipca$ipca) # substitui vírgulas por pontos
dados_ipca$ipca <- as.numeric(dados_ipca$ipca) # converte para numérico

dados_ipca$dtref <- substr(dados_ipca$dtref, start = nchar(dados_ipca$dtref) - 6, stop = nchar(dados_ipca$dtref))


df_aux <- fato_pedidos[, c("Dt_Pedido", "Preco_Venda")]
df_aux$Dt_Pedido <- as.Date(df_aux$Dt_Pedido, format = "%m/%d/%Y")
df_aux$Dt_Pedido <- as.character(df_aux$Dt_Pedido)
df_aux$Dt_Pedido <- sub("^(\\d{4})-(\\d{2})-(\\d{2})$", "\\3-\\2-\\1", df_aux$Dt_Pedido)
df_aux$Dt_Pedido <- substr(df_aux$Dt_Pedido, start = nchar(df_aux$Dt_Pedido) - 6, stop = nchar(df_aux$Dt_Pedido))
df_aux$dtref <- df_aux$Dt_Pedido
#juntando os dados
df_aux <- merge(dados_ipca, df_aux, by = "dtref")
df_aux <- distinct(df_aux)

modelo_ipcavenda <- lm(formula= Preco_Venda ~ ipca ,
                         data = df_aux)
summary(modelo_ipcavenda)

summ(modelo_ipcavenda, confint = T, digits = 4, ci.width = .95)
export_summs(modelo_ipcavenda, scale = F, digits = 4)

ggplotly(
  ggplot(df_aux, aes(x = ipca, y = Preco_Venda)) +
    geom_point(color ="blue" , size = 0.7) +
    geom_smooth(aes(color = "fitted values"),
                method = "lm", formula = y ~ x, se = F, size = 0.5 ) +
    labs(x = "ipca", 
         y = "Venda",
         title = paste("R2:",
                       round(((cor(df_aux$Preco_Venda, df_aux$ipca))^2),4))) +
    scale_color_manual("legenda:",
                       values = "gray70") +
    theme_classic()
)

#####################################################################
# FIM => OLS linear regression
#####################################################################




#Exploração visual do desempenho médio 
#Corrigir o data frame com distinct ... esse gráfico está apresentando valores
#duplicados por o data frame está com as informações de agrupamentos em cada linha

fato_pedidos %>%
  ggplot() +
  geom_point(aes(x = id_cidade, y = Preco_Pedido),color = "orange", alpha = 0.5, size = 4) +
  geom_line(aes(x = id_cidade, y = preco_medio_pedido_por_cliente, 
                group = 1, color = "Preço Médio do Pedido por Cliente"), size = 1.5) +
  scale_colour_viridis_d() +
  labs(x = "Cidade",
       y = "Preço Médio") +
  theme(legend.title = element_blank(),
        panel.border = element_rect(NA),
        panel.grid = element_line("grey"),
        panel.background = element_rect("white"),
        legend.position = "bottom",
        axis.text.x = element_text(angle = 90))

#Boxplot da variável dependente (preço do pedido)
dados_preco_pedidos <- distinct(fato_pedidos[, c("id_cliente", "Preco_Pedido")])
ggplotly(
  ggplot(dados_preco_pedidos, aes(x = "", y = Preco_Pedido)) +
    geom_boxplot(fill = "deepskyblue",    # cor da caixa
                 alpha = 0.7,             # transparência
                 color = "black",         # cor da borda
                 outlier.colour = "red",  # cor dos outliers
                 outlier.shape = 15,      # formato dos marcadores dos outliers
                 outlier.size = 2.5) +    # tamanho dos marcadores dos outliers
    geom_jitter(width = 0.1, alpha = 0.5, size = 1.3, color = "darkorchid") +
    labs(y = "Preço dos Pedidos") +
    theme(panel.background = element_rect("white"),
          panel.grid = element_line("grey95"),
          panel.border = element_rect(NA),
          legend.position="none",
          plot.title = element_text(size=15)) +
    ggtitle("Box Plot dos Preços dos Pedidos - todos os pedidos") +
    xlab("")
)

#Boxplot da variável dependente (preço do pedido) por cidade
dados_preco_pedidos_cidade <- distinct(fato_pedidos[, c("id_cidade", "Preco_Pedido")])
ggplotly(
  ggplot(dados_preco_pedidos_cidade, aes(x = id_cidade, y = Preco_Pedido)) +
    geom_boxplot(fill = "deepskyblue",    # cor da caixa
                 alpha = 0.7,             # transparência
                 color = "black",         # cor da borda
                 outlier.colour = "red",  # cor dos outliers
                 outlier.shape = 15,      # formato dos marcadores dos outliers
                 outlier.size = 2.5) +    # tamanho dos marcadores dos outliers
    geom_jitter(width = 0.1, alpha = 0.5, size = 1.3, color = "darkorchid") +
    labs(y = "Desempenho") +
    theme(panel.background = element_rect("white"),
          panel.grid = element_line("grey95"),
          panel.border = element_rect(NA),
          legend.position="none",
          plot.title = element_text(size=15)) +
    ggtitle("Boxplot da variável 'Preço' por 'Cidade'") +
    xlab("")
)

#Kernel density estimation (KDE) - função densidade de probabilidade da
#variável dependente (preço do pedido), com histograma
ggplotly(
  ggplot(fato_pedidos, aes(x = Preco_Pedido)) +
    geom_density(aes(x = Preco_Pedido), 
                 position = "identity", color = "black", size = 1) +
    geom_histogram(aes(y = ..density..), color = "white", fill = "deepskyblue",
                   bins = 30) +
    theme_classic()
)

#Exploração visual do desempenho médio
fato_pedidos %>%
  group_by(id_cidade) %>%
  mutate(preco_medio_i = mean(Preco_Pedido, na.rm = TRUE)) %>% 
  ggplot() +
  geom_point(aes(x = id_cidade, y = preco_medio_i),color = "orange", alpha = 0.5, size = 4) +
  geom_line(aes(x = id_cidade, y = preco_medio_i, 
                group = 1, color = "Preço Médio por Cidade"), size = 1.5) +
  scale_colour_viridis_d() +
  labs(x = "Cidade",
       y = "Preço Médio") +
  theme(legend.title = element_blank(),
        panel.border = element_rect(NA),
        panel.grid = element_line("grey"),
        panel.background = element_rect("white"),
        legend.position = "bottom",
        axis.text.x = element_text(angle = 90))

