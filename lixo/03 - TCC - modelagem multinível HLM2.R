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

# Carregamento da base de dados
# Montar o dataframe com as observações abaixo
#
# Para o modelo OLS (valor médio pedido ~ sal_med)
# id_cliente | id_cidade | valor pedido | habitantes | IDHM | PIB | sal_med 
#
# Para o modelo HLM2
# id_cliente | id_cidade | valor pedido | habitantes | IDHM | PIB | sal_med

df_cidade <- read.csv("cidade.csv", header = TRUE)

df_cidade$habitantes <- gsub("[^0-9\\.\\,]", "", df_cidade$habitantes) # remove caracteres não numéricos
df_cidade$habitantes <- gsub(",", ".", df_cidade$habitantes) # substitui vírgulas por pontos
df_cidade$habitantes <- as.numeric(df_cidade$habitantes) # converte para numérico
df_cidade$pib <- gsub("[^0-9\\.\\,]", "", df_cidade$pib) # remove caracteres não numéricos
df_cidade$pib <- gsub(",", ".", df_cidade$pib) # substitui vírgulas por pontos
df_cidade$pib <- as.numeric(df_cidade$pib) # converte para numérico
df_cidade$idhm <- gsub("[^0-9\\.\\,]", "", df_cidade$idhm) # remove caracteres não numéricos
df_cidade$idhm <- gsub(",", ".", df_cidade$idhm) # substitui vírgulas por pontos
df_cidade$idhm <- as.numeric(df_cidade$idhm) # converte para numérico
df_cidade$sal_med <- gsub("[^0-9\\.\\,]", "", df_cidade$sal_med) # remove caracteres não numéricos
df_cidade$sal_med <- gsub(",", ".", df_cidade$sal_med) # substitui vírgulas por pontos
df_cidade$sal_med <- as.numeric(df_cidade$sal_med) # converte para numérico

pedidos <- read.csv("TCC_Fato_Pedido.csv", header = TRUE)
dados_f_pedidos <- pedidos[,c("id_pedido", "id_cliente", "id_cidade", "id_produto", "Dt_Pedido", "Preço_Venda", "Total_Pedido")]
dados_f_pedidos$id_produto <- as.factor(dados_f_pedidos$id_produto)
dados_f_pedidos$id_cidade <- as.factor(dados_f_pedidos$id_cidade)
dados_f_pedidos$id_cliente <- as.factor(dados_f_pedidos$id_cliente)
colnames(dados_f_pedidos)[6] <- "vl_produto"
colnames(dados_f_pedidos)[7] <- "vl_pedido"

dados_f_pedidos$vl_produto <- gsub("[^0-9\\.\\,]", "", dados_f_pedidos$vl_produto) # remove caracteres não numéricos
dados_f_pedidos$vl_produto <- gsub(",", ".", dados_f_pedidos$vl_produto) # substitui vírgulas por pontos
dados_f_pedidos$vl_produto <- as.numeric(dados_f_pedidos$vl_produto) # converte para numérico
dados_f_pedidos$vl_pedido <- gsub("[^0-9\\.\\,]", "", dados_f_pedidos$vl_pedido) # remove caracteres não numéricos
dados_f_pedidos$vl_pedido <- gsub(",", ".", dados_f_pedidos$vl_pedido) # substitui vírgulas por pontos
dados_f_pedidos$vl_pedido <- as.numeric(dados_f_pedidos$vl_pedido) # converte para numérico

dados_f_pedidos <- dados_f_pedidos %>%
  group_by(id_pedido) %>%
  mutate(qtd_produtos = n())

#NOVO REGRESSAO PRODUTOS
#=======================
dados_f_produtos <- pedidos %>%
  group_by(id_produto,id_sexo) %>%
  mutate(qtd_clientes = n(), qtd_total_prod = sum(QTD))


dados_d_cliente <- read.csv("Dominio_Cliente.csv", header = TRUE)
dados_d_cliente <- dados_d_cliente[, c("id_cliente", "id_sexo")]
clientes <- dados_d_cliente
clientes$id_cliente <- as.character(clientes$id_cliente)

#ajustando os clientes com sexo = ? no dataset
clientes <- clientes %>% 
  mutate(id_sexo = ifelse(id_sexo == "?", "F", id_sexo))
clientes$id_sexo <- as.factor(clientes$id_sexo)
summary(clientes)

#juntando os dados
fato_pedidos <- merge(dados_f_pedidos, df_cidade, by = "id_cidade")
summary(fato_pedidos)
fato_pedidos <- merge(fato_pedidos, clientes, by = "id_cliente")
summary(fato_pedidos)
fato_pedidos$id_pedido <- as.factor(fato_pedidos$id_pedido)
summary(fato_pedidos)

#juntando os dados NOVO REGRESSAO PRODUTO
#========================================
fato_produtos <- merge(dados_f_produtos, df_cidade, by = "id_cidade")
summary(fato_produtos)
fato_produtos <- merge(fato_produtos, clientes, by = "id_cliente")
summary(fato_produtos)
fato_produtos$id_produtos <- as.factor(fato_produtos$id_produto)
summary(fato_produtos)

# data frame final
fato_pedidos <- fato_pedidos[,c("id_pedido", "id_cliente", "id_sexo", "id_cidade", "vl_pedido", "qtd_produtos", "habitantes", "pib", "idhm", "sal_med")]
fato_pedidos <- distinct(fato_pedidos)

# Substituir os valores da id_cidade
fato_pedidos <- fato_pedidos %>%
  mutate(id_cidade = case_when(
    id_cidade == "CPV" ~ "Cidade 1",
    id_cidade == "SJC" ~ "Cidade 2",
    id_cidade == "TTE" ~ "Cidade 3",
    id_cidade == "TMB" ~ "Cidade 4",
    TRUE ~ id_cidade
  ))

#Visualização da base de dados 
fato_pedidos %>%
  kable() %>%
  kable_styling(bootstrap_options = "striped", 
                full_width = F, 
                font_size = 12)

fato_pedidos <- na.omit(fato_pedidos)
#Estatísticas descritivas
summary(fato_pedidos)

#Estudo sobre o desbalanceamento dos dados
fato_pedidos %>% 
  group_by(id_cidade) %>% 
  summarise(qtd_pedidos = n()) %>% 
  kable() %>%
  kable_styling(bootstrap_options = "striped", 
                full_width = F, 
                font_size = 12)

#Estudo sobre o desbalanceamento dos dados
fato_pedidos %>% 
  group_by(qtd_produtos) %>% 
  summarise(num_pedidos = n()) %>% 
  kable() %>%
  kable_styling(bootstrap_options = "striped", 
                full_width = F, 
                font_size = 12)

#Pedido médio dos clientes por cidade
fato_pedidos %>%
  group_by(id_cidade) %>%
  summarise(`Valor médio pedido` = mean(vl_pedido, na.rm = T)) %>%
  kable() %>%
  kable_styling(bootstrap_options = "striped",
                full_width = F,
                font_size = 12)

#Qtd média de produtos no pedido por cidade
fato_pedidos %>%
  group_by(id_cidade) %>%
  summarise(`qtd média de produtos no pedido` = mean(qtd_produtos, na.rm = T)) %>%
  kable() %>%
  kable_styling(bootstrap_options = "striped",
                full_width = F,
                font_size = 12)

#Exploração visual do valor médio do pedido por cidade
fato_pedidos %>%
  group_by(id_cidade) %>%
  mutate(vl_medio_pedido = mean(vl_pedido, na.rm = TRUE)) %>% 
  ggplot() +
  geom_point(aes(x = id_cidade, y = vl_pedido),color = "orange", alpha = 0.5, size = 4) +
  geom_line(aes(x = id_cidade, y = vl_medio_pedido, 
                group = 1, color = "Valor médio do pedido"), size = 1.5) +
  scale_colour_viridis_d() +
  labs(x = "Cidade",
       y = "Valor médio do pedido") +
  theme(legend.title = element_blank(),
        panel.border = element_rect(NA),
        panel.grid = element_line("grey"),
        panel.background = element_rect("white"),
        legend.position = "bottom",
        axis.text.x = element_text(angle = 90))

#Exploração visual do valor médio do pedido pela qtd de produtos no pedido
fato_pedidos %>%
  group_by(qtd_produtos) %>%
  mutate(pedido_medio = mean(vl_pedido, na.rm = TRUE)) %>% 
  ggplot() +
  geom_point(aes(x = qtd_produtos, y = vl_pedido),color = "orange", alpha = 0.5, size = 4) +
  geom_line(aes(x = qtd_produtos, y = pedido_medio, 
                group = 1, color = "Valor médio do pedido"), size = 1.5) +
  scale_colour_viridis_d() +
  labs(x = "Quantidade de produtos",
       y = "Valor médio do pedido") +
  theme(legend.title = element_blank(),
        panel.border = element_rect(NA),
        panel.grid = element_line("grey"),
        panel.background = element_rect("white"),
        legend.position = "bottom",
        axis.text.x = element_text(angle = 90))


#Boxplot da variável dependente (pedido médio)
ggplotly(
  ggplot(fato_pedidos, aes(x = "", y = vl_pedido)) +
    geom_boxplot(fill = "deepskyblue",    # cor da caixa
                 alpha = 0.7,             # transparência
                 color = "black",         # cor da borda
                 outlier.colour = "red",  # cor dos outliers
                 outlier.shape = 15,      # formato dos marcadores dos outliers
                 outlier.size = 2.5) +    # tamanho dos marcadores dos outliers
    geom_jitter(width = 0.1, alpha = 0.5, size = 1.3, color = "darkorchid") +
    labs(y = "Valor do Pedido") +
    theme(panel.background = element_rect("white"),
          panel.grid = element_line("grey95"),
          panel.border = element_rect(NA),
          legend.position="none",
          plot.title = element_text(size=15)) +
    ggtitle("Boxplot da variável 'Valor do Pedido'") +
    xlab("")
)

#Kernel density estimation (KDE) - função densidade de probabilidade da
#variável dependente (valor do pedido), com histograma
ggplotly(
  ggplot(fato_pedidos, aes(x = vl_pedido)) +
    geom_density(aes(x = vl_pedido), 
                 position = "identity", color = "black", size = 1) +
    geom_histogram(aes(y = ..density..), color = "white", fill = "deepskyblue",
                   bins = 30) +
    theme_classic()
)

#Kernel density estimation (KDE) - função densidade de probabilidade da
#variável dependente (qtd de produtos), com histograma
ggplotly(
  ggplot(fato_pedidos, aes(x = qtd_produtos)) +
    geom_density(aes(x = qtd_produtos), 
                 position = "identity", color = "black", size = 1) +
    geom_histogram(aes(y = ..density..), color = "white", fill = "deepskyblue",
                   bins = 30) +
    theme_classic()
)

#Boxplot da variável dependente (Valor do Pedido) por cidade
ggplotly(
  ggplot(fato_pedidos, aes(x = id_cidade,y = vl_pedido)) +
    geom_boxplot(aes(fill = id_cidade)) +
    geom_jitter(width = 0.1, alpha = 0.5, size = 1.3, color = "darkorchid") +
    scale_fill_viridis_d() +
    labs(y = "vl_pedido") +
    theme_classic() +
    ggtitle("Boxplots da variável 'Valor do Pedido' por cidade")
)

#Boxplot da variável dependente (Valor do Pedido) por qtd de produtos no pedido
ggplotly(
  ggplot(fato_pedidos, aes(x = qtd_produtos ,y = vl_pedido)) +
    geom_boxplot(aes(fill = id_cidade, alpha = 0.7)) +
    geom_jitter(width = 0.1, alpha = 0.3, size = 0.3, color = "darkorchid") +
    scale_fill_viridis_d() +
    labs(y = "vl_pedido") +
    theme_classic() +
    ggtitle("Boxplots da variável 'Valor do Pedido' por Qtd de produtos no pedido")
)

#Distribuições da variável 'Valor Pedido' para as cidades, com boxplots
#pacote 'ggdist'
fato_pedidos %>%
  ggplot(aes(x = id_cidade, y = vl_pedido, fill = id_cidade)) +
  ggdist::stat_halfeye(
    adjust = 0.5,
    justification = -.2,
    .width = 0,
    point_colour = NA
  ) +
  geom_boxplot(
    width = .12,
    outlier.color = NA,
    alpha = 0.5
  ) +
  ggdist::stat_dots(
    side = "left",
    justification = 1.1,
    binwidth = .25
  ) +
  scale_fill_viridis_d() +
  theme_tq() +
  labs(
    title = "Distribuições da variável 'Valor Pedidos'",
    subtitle = "com boxplots",
    x = "Cidade",
    y = "Valor do Pedido") +
  coord_flip()

##############################################################
# não usei no trabalho 
#Gráfico de vl_pedido x Salário Médio (OLS)
ggplotly(
  fato_pedidos %>%
    ggplot(aes(x = sal_med, y = vl_pedido)) +
    geom_smooth(method = "lm", formula = y ~ x, se = F) +
    geom_point() +
    scale_colour_viridis_d() +
    labs(x = "Salário Médio Mensal da cidade onde o pedido foi realizado",
         y = "Valores dos pedidos") +
    theme_bw()
)

#Gráfico de Total Pedido x Salário Médio (OLS)
ggplotly(
  fato_pedidos %>%
    ggplot(aes(x = idhm, y = vl_pedido)) +
    geom_smooth(method = "lm", formula = y ~ x, se = F) +
    geom_point() +
    scale_colour_viridis_d() +
    labs(x = "IDHM da cidade onde o pedido foi realizado",
         y = "Valores dos pedidos") +
    theme_bw()
)

##############################################################
#Gráfico alternativo com distribuições da variável 'Total Pedido' para as cidades
#função 'geom_density_ridges_gradient' do pacote 'ggridges'
ggplot(fato_pedidos, aes(x = vl_pedido, y = id_cidade, fill = ..x..)) +
  geom_density_ridges_gradient(scale = 3, rel_min_height = 0.01) +
  scale_fill_viridis(name = "vl_pedido", option = "turbo", direction = -1) +
  labs(
    title = "Distribuições da variável 'Valor Pedido' para as cidades",
    x = "Valor Pedido",
    y = "Cidade") +
  theme_ipsum() +
  theme(
    legend.position="none",
    panel.spacing = unit(0.1, "lines"),
    strip.text.x = element_text(size = 5)
  )

#Gráfico alternativo com distribuições da variável 'Total Pedido' para as cidades
#função 'geom_density_ridges_gradient' do pacote 'ggridges'
ggplot(fato_pedidos, aes(x = qtd_produtos, y = id_cidade, fill = ..x..)) +
  geom_density_ridges_gradient(scale = 3, rel_min_height = 0.01) +
  scale_fill_viridis(name = "Qtd Produtos", option = "turbo", direction = -1) +
  labs(
    title = "Distribuições da variável 'Qtd de Produdos' por cidade",
    x = "Qtd Produtos nos pedidos",
    y = "Cidade") +
  theme_ipsum() +
  theme(
    legend.position="none",
    panel.spacing = unit(0.1, "lines"),
    strip.text.x = element_text(size = 5)
  )

#Kernel density estimation (KDE) - função densidade de probabilidade da
#variável dependente (Vl Pedido) por cidade
ggplotly(
  ggplot(fato_pedidos, aes(x = vl_pedido)) +
    geom_density(aes(x = vl_pedido), 
                 position = "identity", color = "black" , size = 1) +
    geom_histogram( aes(y= ..density.. ), color = "white", fill = "deepskyblue", 
                    bins = 30) +
    theme_classic()
)

#Kernel density estimation (KDE) - função densidade de probabilidade da
#variável dependente (Vl Pedido) por cidade
ggplotly(
  ggplot(fato_pedidos, aes(x = vl_pedido)) +
    geom_density(aes(color = id_cidade, fill = id_cidade), 
                 position = "identity", alpha = 0.3) +
    scale_color_viridis_d() +
    scale_fill_viridis_d() +
    theme_classic()
)

#Kernel density estimation (KDE) - função densidade de probabilidade da
#variável dependente (Vl Pedido) por cidade
ggplotly(
  ggplot(fato_pedidos, aes(x = qtd_produtos)) +
    geom_density(aes(color = id_cidade, fill = id_cidade), 
                 position = "identity", alpha = 0.3) +
    scale_color_viridis_d() +
    scale_fill_viridis_d() +
    theme_classic()
)

#Kernel density estimation (KDE) - função densidade de probabilidade da
#variável dependente (Total Pedido), com histograma e por cidade separadamente
#(função facet_wrap)
fato_pedidos %>% 
  group_by(id_cidade) %>% 
  mutate(linhas = 1:n()) %>% 
  mutate(x = unlist(density(vl_pedido, n = max(linhas))["x"]),
         y = unlist(density(vl_pedido, n = max(linhas))["y"])) %>%
  ggplot() +
  geom_area(aes(x = x, y = y, group = id_cidade, fill = id_cidade), color = "black", alpha = 0.3) +
  geom_histogram(aes(x = vl_pedido, y = ..density.., fill = id_cidade), 
                 color = "black", position = 'identity', alpha = 0.1) +
  facet_wrap(~ id_cidade) +
  scale_fill_viridis_d() +
  scale_color_viridis_d() +
  theme_bw()


#Kernel density estimation (KDE) - função densidade de probabilidade da
#variável dependente (Qtd de Produtos), com histograma e por cidade separadamente
#(função facet_wrap)
fato_pedidos %>% 
  group_by(id_cidade) %>% 
  mutate(linhas = 1:n()) %>% 
  mutate(x = unlist(density(qtd_produtos, n = max(linhas))["x"]),
         y = unlist(density(qtd_produtos, n = max(linhas))["y"])) %>%
  ggplot() +
  geom_area(aes(x = x, y = y, group = id_cidade, fill = id_cidade), color = "black", alpha = 0.3) +
  geom_histogram(aes(x = qtd_produtos, y = ..density.., fill = id_cidade), 
                 color = "black", position = 'identity', alpha = 0.1) +
  facet_wrap(~ id_cidade) +
  scale_fill_viridis_d() +
  scale_color_viridis_d() +
  theme_bw()


#Gráfico de Total Pedido x Salário Médio (OLS)
#ggplotly(
#  fato_pedidos %>%
#    ggplot(aes(x = sal_med, y = vl_pedido)) +
#    geom_smooth(method = "lm", formula = y ~ x, se = F) +
#    geom_point() +
#    scale_colour_viridis_d() +
#    labs(x = "Salário Médio Mensal da cidade onde o pedido foi realizado",
#         y = "Valores dos pedidos") +
#    theme_bw()
#)

#Gráfico de Qtd de produtos x Valor do pedido (OLS)
ggplotly(
  fato_pedidos %>%
    ggplot(aes(x = qtd_produtos, y = vl_pedido)) +
    geom_smooth(method = "lm", formula = y ~ x, se = F) +
    geom_point() +
    scale_colour_viridis_d() +
    labs(x = "Quantidade de produtos x Valor do pedido",
         y = "Valores dos pedidos") +
    theme_bw()
)


#Gráfico de Total Pedido x sal_med por cidade (visualização do contexto)
#NOTE QUE A PERSPECTIVA MULTINÍVEL NATURALMENTE CONSIDERA O COMPORTAMENTO
#HETEROCEDÁSTICO NOS DADOS!
ggplotly(
  fato_pedidos %>%
    ggplot(aes(x = qtd_produtos, y = vl_pedido, color = id_cidade)) +
    geom_smooth(method = "lm", formula = y ~ x, se = F) +
    geom_point() +
    guides(color = "none") +
    scale_colour_viridis_d() +
    labs(x = "Quantidade de produtos x Valor do pedido x Cidade",
         y = "Valores dos Pedidos") +
    theme_bw()
)

#O gráfico a seguir apresenta uma plotagem sob a perspectiva de um modelo
#com equação única (ex.: OLS)
base_exemplo <- fato_pedidos %>%
  filter(id_cidade %in% c("Cidade 1","Cidade 2","Cidade 3","Cidade 4")) %>%
  mutate(id_cidade = as.numeric(id_cidade))

scatter3d(vl_pedido ~ qtd_produtos + id_cidade, #função scatter3d do pacote car
          data = base_exemplo,
          fit = "linear")

################################################################################
################################################################################
#                         ESTIMAÇÃO DO MODELO NULO HLM2                        #
################################################################################

#Estimação do modelo nulo (função lme do pacote nlme)
modelo_nulo_hlm2 <- lme(fixed = vl_pedido ~ 1, 
                        random = ~ 1 | qtd_produtos,
                        data = fato_pedidos,
                        method = "REML") #restricted estimation of maximum likelihood (Gelman)

#Parâmetros do modelo
summary(modelo_nulo_hlm2)

#Verificando a funcionalidade da função 'stderr_nlme' desenvolvida
stderr_nlme(modelo_nulo_hlm2)

################################################################################
#                    COMPARAÇÃO DO HLM2 NULO COM UM OLS NULO                   #
################################################################################
#Para estimarmos o modelo OLS nulo, podemos comandar o seguinte
modelo_ols_nulo <- lm(formula = vl_pedido ~ 1, 
                      data = fato_pedidos)

#Parâmetros do modelo OLS nulo
summary(modelo_ols_nulo)

#Para comparar os LLs dos modelos, vamos utilizar a função lrtest do pacote lmtest
lrtest(modelo_ols_nulo, modelo_nulo_hlm2)

#Comparação entre os LLs dos modelos
data.frame(OLS_Nulo = logLik(modelo_ols_nulo),
           HLM2_Nulo = logLik(modelo_nulo_hlm2)) %>%
  rename(`OLS Nulo` = 1,
         `HLM2 Nulo` = 2) %>%
  melt() %>%
  ggplot(aes(x = variable, y = (abs(-value)), fill = factor(variable))) +
  geom_bar(stat = "identity") +
  geom_label(aes(label = (round(value,3))), hjust = 1.1, color = "white", size = 7) +
  labs(title = "Comparação do LL", 
       y = "LogLik", 
       x = "Modelo Proposto") +
  coord_flip() +
  scale_fill_manual("Legenda:",
                    values = c("grey25","grey45")) +
  theme(legend.title = element_blank(), 
        panel.background = element_rect("white"),
        legend.position = "none",
        axis.line = element_line())

################################################################################
#            ESTIMAÇÃO DO MODELO COM INTERCEPTOS ALEATÓRIOS HLM2               #
################################################################################

#Estimação do modelo com Interceptos Aleatórios
modelo_intercept_hlm2 <- lme(fixed = vl_pedido ~ qtd_produtos,
                             random = ~ 1 | id_cidade,
                             data = fato_pedidos,
                             method = "REML")

#Parâmetros do modelo
summary(modelo_intercept_hlm2)

#Erros-padrão por meio da função 'stderr_nlme' desenvolvida
stderr_nlme(modelo_intercept_hlm2)

#Comparação entre os LLs dos modelos
data.frame(OLS_Nulo = logLik(modelo_ols_nulo),
           HLM2_Nulo = logLik(modelo_nulo_hlm2),
           HLM2_Intercept_Aleat = logLik(modelo_intercept_hlm2)) %>%
  rename(`OLS Nulo` = 1,
         `HLM2 Nulo` = 2,
         `HLM2 com Interceptos Aleatórios` = 3) %>%
  melt() %>%
  ggplot(aes(x = variable, y = (abs(-value)), fill = factor(variable))) +
  geom_bar(stat = "identity") +
  geom_label(aes(label = (round(value,4))), hjust = 1.1, color = "white", size = 7) +
  labs(title = "Comparação do LL", 
       y = "LogLik", 
       x = "Modelo Proposto") +
  coord_flip() +
  scale_fill_manual("Legenda:",
                    values = c("grey25","grey45","bisque4")) +
  theme(legend.title = element_blank(), 
        panel.background = element_rect("white"),
        legend.position = "none",
        axis.line = element_line())

################################################################################
#      ESTIMAÇÃO DO MODELO COM INTERCEPTOS E INCLINAÇÕES ALEATÓRIOS HLM2       #
################################################################################

#Estimação do modelo com Interceptos e Inclinações Aleatórios
modelo_intercept_inclin_hlm2 <- lme(fixed = vl_pedido ~ qtd_produtos,
                                    random = ~ qtd_produtos | id_cidade,
                                    data = fato_pedidos,
                                    method = "REML")

#Parâmetros do modelo
summary(modelo_intercept_inclin_hlm2)

#Erros-padrão por meio da função 'stderr_nlme' desenvolvida
stderr_nlme(modelo_intercept_inclin_hlm2)

#Comparação entre os LLs do modelos
data.frame(OLS_Nulo = logLik(modelo_ols_nulo),
           HLM2_Nulo = logLik(modelo_nulo_hlm2),
           HLM2_Intercept_Aleat = logLik(modelo_intercept_hlm2),
           HLM2_Intercept_Inclin_Aleat = logLik(modelo_intercept_inclin_hlm2)) %>%
  rename(`OLS Nulo` = 1,
         `HLM2 Nulo` = 2,
         `HLM2 com Interceptos Aleatórios` = 3,
         `HLM2 com Interceptos e Inclinações Aleatórios` = 4) %>%
  melt() %>%
  ggplot(aes(x = variable, y = (abs(-value)), fill = factor(variable))) +
  geom_bar(stat = "identity") +
  geom_label(aes(label = (round(value,4))), hjust = 1.1, color = "white", size = 6) +
  labs(title = "Comparação do LL", 
       y = "LogLik", 
       x = "Modelo Proposto") +
  coord_flip() +
  scale_fill_manual("Legenda:",
                    values = c("grey25","grey45","bisque4","bisque3")) +
  theme(legend.title = element_blank(), 
        panel.background = element_rect("white"),
        legend.position = "none",
        axis.line = element_line())

################################################################################
#                       ESTIMAÇÃO DO MODELO FINAL HLM2                         #
################################################################################

#Estimação do modelo final
modelo_final_hlm2 <- lme(fixed = vl_pedido ~ qtd_produtos + idhm + qtd_produtos:vl_pedido,
                         random = ~ qtd_produtos | id_cidade,
                         data = fato_pedidos,
                         method = "REML")

#Parâmetros do modelo
summary(modelo_final_hlm2)

#Erros-padrão por meio da função 'stderr_nlme' desenvolvida
stderr_nlme(modelo_final_hlm2)

#Comparação entre os LLs do modelos
data.frame(OLS_Nulo = logLik(modelo_ols_nulo),
           HLM2_Nulo = logLik(modelo_nulo_hlm2),
           HLM2_Intercept_Aleat = logLik(modelo_intercept_hlm2),
           HLM2_Intercept_Inclin_Aleat = logLik(modelo_intercept_inclin_hlm2),
           HLM2_Modelo_Final = logLik(modelo_final_hlm2)) %>%
  rename(`OLS Nulo` = 1,
         `HLM2 Nulo` = 2,
         `HLM2 com Interceptos Aleatórios` = 3,
         `HLM2 com Interceptos e Inclinações Aleatórios` = 4,
         `HLM2 Modelo Final` = 5) %>%
  melt() %>%
  ggplot(aes(x = variable, y = (abs(-value)), fill = factor(variable))) +
  geom_bar(stat = "identity") +
  geom_label(aes(label = (round(value,4))), hjust = 1.1, color = "white", size = 6) +
  labs(title = "Comparação do LL", 
       y = "LogLik", 
       x = "Modelo Proposto") +
  coord_flip() +
  scale_fill_manual("Legenda:",
                    values = c("grey25","grey45","bisque4","bisque3",
                               "deepskyblue1")) +
  theme(legend.title = element_blank(), 
        panel.background = element_rect("white"),
        legend.position = "none",
        axis.line = element_line())

#Melhor visualização dos interceptos e das inclinações aleatórios por cidade,
#para o modelo final HLM2

v_final <- data.frame(modelo_final_hlm2[["coefficients"]][["random"]][["id_cidade"]]) %>%
  rename(v00 = 1,
         v10 = 2)
v_final$id_cidade <- c(1:4)
v_final$id_cidade <- as.factor(v_final$id_cidade)

v_final %>% 
  select(id_cidade, everything()) %>%
  kable() %>%
  kable_styling(bootstrap_options = "striped", 
                full_width = F, 
                font_size = 25)

#Para observarmos graficamente o comportamento dos valores de v0j, ou seja,
#dos interceptos aleatórios por cidade, podemos comandar
###################### corrigir .... não sei o que
random.effects(modelo_final_hlm2) %>% 
  rename(v0j = 1) %>% 
  rownames_to_column("id_cidade") %>% 
  mutate(color_v0j = ifelse(v0j < 0, "A", "B"),
         hjust_v0j = ifelse(v0j > 0, 1.15, -0.15)) %>% 
  arrange(id_cidade) %>% 
  ggplot(aes(label = format(v0j, digits = 2), 
             hjust = hjust_v0j)) +
  geom_bar(aes(x = fct_rev(id_cidade), y = v0j, fill = color_v0j),
           stat = "identity", color = "black") +
  geom_text(aes(x = id_cidade, y = 0), size = 4.1, color = "black") +
  coord_flip() +
  labs(x = "Cidade",
       y = expression(nu[0][j])) +
  scale_fill_manual("oi", values = c("firebrick1","green1")) +
  theme(panel.background = element_rect("white"),
        panel.border = element_rect(NA),
        panel.grid = element_line("grey95"),
        legend.position = "none")

#Para observarmos graficamente o comportamento dos valores de v1j, ou seja
#das inclinações aleatórias por cidade, podemos comandar

random.effects(modelo_final_hlm2) %>% 
  rename(v1j = 2) %>% 
  rownames_to_column("id_cidade") %>% 
  mutate(color_v1j = ifelse(v1j < 0, "A", "B"),
         hjust_v1j = ifelse(v1j > 0, 1.15, -0.15)) %>% 
  arrange(id_cidade) %>% 
  ggplot(aes(label = format(v1j, digits = 2), 
             hjust = hjust_v1j)) +
  geom_bar(aes(x = fct_rev(id_cidade), y = v1j, fill = color_v1j),
           stat = "identity", color = "black") +
  geom_text(aes(x = id_cidade, y = 0), size = 4.1, color = "black") +
  coord_flip() +
  labs(x = "Cidade",
       y = expression(nu[1][j])) +
  scale_fill_manual("oi", values = c("firebrick1","green1")) +
  theme(panel.background = element_rect("white"),
        panel.border = element_rect(NA),
        panel.grid = element_line("grey95"),
        legend.position = "none")

#Gerando os fitted values do modelo HLM2 Final
fato_pedidos$hlm2_fitted <- predict(modelo_final_hlm2,
                                        fato_pedidos)

# Visualizando os fitted values do modelo
#Visualizando os fitted values por estudante e por escola
############### corrigir a conversão id_cidade em número
predict(modelo_final_hlm2, level = 0:1) %>% 
  mutate(id_cidade = gsub("^.*?\\/","",id_cidade),
         id_cidade = as.factor(as.numeric(id_cidade)),
         Total_Pedido = fato_pedidos$Total_Pedido,
         etjk = resid(modelo_final_hlm2)) %>% #função resid gera os termos etjk
  select(id_cidade, Total_Pedido, everything()) %>%
  kable() %>%
  kable_styling(bootstrap_options = "striped", 
                full_width = F, 
                font_size = 18)


################################################################################
################################## parei aqui ##################################
################################################################################
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

