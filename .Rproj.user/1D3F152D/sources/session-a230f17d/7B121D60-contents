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
#                      Cronstrução e EXPLORAÇÃO DO DATASET                       #
################################################################################

# Carregamento da base de dados
# Montar os dataframes com as observações abaixo

df_cidade <- read.csv("cidades.csv", header = TRUE)

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

###################################################################
pedidos <- read.csv("TCC_Fato_Pedido_Novo.csv", header = TRUE)
dados_f_pedidos <- pedidos[,c("id_pedido", "id_cliente", "id_cidade", "id_produto", "qtd_tot_gr", "Preço_Venda", "Total_Pedido")]
dados_f_pedidos$id_produto <- as.factor(dados_f_pedidos$id_produto)
dados_f_pedidos$id_cidade <- as.factor(dados_f_pedidos$id_cidade)
dados_f_pedidos$id_cliente <- as.factor(dados_f_pedidos$id_cliente)
colnames(dados_f_pedidos)[6] <- "vl_produto"
colnames(dados_f_pedidos)[7] <- "vl_pedido"

dados_f_pedidos <- subset(dados_f_pedidos, qtd_tot_gr != 0)

dados_f_pedidos$vl_produto <- gsub("[^0-9\\.\\,]", "", dados_f_pedidos$vl_produto) # remove caracteres não numéricos
dados_f_pedidos$vl_produto <- gsub(",", ".", dados_f_pedidos$vl_produto) # substitui vírgulas por pontos
dados_f_pedidos$vl_produto <- as.numeric(dados_f_pedidos$vl_produto) # converte para numérico
dados_f_pedidos$vl_pedido <- gsub("[^0-9\\.\\,]", "", dados_f_pedidos$vl_pedido) # remove caracteres não numéricos
dados_f_pedidos$vl_pedido <- gsub(",", ".", dados_f_pedidos$vl_pedido) # substitui vírgulas por pontos
dados_f_pedidos$vl_pedido <- as.numeric(dados_f_pedidos$vl_pedido) # converte para numérico

dados_f_pedidos <- dados_f_pedidos %>%
  group_by(id_pedido) %>%
  mutate(qtd_produtos = n(), qtd_tot_gr = sum(qtd_tot_gr) )

###################################################################

#pedidos <- read.csv("TCC_Fato_Pedido.csv", header = TRUE)
#dados_f_pedidos <- pedidos[,c("id_pedido", "id_cliente", "id_cidade", "id_produto", "Dt_Pedido", "Preço_Venda", "Total_Pedido")]
#dados_f_pedidos$id_produto <- as.factor(dados_f_pedidos$id_produto)
#dados_f_pedidos$id_cidade <- as.factor(dados_f_pedidos$id_cidade)
#dados_f_pedidos$id_cliente <- as.factor(dados_f_pedidos$id_cliente)
#colnames(dados_f_pedidos)[6] <- "vl_produto"
#colnames(dados_f_pedidos)[7] <- "vl_pedido"

#dados_f_pedidos$vl_produto <- gsub("[^0-9\\.\\,]", "", dados_f_pedidos$vl_produto) # remove caracteres não numéricos
#dados_f_pedidos$vl_produto <- gsub(",", ".", dados_f_pedidos$vl_produto) # substitui vírgulas por pontos
#dados_f_pedidos$vl_produto <- as.numeric(dados_f_pedidos$vl_produto) # converte para numérico
#dados_f_pedidos$vl_pedido <- gsub("[^0-9\\.\\,]", "", dados_f_pedidos$vl_pedido) # remove caracteres não numéricos
#dados_f_pedidos$vl_pedido <- gsub(",", ".", dados_f_pedidos$vl_pedido) # substitui vírgulas por pontos
#dados_f_pedidos$vl_pedido <- as.numeric(dados_f_pedidos$vl_pedido) # converte para numérico

#dados_f_pedidos <- dados_f_pedidos %>%
#  group_by(id_pedido) %>%
#  mutate(qtd_produtos = n() )

 
#visão PRODUTOS
#=======================
#dados_f_produtos <- pedidos %>%
#  group_by(id_produto,id_sexo) %>%
#  mutate(qtd_clientes = n(), qtd_total_prod = sum(QTD))

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

#juntando os dados visão PRODUTO
#========================================
#fato_produtos <- merge(dados_f_produtos, df_cidade, by = "id_cidade")
#summary(fato_produtos)
#fato_produtos <- merge(fato_produtos, clientes, by = "id_cliente")
#summary(fato_produtos)
#fato_produtos$id_produtos <- as.factor(fato_produtos$id_produto)
#summary(fato_produtos)

# data frame final
fato_pedidos <- fato_pedidos[,c("id_pedido", "id_cliente", "id_sexo", "id_cidade", "vl_pedido", "qtd_produtos", "qtd_tot_gr", "habitantes", "pib", "idhm", "sal_med")]
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

#############################################################################
#############################################################################
#############################################################################

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
#fato_pedidos %>%
#  group_by(qtd_produtos) %>%
#  mutate(pedido_medio = mean(vl_pedido, na.rm = TRUE)) %>% 
#  ggplot() +
#  geom_point(aes(x = qtd_produtos, y = vl_pedido),color = "orange", alpha = 0.5, size = 4) +
#  geom_line(aes(x = qtd_produtos, y = pedido_medio, 
#                group = 1, color = "Valor médio do pedido"), size = 1.5) +
#  scale_colour_viridis_d() +
#  labs(x = "Quantidade de produtos",
#       y = "Valor médio do pedido") +
#  theme(legend.title = element_blank(),
#        panel.border = element_rect(NA),
#        panel.grid = element_line("grey"),
#        panel.background = element_rect("white"),
#        legend.position = "bottom",
#        axis.text.x = element_text(angle = 90))


#Boxplot do valor do pedido
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
#do valor do pedido, com histograma
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
#Kernel density estimation (KDE) - função densidade de probabilidade da
#variável dependente (qtd de produtos), com histograma
ggplotly(
  ggplot(fato_pedidos, aes(x = qtd_tot_gr)) +
    geom_density(aes(x = qtd_tot_gr), 
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
#Gráfico de vl_pedido x qtd_tot_gr (OLS)
ggplotly(
  fato_pedidos %>%
    ggplot(aes(x = vl_pedido, y = qtd_tot_gr)) +
    geom_smooth(method = "lm", formula = y ~ x, se = F) +
    geom_point(size = 0.15) +
    scale_colour_viridis_d() +
    labs(x = "Valor do Pedido",
         y = "Quantidade em grama(s)") +
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
##############################################################
#Gráfico alternativo com distribuições da variável 'QTD TOTAL GRAMAS' para as cidades
#função 'geom_density_ridges_gradient' do pacote 'ggridges'
ggplot(fato_pedidos, aes(x = qtd_tot_gr, y = id_cidade, fill = ..x..)) +
  geom_density_ridges_gradient(scale = 3, rel_min_height = 0.01) +
  scale_fill_viridis(name = "qtd_tot_gr", option = "turbo", direction = -1) +
  labs(
    title = "Distribuições da variável 'qtd_tot_gr' para as cidades",
    x = "qtd_tot_gr",
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
#variável dependente (qtd_tot_gr) por cidade
ggplotly(
  ggplot(fato_pedidos, aes(x = qtd_tot_gr)) +
    geom_density(aes(x = qtd_tot_gr), 
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
#variável qtd_tot_gr por cidade
ggplotly(
  ggplot(fato_pedidos, aes(x = qtd_tot_gr)) +
    geom_density(aes(color = id_cidade, fill = id_cidade), 
                 position = "identity", alpha = 0.3) +
    scale_color_viridis_d() +
    scale_fill_viridis_d() +
    theme_classic()
)

#Kernel density estimation (KDE) - função densidade de probabilidade da
#variável qtd_produtos) por cidade
ggplotly(
  ggplot(fato_pedidos, aes(x = qtd_produtos)) +
    geom_density(aes(color = id_cidade, fill = id_cidade), 
                 position = "identity", alpha = 0.3) +
    scale_color_viridis_d() +
    scale_fill_viridis_d() +
    theme_classic()
)

#Kernel density estimation (KDE) - função densidade de probabilidade da
#variável vl_pedido), com histograma e por cidade separadamente
#(função facet_wrap)
fato_pedidos %>% 
  group_by(id_cidade) %>% 
  mutate(linhas = 1:n()) %>% 
  mutate(x = unlist(density(vl_pedido, n = max(linhas))["x"]),
         y = unlist(density(vl_pedido, n = max(linhas))["y"])) %>%
  ggplot() +
  geom_area(aes(x = x, y = y, group = id_cidade, fill = id_cidade), color = "black", alpha = 0.3) +
  geom_histogram(aes(x = vl_pedido, y = ..density.., fill = id_cidade), 
                 color = "black", position = 'identity', alpha = 0.1, bins = 30) +
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

#Gráfico de Qtd de produtos x Valor do pedido (OLS)
ggplotly(
  fato_pedidos %>%
    ggplot(aes(x = qtd_produtos, y = vl_pedido)) +
    geom_smooth(method = "lm", formula = y ~ x, se = F) +
    geom_point() +
    scale_colour_viridis_d() +
    labs(x = "Quantidade de produtos",
         y = "Valores dos pedidos") +
    theme_bw()
)
#Gráfico de qtd_tot_gr x Valor do pedido (OLS)
ggplotly(
  fato_pedidos %>%
    ggplot(aes(x = vl_pedido, y = qtd_tot_gr)) +
    geom_smooth(method = "lm", formula = y ~ x, se = F) +
    geom_point(size = 0.1) +
    scale_colour_viridis_d() +
    labs(x = "vl_pedido",
         y = "qtd_tot_gr") +
    theme_bw()
)


#Gráfico vl_pedidos x qtd_produtos
ggplotly(
  fato_pedidos %>%
    ggplot(aes(x = vl_pedido , y = qtd_produtos, color = id_cidade)) +
    geom_smooth(method = "lm", formula = y ~ x, se = F) +
    geom_point(size = 0.15) +
    guides(color = "none") +
    scale_colour_viridis_d() +
    labs(x = "vl_pedido",
         y = "qtd_produtos") +
    theme_bw()
)

#Gráfico de qtd_tot_gr x vl_pedido por cidade (visualização do contexto)
#NOTE QUE A PERSPECTIVA MULTINÍVEL NATURALMENTE CONSIDERA O COMPORTAMENTO
#HETEROCEDÁSTICO NOS DADOS!
ggplotly(
  fato_pedidos %>%
    ggplot(aes(x = vl_pedido, y = qtd_tot_gr, color = id_cidade)) +
    geom_smooth(method = "lm", formula = y ~ x, se = F) +
    geom_point(size = 0.15) +  # Ajuste o tamanho dos pontos conforme necessário
    guides(color = "none") +
    scale_colour_viridis_d() +
    labs(x = "Valores dos Pedidos",
         y = "qtd_tot_gr") +
    theme_bw()
)

################################################################################
################################################################################
#                         ESTIMAÇÃO DO MODELO NULO HLM2                        #
################################################################################

#Estimação do modelo nulo (função lme do pacote nlme)
modelo_nulo_hlm2 <- lme(fixed = qtd_tot_gr ~ 1, 
                        random = ~ 1 | id_cidade,
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
modelo_ols_nulo <- lm(formula = qtd_tot_gr ~ 1, 
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

#############################################################
#grava o data frame em arquivo
save(fato_pedidos, file = "TCC_dados.RData")

