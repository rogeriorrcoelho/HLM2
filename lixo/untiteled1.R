#Para observarmos graficamente o comportamento dos valores de v0j, ou seja,
#dos interceptos aleatórios por escola, podemos comandar
random.effects(modelo_final_hlm2) %>% 
  rename(v0j = 1) %>% 
  rownames_to_column("Escola") %>% 
  mutate(color_v0j = ifelse(v0j < 0, "A", "B"),
         hjust_v0j = ifelse(v0j > 0, 1.15, -0.15)) %>% 
  arrange(Escola) %>% 
  ggplot(aes(label = format(v0j, digits = 2), 
             hjust = hjust_v0j)) +
  geom_bar(aes(x = fct_rev(Escola), y = v0j, fill = color_v0j),
           stat = "identity", color = "black") +
  geom_text(aes(x = Escola, y = 0), size = 4.1, color = "black") +
  coord_flip() +
  labs(x = "Escola",
       y = expression(nu[0][j])) +
  scale_fill_manual("oi", values = c("firebrick1","green1")) +
  theme(panel.background = element_rect("white"),
        panel.border = element_rect(NA),
        panel.grid = element_line("grey95"),
        legend.position = "none")


#Para observarmos graficamente o comportamento dos valores de v1j, ou seja
#das inclinações aleatórias por escola, podemos comandar
random.effects(modelo_final_hlm2) %>% 
  rename(v1j = 2) %>% 
  rownames_to_column("Escola") %>% 
  mutate(color_v1j = ifelse(v1j < 0, "A", "B"),
         hjust_v1j = ifelse(v1j > 0, 1.15, -0.15)) %>% 
  arrange(Escola) %>% 
  ggplot(aes(label = format(v1j, digits = 2), 
             hjust = hjust_v1j)) +
  geom_bar(aes(x = fct_rev(Escola), y = v1j, fill = color_v1j),
           stat = "identity", color = "black") +
  geom_text(aes(x = Escola, y = 0), size = 4.1, color = "black") +
  coord_flip() +
  labs(x = "Escola",
       y = expression(nu[1][j])) +
  scale_fill_manual("oi", values = c("firebrick1","green1")) +
  theme(panel.background = element_rect("white"),
        panel.border = element_rect(NA),
        panel.grid = element_line("grey95"),
        legend.position = "none")

#Gerando os fitted values do modelo HLM2 Final
estudante_escola$hlm2_fitted <- predict(modelo_final_hlm2,
                                        estudante_escola)

# Visualizando os fitted values do modelo
#Visualizando os fitted values por estudante e por escola
predict(modelo_final_hlm2, level = 0:1) %>% 
  mutate(escola = gsub("^.*?\\/","",escola),
         escola = as.factor(as.numeric(escola)),
         desempenho = estudante_escola$desempenho,
         etjk = resid(modelo_final_hlm2)) %>% #função resid gera os termos etjk
  select(escola, desempenho, everything()) %>%
  kable() %>%
  kable_styling(bootstrap_options = "striped", 
                full_width = F, 
                font_size = 18)

#Efetuando predições
#Exemplo: Quais os valores previstos de desempenho escolar, para dado
#aluno que estuda na escola "1", sabendo-se que ele estuda 11h semanais,
#e que a escola oferece tempo médio de experiência de seus professores
#igual a 3.6 anos?
predict(modelo_final_hlm2, level = 0:1,
        newdata = data.frame(escola = "1",
                             horas = 11,
                             texp = 3.6))

#Valores previstos do desempenho escolar em função da variável horas para o 
#modelo final HLM2 com interceptos e inclinações aleatórios
estudante_escola %>%
  mutate(fitted_escola = predict(modelo_final_hlm2, level = 1)) %>% 
  ggplot() +
  geom_point(aes(x = horas, y = fitted_escola)) +
  geom_smooth(aes(x = horas, y = fitted_escola, color = factor(escola)), 
              method = "lm", se = F) +
  scale_colour_viridis_d() +
  labs(x = "Quantidade Semanal de Horas de Estudo do Aluno",
       y = "Desempenho Escolar (Fitted Values)") +
  theme_bw()


################################################################################
#                       COMPARAÇÃO COM UM MODELO OLS                           #
################################################################################

#Elaborando um modelo OLS para fins de comparação
modelo_ols <- lm(formula = desempenho ~ horas + texp,
                 data = estudante_escola)

#Parâmetros
summary(modelo_ols)

#Comparando os LL dos modelos elaborados
data.frame(OLS = logLik(modelo_ols),
           HLM2_Modelo_Final = logLik(modelo_final_hlm2)) %>%
  rename(`OLS` = 1,
         `HLM2 Modelo Final` = 2) %>%
  melt() %>%
  ggplot(aes(x = variable, y = (abs(-value)), fill = factor(variable))) +
  geom_bar(stat = "identity") +
  geom_label(aes(label = (round(value,3))), hjust = 1.1, color = "white", size = 7) +
  labs(title = "Comparação do LL", 
       y = "LogLik", 
       x = "Modelo Proposto") +
  coord_flip() +
  scale_fill_manual("Legenda:",
                    values = c("darkorchid","deepskyblue1")) +
  theme(legend.title = element_blank(), 
        panel.background = element_rect("white"),
        legend.position = "none",
        axis.line = element_line())

#LR Test
lrtest(modelo_ols, modelo_final_hlm2)

#Comparando a aderência dos fitted values dos modelos estimados
#Gerando os fitted values do modelo OLS
estudante_escola$ols_fitted <- modelo_ols$fitted.values

#Plotagem
estudante_escola %>%
  ggplot() +
  geom_smooth(aes(x = desempenho, y = ols_fitted, color = "OLS"),
              method = "lm", se = F, formula = y ~ splines::bs(x, df = 5),
              size = 1.5) +
  geom_smooth(aes(x = desempenho, y= hlm2_fitted, color = "HLM2 Final"),
              method = "lm", se = F, formula = y ~ splines::bs(x, df = 5),
              size = 1.5) +
  geom_smooth(aes(x = desempenho, y = desempenho), method = "lm", 
              color = "gray44", size = 1.05,
              linetype = "longdash") +
  geom_point(aes(x = desempenho, y = ols_fitted,
                 color = "OLS")) +
  geom_point(aes(x = desempenho, y = hlm2_fitted,
                 color = "HLM2 Final"))  +
  scale_color_manual("Modelos:", 
                     values = c("deepskyblue1","darkorchid")) +
  labs(x = "Desempenho", y = "Fitted Values") +
  theme_bw()


################################################################################
#                 COMPARAÇÃO COM UM MODELO OLS COM DUMMIES                     #
################################################################################

#Procedimento n-1 dummies para o contexto
estudante_escola_dummies <- dummy_cols(.data = estudante_escola,
                                       select_columns = "escola",
                                       remove_first_dummy = TRUE,
                                       remove_selected_columns = TRUE)

#Visualizando as dummies na nova base de dados 'estudante_escola_dummies'
estudante_escola_dummies %>%
  select(-hlm2_fitted,-ols_fitted, everything()) %>%
  kable() %>%
  kable_styling(bootstrap_options = "striped",
                full_width = F,
                font_size = 19)

#Modelo OLS com dummies
modelo_ols_dummies <- lm(formula = desempenho ~ horas + texp + escola_2 +
                           escola_3 + escola_4 + escola_5 + escola_6 +
                           escola_7 + escola_8 + escola_9 + escola_10,
                         data = estudante_escola_dummies)

#Parâmetros
summary(modelo_ols_dummies)

#Procedimento stepwise
modelo_ols_dummies_step <- step(object = modelo_ols_dummies,
                                step = qchisq(p = 0.05, df = 1,
                                              lower.tail = FALSE))

#Parâmetros do modelo OLS estimado com dummies por escola a partir do
#procedimento Stepwise
summary(modelo_ols_dummies_step)

#Comparando os LL dos modelos HLM2 Final, OLs e OLS com Dummies e Stepwise
data.frame(OLS = logLik(modelo_ols),
           OLS_Dummies_Step = logLik(modelo_ols_dummies_step),
           HLM2_Modelo_Final = logLik(modelo_final_hlm2)) %>%
  rename(`OLS` = 1,
         `OLS com Dummies e Stepwise` = 2,
         `HLM2 Modelo Final` = 3) %>%
  melt() %>%
  ggplot(aes(x = variable, y = (abs(-value)), fill = factor(variable))) +
  geom_bar(stat = "identity") +
  geom_label(aes(label = (round(value,3))), hjust = 1.1, color = "white", size = 7) +
  labs(title = "Comparação do LL", 
       y = "LogLik", 
       x = "Modelo Proposto") +
  coord_flip() +
  scale_fill_manual("Legenda:",
                    values = c("darkorchid","maroon1","deepskyblue1")) +
  theme(legend.title = element_blank(), 
        panel.background = element_rect("white"),
        legend.position = "none",
        axis.line = element_line())

#LR Test
lrtest(modelo_ols_dummies_step, modelo_final_hlm2)

#Comparação entre os parãmetros dos modelos (atente-se para a quantidade de
#parâmetros estimados em cada um deles!)
export_summs(modelo_ols_dummies_step, modelo_final_hlm2,
             model.names = c("OLS com Dummies", "HLM2 Final"))


#Comparando a aderência dos fitted values dos modelos HLM2 Final, OLS e
#OLS com Dummies e Stepwise
#Gerando os fitted values do modelo OLS com Dummies e Stepwise
estudante_escola$ols_step_fitted <- modelo_ols_dummies_step$fitted.values

#Gráfico para a comparação entre os fitted values dos modelos HLM2 Final, OLs e
#OLS com Dummies e Procedimento Stepwise
estudante_escola %>%
  ggplot() +
  geom_smooth(aes(x = desempenho, y = ols_step_fitted, color = "OLS com Dummies"),
              method = "lm", se = F, formula = y ~ splines::bs(x, df = 5),
              size = 1.5) +
  geom_smooth(aes(x = desempenho, y= hlm2_fitted, color = "HLM2 Final"),
              method = "lm", se = F, formula = y ~ splines::bs(x, df = 5),
              size = 1.5) +
  geom_smooth(aes(x = desempenho, y= ols_fitted, color = "OLS"),
              method = "lm", se = F, formula = y ~ splines::bs(x, df = 5),
              size = 1.5) +
  geom_smooth(aes(x = desempenho, y = desempenho), method = "lm", 
              color = "gray44", size = 1.05,
              linetype = "longdash") +
  scale_color_manual("Modelos:", 
                     values = c("deepskyblue1", "maroon1", "darkorchid")) +
  labs(x = "Desempenho", y = "Fitted Values") +
  theme_bw()


#Comparação entre os LLs de todos os modelos estimados neste exemplo
data.frame(OLS_Nulo = logLik(modelo_ols_nulo),
           HLM2_Nulo = logLik(modelo_nulo_hlm2),
           OLS = logLik(modelo_ols),
           HLM2_Intercept_Aleat = logLik(modelo_intercept_hlm2),
           OLS_Dummies_step = logLik(modelo_ols_dummies_step),
           HLM2_Intercept_Inclin_Aleat = logLik(modelo_intercept_inclin_hlm2),
           HLM2_Modelo_Final = logLik(modelo_final_hlm2)) %>%
  rename(`OLS Nulo` = 1,
         `HLM2 Nulo` = 2,
         `OLS` = 3,
         `HLM2 com Interceptos Aleatórios` = 4,
         `OLS com Dummies e Stepwise` = 5,
         `HLM2 com Interceptos e Inclinações Aleatórios` = 6,
         `HLM2 Modelo Final` = 7) %>%
  melt() %>%
  ggplot(aes(x = variable, y = (abs(-value)), fill = factor(variable))) +
  geom_bar(stat = "identity") +
  geom_label(aes(label = (round(value,3))), hjust = 1.1, color = "white", size = 5) +
  labs(title = "Comparação do LL", 
       y = "LogLik", 
       x = "Modelo Proposto") +
  coord_flip() +
  scale_fill_manual("Legenda:",
                    values = c("grey25","grey45","darkorchid","bisque4",
                               "maroon1","bisque3","deepskyblue1")) +
  theme(legend.title = element_blank(), 
        panel.background = element_rect("white"),
        legend.position = "none",
        axis.line = element_line())


##################################### FIM ######################################
