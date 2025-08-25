setwd("C:/Users/rogerio/Desktop/Power BI")
meu_data_frame <- read.csv("NYPD_Complaint_Data_Historic.csv")

df_nypd <- meu_data_frame %>%
  filter(OFNS_DESC == "ASSAULT 3 & RELATED OFFENSES")

df_nypd <- subset(meu_data_frame, select = c("BORO_NM", "OFNS_DESC", "PREM_TYP_DESC", "SUSP_RACE", "SUSP_SEX", "VIC_RACE"))
resultado <- aggregate(BORO_NM ~ BORO_NM + OFNS_DESC , data = df_nypd, FUN = length)

resultado <- df_nypd %>%
  group_by(BORO_NM, OFNS_DESC) %>%
  summarise(ocorrencias = n())
#gráfico boxplot das ocorrencias por bairro/vizinhança
#Boxplot da variável dependente (ocorrencias)
ggplotly(
  ggplot(resultado, aes(x = BORO_NM, y = ocorrencias)) +
    geom_boxplot(fill = "blue",    # cor da caixa
                 alpha = 0.7,             # transparência
                 color = "black",         # cor da borda
                 outlier.colour = "orange",  # cor dos outliers
                 outlier.shape = 15,      # formato dos marcadores dos outliers
                 outlier.size = 2.5) +    # tamanho dos marcadores dos outliers
    geom_jitter(width = 0.1, alpha = 0.5, size = 1.3, color = "red") +
    labs(y = "ocorrencias") +
    theme(panel.background = element_rect("white"),
          panel.grid = element_line("grey95"),
          panel.border = element_rect(NA),
          legend.position="none",
          plot.title = element_text(size=15)) +
    ggtitle("Boxplot da variável 'ocorrências' por Bairro") +
    xlab("")
)
