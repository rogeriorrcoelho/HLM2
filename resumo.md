# Regressão Multinível Tipo 2 (HLM2) - Revenda de Produtos Orgânicos

Análise dos preços dos pedidos por meio da modelagem multinível.

---

## Objetivo

Este projeto visa estudar os relacionamentos, padrões e comportamentos do consumo de produtos orgânicos, utilizando técnicas de regressão multinível para construção de um modelo preditivo de machine learning. Usou-se um dataset de um microempreendedor do segmento de revenda de produtos orgânicos do interior de São Paulo para construir o modelo linear hierárquico do tipo 2 (HLM2)

<p align="center">
  <img src="/imagens/estrutura_agrupamento.jpeg" />
</p>

---

## Análise Exploratória

- Foram analisados o balanceamento dos dados, média dos pedidos por cidade, densidade dos pedidos e características dos municípios onde a microempresa atua.
- A análise multinível permite considerar contextos sociais e atores individuais, com o valor dos pedidos como nível mais baixo e a cidade como primeiro nível.
- A maioria dos clientes ativos é do sexo feminino. O número de pedidos por sexo acompanhou a razão de 9 para 1 em favor das mulheres, sugerindo papel fundamental delas na promoção de hábitos alimentares mais saudáveis.

<p align="center">
  <img src="/imagens/tab_cad_cliente.jpeg" />
</p>

---

## Descobertas

- A cidade com a maior quantidade de pedidos não é a mais populosa, nem a cidade com menos pedidos é a de menor população.
- Não há correspondência clara entre tamanho populacional e número de pedidos, mas há entre valor médio do pedido e quantidade média de produtos.
- A visualização por boxplot revelou outliers, que refletem clientes recorrentes e não devem ser removidos.
- A sazonalidade e datas comemorativas influenciam os valores dos pedidos, produzindo outliers valiosos para análise.

<p align="center">
  <img src="/imagens/tab_desbalanceamento.jpeg" />
</p>

<p align="center">
  <img src="/imagens/vl_medio_pedido.jpeg" />
</p>

---

## Distribuição dos Pedidos

- O gráfico de densidade revelou assimetria e a existência de múltiplos picos nos valores dos pedidos, sugerindo subpopulações de clientes.
- Diversos picos foram identificados (~58, ~143, ~220, ~258), indicando diferentes perfis
- Os valores variam de próximos a zero até 300, revelando grande variabilidade nos pedidos.
- Cidade 1 tem mais pedidos de baixo valor; Cidade 4 se destaca por pedidos de maior valor devido à distância da sede e à política de entregas.

---

## Análise Multinível de Grupos de Pedido

- No histograma da quantidade de produtos por pedido surgiram descontinuidades, sugerindo subgrupos como solteiros, famílias pequenas, médias e grandes.
- O gráfico de densidade por cidade mostra que cidades diferentes possuem diferentes comportamentos de concentração e dispersão dos valores dos pedidos
- Segundo literatura, escolaridade e acesso à informação influenciam o consumo de orgânicos nessas localidades.

---

## Modelagem e Regressores

- A dispersão entre valor do pedido (vl_pedido) e peso total do pedido (qtd_tot_gr) há relação proporcional.
- A regressão linear apresenta boa aproximação para a variável resposta. Cada cidade pode ter um intercepto e inclinação distintos na regressão linear individual
- O uso da modelagem multinível permite interceptos diferenciados por cidade, capturando características contextuais.

---

## Estratégia e Seleção de Variáveis

- No modelo multinível a seleção step-up inicia pelo modelo nulo, adicionando progressivamente variáveis e avaliando impacto por critérios estatísticos (valor-p, AIC, BIC).
- O modelo nulo HLM2 não apresentou resultados estatisticamente significativos para efeito aleatório de intercepto (p-valor = 0,257).
- O ICC (“Intraclass Correlation Coefficient”) apontou que 8,9% da variação do valor do pedido se deve à cidade, sem significância estatística notable.

---

## Conclusões

- O modelo final HLM2 teve ajuste levemente melhor, mas sem ganho estatístico frente ao modelo de regressão linear simples (OLS)
- A ausência de informações complementares sobre clientes e vendas limitou agrupamentos hierárquicos adicionais e a entrada de mais variáveis explicativas.
- A maturidade operacional do microempreendedor em gestão de informações ainda tem espaço para evoluir
- Recomenda-se, neste caso, o uso de regressão linear simples para predição dos valores de pedidos

---

