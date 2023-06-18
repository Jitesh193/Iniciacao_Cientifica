# Iniciação Científica - Estimativa de crescimento e idade de tumores através da técnica de Efeitos Mistos da população e indivíduo

# Resumo

O tema deste projeto é a obtenção de estimativas da forma de crescimento e da idade de tumores utilizando modelos de crescimentos agregados de tumores, como os Modelos Logístico ou o de Gompertz, por exemplo. Utilizaremos uma estrutura de Estimativa Bayesiana conhecida como estimativa de Efeitos Mistos não-linear (nonlinear Mixed-Effect Estimation) com distribuição Fat-tail,  a partir de poucos dados disponíveis para cada indivíduo, porém em número razoável quando agregados de populações semelhantes sujeitas ao mesmo tipo de carcinoma.

Uma distribuição mista Laplaciana-Gaussiana é usada para alcançar maior robustez em vista da dispersão de dados, em geral mal representados pelos modelos convencionais Gaussianos e métodos associados como os de quadrado mínimo e suas variações. Modelos estatísticos mais robustos e flexíveis são críticos neste tipo de estimativas nas quais poucos dados para um indivíduo estão disponíveis. 

O escopo deste projeto consiste em, a partir dos modelos de crescimentos estudados, utilizar-se da estrutura de Efeitos Mistos não-linear --- compondo os dados individuais com a abordagem populacional obtida através dos dados de indivíduos semelhantes --- para se obter um modelo \emph{a priori} da estimativa da curva de crescimento e da idade de um tumor. O efeito misto se obtém compondo a informação da população como a priori com as poucas medidas disponíveis do indivíduo. A adoção de distribuição Fat-tail permite estimar os valores da curva de crescimento e da idade do tumor para cada indivíduo com uma boa acuidade, vide [1].  Uma comparação será feita entre os modelos de crescimento citados ou outros e também com os ajustes obtidos através com a distribuição Gaussiana tradicional. Pretende-se verificar se a maior variabilidade de eventos capturada por essas distribuições e modelos irá gerar resultados dos parâmetros de crescimento e da idade do tumor mais precisos e que acenem para um possível uso clínico.


# Descrição do projeto

# Conteúdo do repositório

# Referências
[1] R. F. S. Marcos R. Fernandes and J. B. R. do Val. Robust Mixed-Effect Estimation of Tumor Growth and Age based on Gompertz Model. IEEE Control Systems Letters, 7:31–36, 2023
