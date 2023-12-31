# Iniciação Científica - Estimativa de crescimento e idade de tumores através da técnica de Efeitos Mistos da população e indivíduo

# Resumo

O tema deste projeto é a obtenção de estimativas da forma de crescimento e da idade de tumores utilizando modelos de crescimentos de tumores com modelos agregados representando volume ou número de celúlas, como os Modelos Logístico ou o de Gompertz. Utilizamos uma estrutura de Estimativa Bayesiana conhecida como estimativa de Efeitos Mistos não-linear (nonlinear Mixed-Effect Estimation) com distribuição $\textit{Fat-tail}$. A técnica de Efeitos Mistos permite compor os poucos dados disponíveis para cada indivíduo de medidas quando agrega-se populações semelhantes sujeitas ao mesmo tipo de carcinoma.

Utilizamos uma distribuição mista Laplaciana-Gaussiana para alcançar maior robustez em vista da dispersão de dados, em geral mal representados pelos modelos convencionais Gaussianos e métodos associados como os de mínimos quadrados e suas variações. Modelos estatísticos mais robustos e flexíveis são críticos neste tipo de estimativas nas quais poucos dados para um indivíduo estão disponíveis. 

O escopo deste projeto consistiu em compreender experimentar com a estrutura de Efeitos Mistos não-linear, a partir dos modelos de crescimentos estudados, compondo os dados individuais com a abordagem populacional obtida através dos dados de indivíduos semelhantes, Os dados populacionais agregados permitem se obter o modelo estatístico a priori da estimativa da curva de crescimento e da idade de um tumor. O efeito misto se obtém compondo a informação da população como a priori com as poucas medidas disponíveis do indivíduo. A adoção de distribuição de caudas pesadas ($\textit{Fat-tail}$) permite estimar os valores da curva de crescimento e da idade do tumor para cada indivíduo com boa acuidade, como será desenvolvido neste relatório a partir de estudos anteriores, vide [1].  Uma comparação é feita entre os modelos de crescimento existentes e também com os ajustes obtidos através da distribuição Gaussiana tradicional. Pretende-se verificar dois aspectos principais: se a maior variabilidade de eventos capturada pela distribuição mista Laplace-Gaussiana e os modelos de crescimento irão ter impacto na melhoria das estimativas dos parâmetros de crescimento e da idade do tumor resultando em modelos mais precisos e aderentes aos dados. O segundo resultado decorrente é ganhar segurança que possam acenar para um possível uso clínico.


# Descrição do projeto

Foi utilizado como base de dados presente no artigo [2], que consiste em medidas de crescimento de tumores em uma população de 66 ratos com o carcinoma LM2-4LUC+ implantadas. O volume inicial foi de $1{mm^3}$ implantados em todos os ratos da população. 

Para caracterizar o crescimento dos tumores, empregou-se o modelo de Richards, mostrado abaixo, com o intuito de verificar a eficiência da aplicação dos conceitos de Efeitos Mistos não-linear e da distribuição Normal-Laplace para a estimação dos parâmetros individuais, comparado com o uso da tradicional distribuição Gaussiana.

$$
    V(t) = \frac{V_{0}.K}{[V_{0}^{\nu}+(K^{\nu}-V_{0}^{\nu})e^{-r\nu t}]^{1/\nu}}
$$

Além da comparação com o modelo da distribuição Gaussiana puro, foi feita uma comparação com os resultados obtidos com o modelo de Gompertz [1], com o intuito de verificar se os resultados do Richards possam ser melhores ou iguais ao modelo de Gompertz.

Para a confecção dos códigos, foi utilizado o Matlab R2017a e foi tomado como base os códigos presentes em [1], mas adaptando os códigos de priori populacional e de estimação individual para o modelo de Richards.

# Conteúdo do repositório

- README.md;
- Teste da matriz Jacobiana;
- Priori Populacional;
- Estimação Individual;
- Função EVIU;
- Função de Estimação Populacional;
- Função do algoritmo "Coordinated Descent";
- Rotina para gerar os gráficos utilizados para a comparação entre os modelos de Gompertz e de Richards;
- Grafícos Obtidos.

# Referências
[1] R. F. S. Marcos R. Fernandes and J. B. R. do Val. Robust Mixed-Effect Estimation of Tumor Growth and Age based on Gompertz Model. IEEE Control Systems Letters, 7:31–36, 2023

[2] C. V. et al. Population modeling of tumor growth curves and the reduced gompertz model improve prediction of the age of experimental tumors. 16:e1007178
