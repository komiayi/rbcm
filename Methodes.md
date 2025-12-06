Dans le cadre de l’analyse de médiation, un chercheur peut, en toute
connaissance de cause, choisir de se concentrer sur un médiateur
d’intérêt particulier parmi plusieurs médiateurs présents entre une
exposition et une réponse. Cette approche, qui revient à mener une
analyse de médiation simple dans un contexte de médiation multiple, vise
à estimer l’effet indirect de l’exposition à travers le médiateur cible,
ainsi que l’effet direct complémentaire, représentant l’ensemble des
autres voies d’influence, y compris celles passant par les autres
médiateurs.

Comme souligné dans l’introduction, dans le contexte de la Figure , où
le médiateur d’intérêt partage une cause commune non mesurée avec
d’autres médiateurs, l’analyse classique de médiation conduit à des
estimateurs biaisés des effets direct et indirect du médiateur
d’intérêt. En effet, cette cause commune agit comme un confondant non
mesuré entre le médiateur d’intérêt et la réponse, violant l’hypothèse
de SIA. Le Chapitre 2 qui suit introduit une méthodologie permettant de
contourner ce problème.

Ce chapitre propose une méthodologie pour effectuer une analyse de
médiation simple pour un médiateur particulier, en tenant compte de la
présence d’un second médiateur présumé non lié causalement au premier.
Dans ce contexte, nous supposons donc que ces deux médiateurs sont
mesurés.

# Mise en contexte

Considérons le vecteur **X**, qui représente les variables de
préexposition, *T* une variable binaire dite d’exposition, où *T* = 1
indique une exposition et *T* = 0 une absence d’exposition, ainsi que
les médiateurs *M*<sub>1</sub>, *M*<sub>2</sub> impliqués dans le
mécanisme causal. La variable *U* représente une cause commune non
mesurée des médiateurs, et *Y* la réponse. La structure causale entre
ces éléments est représentée dans le diagramme causal de la Figure
[Figure 1](#fig1).

![](figure_tikz.png)

Figure 1: DAG d’analyse de médiation avec deux médiateurs disposant
d’une cause commune.

Soit *M*<sub>1</sub> le médiateur auquel on s’intéresse. Nous notons par
**M** le vecteur des deux médiateurs observés. Dans le diagramme causal
de la Figure [Figure 1](#fig1), *M*<sub>1</sub> est influencé par **X**,
*T* et la cause commune non mesurée *U*. Cette dernière affecte
simultanément *M*<sub>1</sub> et le médiateur *M*<sub>2</sub>, ce qui
crée une dépendance entre *M*<sub>1</sub> et *M*<sub>2</sub>.
S’intéresser à un médiateur en particulier conduit à une analyse
restreinte, où nous nous concentrons sur le médiateur *M*<sub>1</sub>,
comme illustré dans le diagramme de la Figure [Figure 2](#fig2),
représentant les relations causales directes et indirectes entre *T*,
*M*<sub>1</sub> et *Y*. Les chemins *T* → *Y* et
*T* → *M*<sub>2</sub> → *Y* sont les chemins contribuant à l’effet
direct de *T* sur *Y* en considérant *M*<sub>1</sub> comme médiateur
ciblé, tandis que *T* → *M*<sub>1</sub> → *Y* est le chemin indirect via
*M*<sub>1</sub>. Outre le lien direct de *M*<sub>1</sub> sur *Y*, un
chemin non causal existe entre *M*<sub>1</sub> et *Y*, qui passe par
*U*, soit *M*<sub>1</sub> ← *U* → *M*<sub>2</sub> → *Y*. Ce chemin
représente un chemin de confusion entre *M*<sub>1</sub> et *Y* par la
présence de la cause commune *U* non mesurée, posant ainsi problème pour
l’utilisation de l’approche de médiation simple standard.

![](figure2_tikz.png)

Figure 2: DAG d’analyse simple de *M*<sub>1</sub>.

Notre objectif est donc d’identifier et d’estimer les effets naturels
simples, direct et indirect de *M*<sub>1</sub>, en contournant la
difficulté liée à la présence de la variable non mesurée *U*.

Dans la littérature, une solution classique pour mitiger le biais de
confusion dans une analyse est d’ajuster ou stratifier sur les
confondants. Cependant, comme *U* n’est pas mesurée, nous ne pouvons pas
adopter cette stratégie. L’ajustement sur *M*<sub>2</sub> pourrait
effectivement éliminer la confusion entre le médiateur d’intérêt et la
réponse qui découle de l’ouverture du chemin
*M*<sub>1</sub> ← *U* → *M*<sub>2</sub> → *Y*. Cependant, cette
stratégie mène à un autre problème. On bloque le chemin causal
*T* → *M*<sub>2</sub> → *Y*, ce qui engendre un biais dans l’estimation
de l’effet direct de *T* sur *Y* avec *M*<sub>1</sub> comme médiateur
d’intérêt. Cela crée un dilemme, car il faut à la fois ajuster et ne pas
ajuster pour *M*<sub>2</sub>.

Pour résoudre ce dilemme, notre approche dans un premier temps consiste
à inclure le médiateur *M*<sub>2</sub> dans l’analyse en adaptant
l’hypothèse de composition de Halpern and Pearl (2005) afin de tenir
compte de cette structure causale complexe.

## Incorporation du médiateur *M*<sub>2</sub> dans l’analyse de médiation simple

### Définition des effets direct et indirect dans l’analyse de médiation simple

Dans cette section, nous incorporons le médiateur *M*<sub>2</sub> dans
l’analyse de médiation simple afin d’évaluer les effets direct et
indirect de *M*<sub>1</sub> en tenant compte de la dépendance entre les
deux médiateurs. Soit **M**(*t*) le vecteur des valeurs potentielles des
médiateurs lorsque l’exposition est fixée à *t*. De même,
*Y*(*t*,**M**(*t*′)) désigne la valeur potentielle de la réponse lorsque
l’exposition est fixée à *t* et que les médiateurs prennent leurs
valeurs potentielles sous une exposition *t*′. Enfin,
*Y*(*t*,*M*<sub>1</sub>(*t*′),*M*<sub>2</sub>(*t*″)) représente la
valeur potentielle de la réponse lorsque *M*<sub>1</sub> et
*M*<sub>2</sub> sont fixés respectivement à leurs valeurs potentielles
sous les expositions *t*′ et *t*″.

L’hypothèse de composition s’exprime ainsi :
**(H**<sub>0</sub>**)**  *Y*(*t*,*M*<sub>1</sub>(*t*′)) = *Y*(*t*,*M*<sub>1</sub>(*t*′),**M**<sub>2</sub>(*t*)), ∀ *t*, *t*′.
Cette hypothèse énonce que la valeur potentielle de *Y* lorsque *T* est
fixé à *t* et *M*<sub>1</sub> à sa valeur potentielle sous une autre
exposition *t*′ est égale à la valeur potentielle de *Y* lorsque *T* est
fixé à *t*, *M*<sub>1</sub> à sa valeur sous *t*′, et le vecteur
*M*<sub>2</sub> à sa valeur potentielle sous *t*. En d’autres termes,
l’intervention sur *M*<sub>2</sub> n’a aucun effet supplémentaire sur la
réponse potentielle.

Sous cette hypothèse, les effets direct et indirect de *M*<sub>1</sub>
sont définis comme suit :
$$
\zeta= \mathbb{E}\left\\Y(1,M\_1(0),M\_2(1))-Y(0,M\_1(0),M\_2(0))\right\\,\\
\delta= \mathbb{E}\left\\Y(1,M\_1(1),M\_2(1))-Y(1,M\_1(0),M\_2(1))\right\\.
$$
L’effet indirect *δ* mesure l’impact de *T* sur *Y* à travers les
variations spécifiques de *M*<sub>1</sub>, tout en maintenant le
médiateur *M*<sub>2</sub> fixé à sa valeur potentielle sous l’exposition
*T* = 1. Il est égal à l’effet indirect individuel de ce médiateur dans
le contexte d’une analyse de médiation multiple. En revanche, l’effet
direct *ζ* avec *M*<sub>1</sub> comme médiateur ciblé capture l’impact
de *T* sur *Y*, tout en contrôlant uniquement *M*<sub>1</sub> à une
valeur fixée. Il est à noter que le changement d’état de l’exposition
*T* affecte aussi *M*<sub>2</sub>, ce qui fait que l’effet direct lié au
médiateur *M*<sub>1</sub> intègre aussi l’effet du médiateur
*M*<sub>2</sub> sur la réponse.

Ainsi, l’effet direct *ζ* peut être décomposé en deux parties
distinctes. Tout d’abord, l’impact de l’exposition sur *Y*, en
maintenant *M*<sub>2</sub> fixé à sa valeur potentielle sous *T* = 0 :
𝔼{*Y*(1,*M*<sub>1</sub>(0),*M*<sub>2</sub>(0)) − *Y*(0,*M*<sub>1</sub>(0),*M*<sub>2</sub>(0))}.
Cette expression correspond à l’effet direct joint dans le cadre d’une
analyse multiple avec les deux médiateurs. Puis, l’impact de
l’exposition sur *Y* par le seul changement d’état de *M*<sub>2</sub> :
𝔼{*Y*(1,*M*<sub>1</sub>(0),*M*<sub>2</sub>(1)) − *Y*(1,*M*<sub>1</sub>(0),*M*<sub>2</sub>(0))},
qui correspond à l’effet indirect à travers *M*<sub>2</sub> uniquement
(Jérolon et al. 2021).

Ainsi, avec l’hypothèse de composition, les effets naturels *δ* et *ζ*
sont déterminés comme étant des fonctions de l’effet direct joint et des
effets indirects individuels d’une analyse multiple. Dans ces effets,
les réponses potentielles *Y*(1,*M*<sub>1</sub>(0),*M*<sub>2</sub>(1))
et *Y*(1,*M*<sub>1</sub>(0),*M*<sub>2</sub>(0)), nécessaires pour
calculer *δ* et *ζ* ne sont pas observables dans les données. Par
conséquent, ces effets naturels ne sont pas identifiables sans
introduire des hypothèses supplémentaires.

### Identification des effets naturels *δ* et *ζ*

Les effets naturels simples *δ* et *ζ*, définis étant comme des
fonctions de l’effet direct joint et des effets indirects individuels
dans une analyse multiple, reposent pour leur identification sur
l’hypothèse SIMMA (Jérolon et al. 2021).

Sous cette l’hypothèse, les effets direct et indirect peuvent être
partiellement identifiés de manière non paramétrique à l’aide des
formules suivantes :
$$
    \zeta=\int\_{\mathcal{X}}\Bigg(\int\_{\mathcal{M}} \mathbb{E}\[Y|T=1,M\_1=m\_1,M\_2=m\_2,\boldsymbol{X}=\boldsymbol{x}\]dF\_{(M\_1(0),M\_2(1))|\boldsymbol{X}=\boldsymbol{x}}(m\_1,m\_2)\nonumber\\
    \qquad\qquad-\int\_{\mathcal{M}} \mathbb{E}\[Y|T=0,M\_1=m\_1,M\_2=m\_2,\boldsymbol{X}=\boldsymbol{x}\]dF\_{(M\_1,M\_2)|T=0,\boldsymbol{X}=\boldsymbol{x}}(m\_1,m\_2)\Bigg)dF\_{\boldsymbol{X}}(\boldsymbol{x}),
$$

$$
\delta=\int\_{\mathcal{X}}\Bigg(\int\_{\mathcal{M}} \mathbb{E}\[Y|T=1,M\_1=m\_1,M\_2=m\_2,\boldsymbol{X}=\boldsymbol{x}\]dF\_{(M\_1,M\_2)|T=1,\boldsymbol{X}=\boldsymbol{x}}(m\_1,m\_2)\nonumber\\
\qquad\qquad-\int\_{\mathcal{M}} \mathbb{E}\[Y|T=1,M\_1=m\_1,M\_2=m\_2,\boldsymbol{X}=\boldsymbol{x}\]dF\_{(M\_1(0),M\_2(1))|\boldsymbol{X}=\boldsymbol{x}}(m\_1,m\_2)\Bigg)dF\_{\boldsymbol{X}}(\boldsymbol{x}),
$$
où *d**F*<sub>**X**</sub> est la mesure de probabilité de la
distribution de **X** et
*d**F*<sub>**M**|*T* = *t*′, **X** = **x**</sub> la mesure de
probabilité conditionnelle de **M** étant donné que *T* = *t*′ et
**X** = **x**.

L’identification complète de ces effets simples nécessite cependant des
hypothèses supplémentaires sur la loi conjointe de
(*M*<sub>1</sub>(0),*M*<sub>2</sub>(1))|**X** = **x**. Dans la section
suivante, nous détaillons une approche par régression pour des variables
continues, et l’identification des estimateurs des effets *ζ* et *δ*.

# Références

Halpern, Joseph Y., and Judea Pearl. 2005. “Causes and Explanations: A
Structural-Model Approach. Part i: Causes.” *British Journal for the
Philosophy of Science* 56 (4): 843–87.
<https://doi.org/10.1093/bjps/axi147>.

Jérolon, Allan, Laura Baglietto, Etienne Birmelé, Flora Alarcon, and
Vittorio Perduca. 2021. “Causal Mediation Analysis in Presence of
Multiple Mediators Uncausally Related.” *The International Journal of
Biostatistics* 17 (2): 191–221.
