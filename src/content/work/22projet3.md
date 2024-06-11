---
title: Projet 3
publishDate: 2020-03-04 00:00:00
img: /assets/stock-3.jpg
img_alt: Pearls of silky soft white cotton, bubble up under vibrant lighting
description: |
  Étude de stratifiés thermo-elastiquement stables sur pâle d'hélicopter
tags:
  - Design
  - Dev
  - Branding
---

<script>
MathJax = {
  	tex: {
    	inlineMath: [['$', '$'], ['\\(', '\\)']]
  	},
  	svg: {
    fontCache: 'global'
  	}
};
</script>

## Description du problème

Dans la conception de structures composites il peut être utile de réaliser des stratifiés couplés élastiquement afin de produire des changements de forme lors d’un chargement mécanique : par exemple, le \textbf{couplage traction-torsion} peut être exploité pour construire des aubes ou pales à pas adaptatif, la torsion intervenant lors de la rotation des pales/aubes par effet de la sollicitation produite par effet centrifuge (plus la vitesse de rotation est élevée et plus important sera l’angle de torsion).
\begin{figure}[H]
    \centering
    \includegraphics[width=60mm]{./images/couplage.png}
    \caption{Couplage traction/torsion}
\end{figure}

Néanmoins, lors de la fabrication, les pièces composites subissent des variations de température et les couplages peuvent produire des déformations/courbures d’origine thermo-élastique, ce qui n’est pas souhaitable pour garantir la stabilité de la géométrie de ces pièces.\medbreak
\noindent
Il est donc nécessaire de concevoir des stratifiés ayant un couplage thermo-élastique nul \\(\textbf{V = O}). On parle dans ce cas de \textbf{\textit{stratifiés thermo-élastiquement stables.}}
On se propose d’étudier ici comment concevoir des stratifiés de ce type et d’analyser leurs propriétés. Ensuite on considérera une famille particulière de stratifiés thermo-élastiquement stables et on en étudiera le comportement.
On considère ici le cas de stratifiés constitués de couches identiques (même matériau et même épaisseur) : les paramètres de conception seront donc le nombre de couches n et la séquence d’angles d’orientation [$\delta$k] (k = 1,..., n).
On utilisera la représentation polaire pour analyser le comportement de stratifiés, partout où il est possible.
On notera « CB » les paramètres qui décrivent les propriétés de la couche de base et par des symboles $\bar{}$, $\hat{}$ et $\tilde{}$ les termes en lien respectivement à $\textbf{A}$,  $\textbf{B}$ et  $\textbf{D}$.

\newpage
\subsection*{Propriétés générales}
\noindent
1. \textit{Expliquer pourquoi si V = O pour un stratifié à couches identiques, alors le tenseur de couplage élastique B est à symétrie du carré.}\smallbreak
\noindent
L'objecif est de montrer que le tenseur \textbf{B} d'ordre 4 de couplage élastique est à symétrie du carré en considérant le cas de figure de stratifié à couches identiques et d'un couplage thermo-élastique nul i.e \textbf{V=O}.\smallbreak
\noindent
\underline{Remarque} : le tenseur \textbf{V} permet de coupler les effets de contributions de courbures due à une élévation uniforme de température et les effets de contributions menbranaires due à un gradient thermique dans l'épaisseur du stratifié.\smallbreak
\noindent
Rappelons d'abord la signification d'un stratifié à couches identiques, il correspond à des plis constitués du même matériau et de même épaisseur.\smallbreak
\noindent
On part alors de l'expression de \textbf{V}:
\begin{equation*}
    \textbf{V}=\frac{1}{2}\sum_{k=1}^{N}{\boldsymbol{\gamma}}(\delta_k)(z^2_k-z^2_{k-1})
\end{equation*}
avec $\boldsymbol{\gamma}(\delta_k)$ le tenseur de contrainte thermique par unité de température (Mpa.$^oC^{-1}$),\\
que l'ont particularise pour des couches identiques :
\begin{equation*}
    \textbf{V}=\frac{1}{2}\frac{h^2}{N^2}\sum_{k=1}^{N}b_k{\boldsymbol{\gamma}}(\delta_k)
\end{equation*}
avec le coefficient $b_k = 2k-N-1$\\
N et h nombre et épaisseur de couche total respectivement.\smallbreak
\noindent
Nous pouvons alors exprimer ce tenseur sous forme cartesienne ou polaire. Pour des raisons de simplifications on écrit le tenseur \textbf{V} sous sa forme polaire à l'aide de la représentation des cercles de Mohr (deux paramètres invariants T,R partie sphérique et déviatorique et un angle $\Phi$ de direction principal de contrainte qui varie avec le repère).\\
La décomposition polaire la plus générale de \textbf{V}  nous donne:
\begin{equation*}
    T_v=\frac{1}{2}\big(\sum_{k=1}^{N}(z^2_k-z^2_{k-1})T_{\gamma}\big) \;\;\;;\;\;\; R_v e^{2i\Phi_v} = \frac{1}{2}\big(\sum_{k=1}^{N}(z^2_k-z^2_{k-1})R_{\gamma}e^{2i(\Phi_{\gamma}+\delta_k)}\big)
\end{equation*}
ont particularise pour des couches identiques :\\
\textbullet \quad $T_{\gamma}$, $R_{\gamma}$ et $e^{2i\Phi_{\gamma}}$ sont les mêmes pour toutes les couches (constantes).\\
\textbullet \quad ($z^2_k-z^2_{k-1}$) $\rightarrow$ $\frac{h^2}{N^2}b_k$ et $b_{ii} = 0$\\
Soit,
\begin{equation*}
    T_v= \frac{1}{2}\frac{h^2}{N^2}T_{\gamma}\big(\sum_{k=1}^{N}b_k\big) = 0 \;\;\;;\;\;\; R_v e^{2i\Phi_v} = \frac{1}{2}\frac{h^2}{N^2}R_{\gamma}e^{2i\Phi_{\gamma}}\big(\sum_{k=1}^{N}b_k e^{2i\delta_k}\big)
\end{equation*}
Ainsi, un stratifié à couche identique est thermo-élastiquement stable, i.e \textbf{V = 0}, si et seulement si 
\begin{equation}
    R_v e^{2i\Phi_v} = 0 \;\;\; \Leftrightarrow \;\;\; \fcolorbox{red}{Tan}{$\sum_{k=1}^{N}b_ke^{2i\delta_k} = 0$} \;\; ou \;\; \fcolorbox{red}{SpringGreen}{$R_{\gamma} = 0$} 
    \label{condition}
\end{equation}
\noindent
\underline{Conclusion} : pour avoir \textbf{V=0} dans un stratifié à couche identique il faut satisfaire au moins l'une des deux conditions [\ref{condition}] ci-dessus, soit:\\
\noindent
\fcolorbox{black}{Tan}{\parbox{14cm}{\textbf{1)} La première équation est une condition géométrique, elle dépend de la combinaison des angles (séquences) donné au plis et du coefficient $b_k$. Cette condition complexe donne en réalité deux conditions, une venant de la partie réel et une autre de la partie imaginaire.}}\smallbreak
\noindent
\fcolorbox{black}{SpringGreen}{\parbox{14cm}{\textbf{2)} La deuxième est une condition matérielle, $R_{\gamma} = 0$ signifie que $\boldsymbol{\gamma}$ est purement sphérique. Pour obtenir cela, il faut supposer que les propriétés du matériau de la couche de base est à symétrie du carré (ou cas trivial i.e isotrope), ainsi $R_1^{CB} = 0$ et impliquant $\hat{R}_1 = \tilde{R}_1 = \bar{R}_1 = 0$, ce qui nous ramène à la propriété de symétrie du carré pour \textbf{A},\textbf{B} et \textbf{D}.}} \smallbreak
\noindent
Considérons dans la suite la première condition afin de satisfaire la propriété de symétrie du carré seulement pour \textbf{B}.\smallbreak
\noindent
Le tenseur \textbf{B} d'ordre 4 possède par ailleurs 4 composante polaire. On effectue le même déroulement que précédemment.
\noindent
On part alors de l'expression de \textbf{B}:
\begin{equation*}
    \textbf{B}=\frac{1}{2}\sum_{k=1}^{N}\textbf{Q}(\delta_k)(z^2_k-z^2_{k-1})
\end{equation*}
avec $\textbf{Q}(\delta_k)$ le tenseur d'ordre 4 d'élasticité en 2D (Mpa)\\
que l'ont particularise pour des couches identiques :
\begin{equation*}
    \textbf{B}=\frac{1}{2}\frac{h^2}{N^2}\sum_{k=1}^{N}b_k\textbf{Q}(\delta_k)
\end{equation*}
\noindent
La décomposition polaire de \textbf{B} est donné par $\hat{T}_0$, $\hat{T}_1$, $\hat{R}_0e^{4i\hat{\Phi}_0}$ et $\hat{R}_1e^{2i\hat{\Phi}_1}$ que l'ont écrit sous forme particularisé pour des couches identiques.\\ Comme pour le calcul de $T_v$ on a $\hat{T}_0$ = $\hat{T}_1$ = 0 car $\sum_{k=1}^{N}b_k = 0$.\\
\noindent
Pour $\hat{R}_0e^{4i\hat{\Phi}_0}= \frac{1}{2}\big(\sum_{k=1}^{N}R_{0k}e^{4i(\Phi_{0k}+\delta_k)}(z^2_k-z^2_{k-1})\big)$  on a,
\begin{align*}
   \hat{R}_0e^{4i\hat{\Phi}_0}&= \frac{1}{2}\frac{h^2}{N^2}\big(\sum_{k=1}^{N}b_k{R}_0^{CB}e^{4i(\Phi_0^{CB}+\delta_k)}\big) \\
   &= \frac{1}{2}\frac{h^2}{N^2}{R}_0^{CB}e^{4i\Phi_0^{CB}}\big(\sum_{k=1}^{N}b_ke^{4i\delta_k}\big) 
\end{align*}
\noindent
et pour $\hat{R}_1e^{2i\hat{\Phi}_1}= \frac{1}{2}\big(\sum_{k=1}^{N}R_{1k}e^{2i(\Phi_{1k}+\delta_k)}(z^2_k-z^2_{k-1})\big)$  on a,
\begin{align*}
   \hat{R}_1e^{2i\hat{\Phi}_1}&= \frac{1}{2}\frac{h^2}{N^2}\big(\sum_{k=1}^{N}b_k{R}_1^{CB}e^{2i(\Phi_1^{CB}+\delta_k)}\big) \\
   &= \frac{1}{2}\frac{h^2}{N^2}{R}_1^{CB}e^{2i\Phi_1^{CB}}\big(\sum_{k=1}^{N}b_ke^{2i\delta_k}\big) 
\end{align*}
\noindent
Or d'après la condition [\ref{condition}] pour que \textbf{V} = 0, il vient:

\begin{equation*}
     \hat{R}_1e^{2i\hat{\Phi}_1}=0 \;\;\; \Leftrightarrow \;\;\; \fcolorbox{red}{Tan}{$\hat{R}_1=0$}
\end{equation*}
Le paramètre polaire $\hat{R}_1$ étant nul, nous avons bien la symétrie du carré seulement pour le tenseur \textbf{B}, en donnant une orientation adéquate de la séquences d'angles. Cela est satisfait par exemple dans la configuration de bi-couche à $0^o$ et $90^o$ (cross-ply).\smallbreak
\noindent
\textit{Quelles sont donc les composantes non nulles de \textbf{B} ?}\smallbreak
\noindent
D'après la question précédente, la composante polaire non nulle de \textbf{B} est
\begin{equation*}
\hat{R}_0e^{4i\hat{\Phi}_0}   
\end{equation*}donc,
\begin{equation*}
\fcolorbox{red}{white}{$\hat{R}_0cos{4\hat{\Phi}_0}$} \;\;\; et \;\;\; \fcolorbox{red}{white}{$\hat{R}_0sin{4\hat{\Phi}_0}$}  
\end{equation*}
\smallbreak
\noindent
\textit{Ecrire dans ce cas la représentation cartésienne de B et analyser les divers types de couplages élastiques qui peuvent se manifester.}\\
\smallbreak
\noindent
Dans le cas générique les composantes polaire de \textit{Verchery} pour un tenseur d'ordre 4 s'écrivent:
\begin{align*}
L_{1111} &= \quad T_0 + 2T_1 + R_0cos4\Phi_0 + 4R_1cos2\Phi_1\\
L_{1112} &= \quad \quad  \quad  \quad \quad \quad R_0sin4\Phi_0 + 2R_1sin2\Phi_1\\
L_{1122} &=-T_0 + 2T_1 - R_0cos4\Phi_0\\
L_{1212} &= \quad T_0  \quad \quad \quad - R_0cos4\Phi_0\\
L_{2212} &= \quad \quad \quad  \quad \quad \;- R_0sin4\Phi_0 + 2R_1sin2\Phi_1\\
L_{2222} &= \quad T_0 + 2T_1 + R_0cos4\Phi_0 - 4R_1cos2\Phi_1
\end{align*}
\noindent
et dans notre cas où $\hat{T}_0=\hat{T}_1=\hat{R}_1=0$ le tenseur \textbf{B} s'écrit alors:
\begin{equation*}
\textbf{B} =
\begin{pmatrix}
\hat{R}_0cos4\hat{\Phi}_0 & -\hat{R}_0cos4\hat{\Phi}_0 & \hat{R}_0sin4\hat{\Phi}_0  \\
-\hat{R}_0cos4\hat{\Phi}_0  & \hat{R}_0cos4\hat{\Phi}_0 & -\hat{R}_0sin4\hat{\Phi}_0 \\
\hat{R}_0sin4\hat{\Phi}_0  & -\hat{R}_0sin4\hat{\Phi}_0 & -\hat{R}_0cos4\hat{\Phi}_0 
\end{pmatrix}
\end{equation*}

\noindent
Pour analyser les divers couplages élastiques, écrivons la loi de comportement thermo-élastique (sans gradient thermique et \textbf{V}=0):

\begin{equation*}
\begin{bmatrix}
\textbf{N}\\
\textbf{M}
\end{bmatrix}=
\begin{bmatrix}
\textbf{A} &\textbf{B} \\
\textbf{B} & \textbf{D} 
\end{bmatrix}
\begin{bmatrix}
\boldsymbol{\mathcal{E}^0}\\
\boldsymbol{\mathcal{K}}
\end{bmatrix}-\Delta T_0
\begin{bmatrix}
\textbf{U}\\
\textbf{0}
\end{bmatrix}
\end{equation*}

\begin{equation*}
\Rightarrow \; \; \left\{
    \begin{array}{ll}
        \textbf{N} = \textbf{A}\boldsymbol{\mathcal{E}}^0 + \textbf{B}\boldsymbol{\mathcal{K}} - T_0\textbf{U} \\
        \textbf{M} = \textbf{B}\boldsymbol{\mathcal{E}}^0 + \textbf{D}\boldsymbol{\mathcal{K}}
    \end{array}
\right.
\end{equation*}
\noindent
Nous analysons le terme non nul $\textbf{B}\boldsymbol{\mathcal{E}}^0$ qui couple les effets d'efforts de membranes et de moments.
\begin{equation*}
\textbf{M}^{coupl} = \textbf{B}\boldsymbol{\mathcal{E}}^0=
\begin{bmatrix}
M_{1}\\
M_{2}\\
M_{6}
\end{bmatrix}=
\begin{bmatrix}
B_{11}\mathcal{E}^0_1 + B_{12}\mathcal{E}^0_2 + B_{16}\mathcal{E}^0_6\\
B_{12}\mathcal{E}^0_1 + B_{22}\mathcal{E}^0_2 + B_{26}\mathcal{E}^0_6\\
B_{16}\mathcal{E}^0_1 + B_{26}\mathcal{E}^0_2 + B_{66}\mathcal{E}^0_6
\end{bmatrix}
\end{equation*}

\noindent
Si on veut obtenir un couplage en "traction-torsion" il faut activer les composantes $B_{16}$ et $B_{26}$ qui lie justement le moment de torsion $M_6$ avec les déformations de membranes longitudinal $\mathcal{E}^0_1$ et $\mathcal{E}^0_2$ (en réalité il faut d'abord inverser la loi pour parler d'un tel couplage cf. question 8). \\
A titre d'exemple, les constructeurs de pales d'hélicoptères exploite cette effet de traction-torsion. En effet, lors de la rotation des pales la force centrifuge exercé de manière radial dans le plan des pales va déclencher une torsion de ces pales, cette effet est d'un grand intérêt aérodynamique.\smallbreak
\noindent
Pour cela on peut imposer $sin4\hat{\Phi}_0 = \pm 1$ avec ${\Phi}_0 = \frac{\pi}{8}$ ce qui implique automatiquement la correspondance $cos4\hat{\Phi}_0 = 0$ 
% Dans le cas général de matériau orthotrope $\Phi_1 = 0^o$ dans la couche de base et que $\Phi_1 = 0^o$-$\Phi_1 = k\.rac{\pi}{4}$, on prend ici $k=0$. Soit, $sin4\hat{\Phi}_0 = \pm 1$ avec ce qui implique automatiquement la correspondance $cos4\hat{\Phi}_0 = 0$ 
\begin{equation*}
\textbf{B}=
\begin{bmatrix}
0 & 0 & \pm\hat{R}_0  \\
0 & 0 & \mp\hat{R}_0  \\
\pm\hat{R}_0 & \mp\hat{R}_0 & 0
\end{bmatrix}
\Rightarrow \textbf{Extension-Twisting and Shearing-Bending;(ETSB)}
\end{equation*}

\noindent
ou alors si on veut annuler entièrement les effets de traction-torsion, il suffit de tourner la séquence des plis dans le stratifié de telle manière à obtenir $cos4\hat{\Phi}_0 = \pm 1$ ce qui implique que $sin4\hat{\Phi}_0 = 0$, en prenant par exemple ${\Phi}_0 = {\Phi}_1 = 0$ (orthotropie ordinaire)

\begin{equation*}
\textbf{B}=
\begin{bmatrix}
\pm\hat{R}_0  & \mp\hat{R}_0  & 0  \\
\mp\hat{R}_0  & \pm\hat{R}_0  & 0  \\
0 & 0 & \mp\hat{R}_0 
\end{bmatrix}
\Rightarrow \textbf{
Extension-Bending and Shearing-Twisting;(EBST)}
\end{equation*}

\noindent
on alors ici un couplage en flexion-traction et cisaillement-torsion.

\noindent
2. \textit{Un autre aspect important de la stabilité thermo-élastique est lié à la forme du tenseur \textbf{U} : celui-ci mesure les efforts thermiques dans le stratifié induits par une élévation uniforme de température. Est-ce que le tenseur \textbf{U} peut être nul ?}\smallbreak
\noindent
On écrit pour cela le tenseur \textbf{U} sous sa forme polaire la plus générale avec $T_u$ et $R_u e^{2i\Phi_u}$ :
\begin{equation*}
    T_u=\big(\sum_{k=1}^{N}(z_k-z_{k-1})T_{\gamma}\big) \;\;\;;\;\;\; R_u e^{2i\Phi_u} = \big(\sum_{k=1}^{N}(z_k-z_{k-1})R_{\gamma}e^{2i(\Phi_{\gamma}+\delta_k)}\big)
\end{equation*}
ont particularise pour des couches identiques :\\
\begin{equation*}
    T_u= \frac{h}{N}\big(\sum_{k=1}^{N}T_{\gamma}\big) = hT_{\gamma}\;\;\;;\;\;\; R_u e^{2i\Phi_u} = \frac{h}{N}\big(\sum_{k=1}^{N}R_{\gamma}e^{2i(\Phi_{\gamma}+\delta_k}\big)=\frac{h}{N}R_{\gamma}e^{2i\Phi_{\gamma}}\big(\sum_{k=1}^{N}e^{2i\delta_k}\big)
\end{equation*}
\noindent
Ainsi, un stratifié à couche identique donne \textbf{U = 0}, si et seulement si,

\begin{equation*}
\left\{
    \begin{array}{ll}
        T_u = 0  \\
        R_u e^{2i\Phi_u} = 0 
    \end{array}
\right.
 \Leftrightarrow \; \; \left\{
    \begin{array}{ll}
        T_{\gamma} = 0  \\
        R_{\gamma} = 0 \quad ou \quad \sum_{k=1}^{N}e^{2i\delta_k}=0 
    \end{array}
\right.
\end{equation*} 
\noindent

% Comme déjà dit, $R_{\gamma} = 0$ signifie que $\boldsymbol{\gamma}$ est purement sphérique, c'est à dire que les propriétés du matériau de la couche de base du stratifié est à symétrie du carré ($R_1^{CB} = 0$) ou cas trivial i.e isotrope.


\noindent
Pour voir comment annuler $T_{\gamma}$ et $R_{\gamma}$, donnons avant tout leurs expressions à l'aide de la représentation polaire sur les cercles de Mohr.\\ En effet $T_{\gamma}$ est l'invariant sphérique qui représente l'abscisse du centre du cercle de Mohr et $R_{\gamma}$ l'invariant déviatorique qui représente le rayon du cercle de Mohr, soit,

\begin{equation*}
T_{\gamma} = \frac{\gamma_1 + \gamma_2}{2} \;\;\; et \;\;\; R_{\gamma} = \sqrt{\Big(\frac{\gamma_1 - \gamma_2}{2}\Big)^2}
\end{equation*}

\noindent
Calculons ensuite les composante du tenseur $\boldsymbol{\gamma}$ tel que $\boldsymbol{\gamma}$=\textbf{Q}$\boldsymbol{\alpha}$, en prenant un comportement au moins de type orthotrope pour le tenseur de dilatation thermique $\boldsymbol{\alpha}$ et le tenseur de rigidité \textbf{Q} exprimé dans la couche de base (i.e dans ces axes d'orthotropie $\rightarrow$ $Q_{16}$ = $Q_{26}$ = 0):
\begin{equation*}
\boldsymbol{\gamma} = 
\begin{bmatrix}
\gamma_1 \\
\gamma_2 \\
\gamma_6 
\end{bmatrix}=
\begin{bmatrix}
Q_{11} & Q_{12} & 0 \\
Q_{12} & Q_{22} & 0 \\
0 & 0 & Q_{66} 
\end{bmatrix}
\begin{bmatrix}
\alpha_1 \\
\alpha_2 \\
0
\end{bmatrix}=
\begin{bmatrix}
Q_{11}\alpha_1 + Q_{12}\alpha_2 \\
Q_{12}\alpha_1 + Q_{22}\alpha_2 \\
0
\end{bmatrix}
\end{equation*}

\noindent
ainsi en remplaçant, 

\begin{equation*}
T_{\gamma} = \frac{1}{2}\big[\alpha_1(Q_{11} + Q_{12}) + \alpha_2(Q_{12} + Q_{22})\big]
\end{equation*}


\begin{equation*}
R_{\gamma} = \frac{1}{2}[\alpha_1(Q_{11} - Q_{12}) + \alpha_2(Q_{12} - Q_{22})] 
\end{equation*}
\noindent
Dans le cas particulier où la couche de base est à symétrie du carré les différentes composantes s'identifient comme suit,

\begin{equation*}
\alpha_1 = \alpha_2 \;\;\; \text{et} \;\;\; Q_{11} = Q_{22} 
\end{equation*}
\noindent
et on obtient 

\begin{equation*}
\fcolorbox{red}{SpringGreen}{$R_{\gamma} = 0$} 
\end{equation*}
\noindent
On alors \textbf{V}=0 (condition matérielle à satisfaire cf. [\ref{condition}]) et aussi \boxed{\textbf{U} \,\text{sphérique}}\smallbreak
\noindent
Et, pour obtenir U nul pour les deux solutions possible il faut que $T_\gamma$ soit nul.\smallbreak

\noindent
Est-il possible d'obtenir  $T_{\gamma} = 0$ ?\\
En utilisant les expressions de la page 4 donnant les composantes polaire pour le tenseur \textbf{B}, il vient,
\begin{align*}
Q_{11} + Q_{12} &= 4T_1 + 4R_1cos2\Phi_1\\ 
Q_{12} + Q_{22} &= 4T_1 - 4R_1cos2\Phi_1 
\end{align*}
Or $\Phi_1$ = 0 pour la couche de base, 
\begin{align*}
Q_{11} + Q_{12} &= 4T_1 + 4R_1\\
Q_{12} + Q_{22} &= 4T_1 - 4R_1
\end{align*}

\noindent
ce qui conduit à,

\begin{equation*}
T_{\gamma} = \frac{1}{2}\big[4\alpha_1(T_1 + R_1) + 4\alpha_2(T_1 - R_1)\big] = 0
\end{equation*}
\noindent
Étant donné que $T_1$, $\alpha_1$, $\alpha_2$ sont strictement positif dans le cas de matériau d'usage commun. De plus la dilatation thermique $\alpha_2$ défini dans le sens de la matrice est supérieur à $\alpha_1$ dans la direction des fibres. Donc la condition à satisfaire est telle que :

\begin{equation*}
\fcolorbox{red}{white}{$R_1 = T_1\frac{(\alpha_1+\alpha_2)}{(\alpha_2R_1 - \alpha_1)} > 0$}
\end{equation*}

\noindent
Les matériaux d'usage commun ne satisfont pas cette condition. Ainsi $T_{\gamma}$ n'est pas nul en générale, même si théoriquement cela est possible, et donc faisable seulement pour des matériaux \textbf{ad-hoc}. Donc \textbf{U} est en générale non nul.

\noindent
\textit{Si on imagine une plaque stratifiée de forme rectangulaire dans un référentiel d’axes x-y, quelle forme du tenseur U garantit que la plaque reste de forme rectangulaire ?}

\noindent
Le tenseur des efforts linéique thermique par unité de température uniforme \textbf{U} mesure l'élévation ou la diminution uniforme de température. Si la composante $\textbf{U}_6$ = $\textbf{U}_{xy}$ est nul alors les directions de \textbf{U} sont alignés dans les directions des axes principales d'orthotropie de la couche.

\begin{equation*}
 \textbf{U}=
   \begin{bmatrix}
U_{xx} \\
U_{yy}\\
0 
\end{bmatrix} 
\Rightarrow \; \; \; \; \textbf{Déformation dans le plan uniquement}
\end{equation*}
Ainsi, la plaque rectangulaire reste rectangulaire.

\noindent
\textit{Quel est le lien avec le tenseur A de rigidité en membrane ?} \smallbreak

\noindent
Le tenseur \textbf{A} d'ordre 4 possède 4 composantes polaires. L'expression de \textbf{A} sous forme polaire est donc donné par $\hat{T}_0$, $\hat{T}_1$, $\hat{R}_0e^{4i\hat{\Phi}_0}$ et $\hat{R}_1e^{2i\hat{\Phi}_1}$ que l'ont écrit sous forme particularisé pour des couches identiques, on a pour $\hat{R}_1e^{2i\hat{\Phi}_1}$:
\begin{equation*}
\hat{R}_1e^{2i\hat{\Phi}_1} =
  \frac{1}{2}\frac{h}{N}{R}_1^{CB}e^{2i\Phi_1^{CB}}\big(\sum_{k=1}^{N}e^{2i\delta_k}\big) 
\end{equation*}
\noindent

\noindent
Le lien entre \textbf{A} et \textbf{U} est qu'ils dépendent tout les deux du terme \fcolorbox{red}{white}{$\sum_{k=1}^{N}e^{2i\delta_k}$}.
\smallbreak
\noindent
\textit{Si on veut que la stabilité thermo-élastique soit étendue au comportement de dilatation thermique dans le plan, quelle condition faut-il imposer ?}

\smallbreak
\noindent
Si maintenant on impose de plus que $U_{xx}$=$U_{yy}$ c'est à dire \textbf{U} sphérique,
\begin{equation*}
 \textbf{U}=
   \begin{bmatrix}
U_{xx} \\
U_{xx}\\
0 
\end{bmatrix} 
\end{equation*}
\noindent
alors la déformation de la plaque rectangulaire est une plaque rectangulaire avec les mêmes proportion de dimension, i.e \textbf{transformation homothétique}.\smallbreak
\noindent
\textit{Quelle conséquence sur le comportement de rigidité en membrane (tenseur A) ?}
\smallbreak
\noindent
Puisque que cette condition de sphéricité pour \textbf{U} est obtenue dans le cas le plus courant comme vue précédemment, lorsque 
\begin{equation*}
\sum_{k=1}^{N}e^{2i\delta_k} = 0
\end{equation*}
\noindent
et d'après l'analyse précèdent le tenseur de rigidité en membrane \textbf{A} \textbf{sera à symétrie du carré}, car $\bar{R}_1$ = 0 dans ce cas.
\smallbreak
\noindent
\textit{A partir des expressions de la CLPT écrites en polaire, écrire l’ensemble des conditions qui doivent être satisfaites afin d’obtenir la stabilité thermo-élastique (à la fois \textbf{V} = O et \textbf{U} sphérique) en fonction des angles de la stratification.}
\smallbreak
\noindent
\ding{226} Pour avoir un stratifié a comportement thermo-élastiquement stable il suffit de prendre une couche de base à $R_1^{CB} = 0$, i.e à symétrie du carré, au quel cas \textbf{V}=0 et \textbf{U} est sphérique quelque soit la séquence, mais A,B et D sont nécessairement aussi à symétrie du carré (limitant) et les tissus sont équilibrés (rare et coûteux par rapport aux UD). \\
\ding{226} L'autre solution est de prendre des couches UD, dans ce cas le stratifié a un comportement thermo-élastiquement stable s'il vérifie la séquence [$\delta_k$] tel que, \textbf{V}=0 $\Leftrightarrow$ $\sum_{k=1}^{N}b_ke^{2i\delta_k} = 0$ et vérifie U sphérique si $\sum_{k=1}^{N}e^{2i\delta_k}=0$.
\smallbreak
\noindent
\textbf{Résumé} :\smallbreak
\noindent
\begin{tabular}{|c|c|c|}
  \hline
          & $1^{ere}$ cas : propriété matérielle  & $2^{eme}$ cas : propriété séquentielle \\
  \hline
  \textbf{V} = 0 & $R_1^{CB} = 0$ &  $\sum_{k=1}^{N}b_ke^{2i\delta_k} = 0$ \\
  \hline
  \textbf{U} sphérique  & $R_1^{CB} = 0$ & $\sum_{k=1}^{N}e^{2i\delta_k} = 0$ \\
  \hline
    \textbf{U} nul  & ad-hoc & $\varnothing$  \\
  \hline
    conséquence & \textbf{A},\textbf{B},\textbf{D} à symétrie du carré & \textbf{B},\textbf{A} à symétrie du carré   \\
  \hline
\end{tabular}

\smallbreak
\noindent
\textit{Quel est le nombre minimal de couches nécessaire pour satisfaire l’ensemble de ces conditions ?}\smallbreak
\noindent
Dans le cadre du $2^{eme}$ cas du tableau nous avons en réalité 4 conditions à satisfaire:\smallbreak
\noindent
\begin{equation*}
\left\{
    \begin{array}{ll}
        \sum_{k=1}^{N}b_kcos{2\delta_k} = 0\\
        \sum_{k=1}^{N}b_ksin{2\delta_k} = 0\\
        \sum_{k=1}^{N}cos{2\delta_k} = 0\\
        \sum_{k=1}^{N}sin{2\delta_k} = 0
    \end{array}
\right.  
\end{equation*}

\noindent
Soit 4 équations à résoudre donc il faut au minimum 4 ou 5 couches pour satisfaire l'ensemble de ces conditions.
\smallbreak
\noindent
\textit{Comment résoudre le problème pour un nombre de couches plus élevé ?}\smallbreak
\noindent
Pour résoudre un problème par exemple à 6 couches, on peut fixer une orientation de l'une des couches puis résoudre le problème pour 5 couches en prenant une couche parmi les cinq définissant le repère. Mais pour beaucoup plus de couches il faut résoudre un problème d'optimisation.
\smallbreak
\noindent
\textit{
Écrire aussi les conditions à satisfaire pour obtenir le découplage élastique et comparer avec les relations établies pour la stabilité thermo-élastique. Commenter}\smallbreak
\noindent
Pour obtenir le découplage élastique il faut que \textbf{B} = \textbf{0}, cette condition est satisfaite si et seulement si :
\begin{equation*}
\left\{
    \begin{array}{ll}
        \hat{T}_0 = 0\\
        \hat{T}_1= 0\\
        \hat{R}_0e^{4i\hat{\Phi}_0} = 0 \\
        \hat{R}_1e^{2i\hat{\Phi}_0} = 0
    \end{array}
        \right.
  \Leftrightarrow
  \left\{
    \begin{array}{ll}
        T_0^{CB}\big(\sum_{k=1}^{N}b_k\big) = 0 \quad (\text{toujours vrai}) \\
        T_1^{CB}\big(\sum_{k=1}^{N}b_k\big) = 0 \quad  (\text{toujours vrai})\\
        {R}_0^{CB}\big(\sum_{k=1}^{N}b_ke^{4i\delta_k}\big) = 0 \\
        {R}_1^{CB}\big(\sum_{k=1}^{N}b_ke^{2i\delta_k}\big)= 0
    \end{array}
        \right.
     \Leftrightarrow
     \left\{
        \begin{array}{ll}
        \sum_{k=1}^{N}b_kcos{4\delta_k} = 0\\
        \sum_{k=1}^{N}b_ksin{4\delta_k} = 0\\
        \sum_{k=1}^{N}b_kcos{2\delta_k} = 0\\
        \sum_{k=1}^{N}b_ksin{2\delta_k} = 0
    \end{array}
    \right.
\end{equation*}
\noindent
On s'aperçoit alors que deux des conditions sur quatre à satisfaire pour obtenir \textbf{B} = \textbf{0} correspondent aux deux conditions à satisfaire pour obtenir \textbf{V} = \textbf{0}

\smallbreak
\noindent
\textit{Une solution bien connue pour la stabilité thermo-élastique a été obtenue de manière empirique par Winckler (1985). Il s’agit de stratifiés à n = 8 couches caractérisés par la séquence :}
\begin{equation*}
[\theta/(\theta - 90)_2/\theta/-\theta/(90 - \theta)_2/-\theta]
\end{equation*}

\subsection*{Etude de stratifiés de Winckler}
\noindent
\textit{Vérifier que les stratifiés de Winckler sont thermo-élastiquement stables quelle que soit la valeur de l’angle $\theta$.}
\smallbreak
\noindent
Calcul de $b_k$ pour n = 8 couches, avec $b_k = 2k-n-1$ :
\begin{equation*}
    b_1 = -7 ; \quad
    b_2 = -5 ; \quad
    b_3 = -3 ; \quad
    b_4 = -1 ; \quad
    b_5 = 1 ; \quad
    b_6 = 3 ; \quad
    b_7 = 5 ; \quad
    b_8 = 7
\end{equation*}
\noindent
Calcul de $\sum_{k=1}^{N}b_kcos{2\delta_k}$ pour la séquence de \textit{Wincker},

\begin{align*}
\sum_{k=1}^{N}b_kcos{2\delta_k} &= -7cos2\theta - 5cos2(\theta - 90) - 3cos2(\theta - 90) - cos2\theta\\
&+ cos2(-\theta) + 3cos2(90 - \theta) + 5cos2(90 - \theta) + 7cos2(-\theta) 
\end{align*}
or $cos(x) = cos(-x)$
\begin{align*}
\sum_{k=1}^{N}b_kcos{2\delta_k} &= -7cos2\theta - 5cos2(\theta - 90) - 3cos2(\theta - 90) - cos2\theta\\
&+ cos2(\theta) + 3cos2(\theta - 90) + 5cos2(\theta - 90 ) + 7cos2(\theta) = 0
\end{align*}
\noindent
donc, 
\begin{equation*}
\fcolorbox{red}{white}{$\sum_{k=1}^{N}b_kcos{2\delta_k} = 0 \quad \forall \theta$}
\end{equation*}
\smallbreak
\noindent
Calcul de $\sum_{k=1}^{N}b_ksin{2\delta_k}$ pour la séquence de \textit{Winckler},
\begin{align*}
\sum_{k=1}^{N}b_ksin{2\delta_k} &= -7sin2\theta - 5sin2(\theta - 90) - 3sin2(\theta - 90) - sin2\theta\\
&+ sin2(-\theta) + 3sin2(90 - \theta) + 5sin2(90 - \theta) + 7sin2(-\theta)
\end{align*}
or $sin(-x) = -sin(x)$, $sin(x-90) = -cos(x)$ et $sin(2x) = 2sin(x)cos(x)$, donc
\begin{align*}
\sum_{k=1}^{N}b_ksin{2\delta_k} &= -16sin2\theta - 16sin2(\theta - 90)\\
&= -32sin\theta cos\theta + 32sin\theta cos\theta = 0 
\end{align*}
\noindent
donc, 
\begin{equation*}
 \fcolorbox{red}{white}{$\sum_{k=1}^{N}b_ksin{2\delta_k} = 0 \quad \forall \theta$}
\end{equation*}
\smallbreak
\noindent
Ainsi les deux conditions à satisfaire pour avoir \textbf{V} = \textbf{0} sont vérifiées.
\smallbreak
\noindent
De la même façon en calculant  $\sum_{k=1}^{N}cos{2\delta_k}$ et  $\sum_{k=1}^{N}sin{2\delta_k}$ pour la séquence de \textit{Winckler} ont trouve que ces sommes s'annule 
\begin{equation*}
  \fcolorbox{red}{white}{$\sum_{k=1}^{N}cos{2\delta_k} = 0$} \quad ; \quad  \fcolorbox{red}{white}{$\sum_{k=1}^{N}sin{2\delta_k} = 0$}
\end{equation*}
\noindent
et satisfont donc les conditions pour que \textbf{U} soit sphérique.
\smallbreak
\noindent
\textbf{Conclusion}: La séquence de \textit{Winckler} vérifie les conditions pour que le stratifié soit thermo-elastiquement stable (\textbf{V} = \textbf{0}) et pour que les dilatations thermique s'effectue sphériquement dans le plan (\textbf{U} sphérique).

\smallbreak
\noindent
5. \textit{Calculer les propriétés de rigidité en membrane et flexion, \textbf{A} et \textbf{D}, ainsi que le couplage élastique \textbf{B} pour ces stratifiés. Commenter sur la forme et les propriétés de ces comportements.}
\smallbreak
\noindent
On a vu que la séquence de \textit{Winckler} permet d'avoir,
\begin{equation*}
  \sum_{k=1}^{N}e^{2i\delta_k} = 0 \quad \text{et} \quad \sum_{k=1}^{N}b_ke^{2i\delta_k} = 0
\end{equation*}
\noindent
or le terme en $R_1$ de la forme polaire de \textbf{A} et \textbf{B} s'écrivent respectivement comme suit,

\begin{equation*}
   \bar{R}_1e^{2i\bar{\Phi}_1} = \frac{h}{N}R_1^{CB}e^{2i\Phi_1^{CB}}\sum_{k=1}^{N}e^{2i\delta_k} \quad ;\quad \hat{R}_1e^{2i\hat{\Phi}_1} = \frac{1}{2}\frac{h^2}{N^2}R_1^{CB}e^{2i\Phi_1^{CB}}\sum_{k=1}^{N}b_ke^{2i\delta_k} 
\end{equation*}
donc,
\begin{equation*}
    \fcolorbox{red}{white}{$\bar{R}_1=0$} \quad ;\quad  \fcolorbox{red}{white}{$\hat{R}_1 = 0$}
\end{equation*}
\noindent
On en conclut que \textbf{A} et \textbf{B} ont au moins comme propriété la symétrie du carré. \\
Pour s'en assurer on calcul $\sum_{k=1}^{N}cos{4\delta_k}$ pour la séquence de \textit{Wincker} qui nous donne,
\begin{equation*}
\sum_{k=1}^{N}cos{4\delta_k} = 8cos4\theta 
\end{equation*}
\smallbreak
\noindent
, le calcul de $\sum_{k=1}^{N}sin{4\delta_k}$ pour la séquence de \textit{Winckler}, nous donne
\begin{equation*}
\sum_{k=1}^{N}sin{4\delta_k} = 0 \quad \forall \theta
\end{equation*}
\noindent
, calcul de $\sum_{k=1}^{N}b_kcos{4\delta_k}$ pour la séquence de \textit{Winckler}, nous donne
\begin{equation*}
\sum_{k=1}^{N}b_kcos{4\delta_k} = 0 \quad \forall \theta 
\end{equation*}
\smallbreak
\noindent
et calcul de $\sum_{k=1}^{N}b_ksin{4\delta_k}$ pour la séquence de \textit{Winckler}, nous donne
\begin{equation*}
\sum_{k=1}^{N}b_ksin{4\delta_k} = -32sin4\theta
\end{equation*}
\noindent
or le terme en $R_0$ de la forme polaire de \textbf{A} et \textbf{B} s'écrivent respectivement comme suit,
\begin{equation*}
   \bar{R}_0e^{4i\bar{\Phi}_0} = \frac{h}{N}R_0^{CB}e^{4i\Phi_0^{CB}}\sum_{k=1}^{N}e^{4i\delta_k} \quad ;\quad \hat{R}_0e^{4i\hat{\Phi}_0} = \frac{1}{2}\frac{h^2}{N^2}R_0^{CB}e^{4i\Phi_0^{CB}}\sum_{k=1}^{N}b_ke^{4i\delta_k} 
\end{equation*}
donc,
\begin{equation*}
   \fcolorbox{red}{white}{$\bar{R}_0e^{4i\bar{\Phi}_0} = \frac{h}{N}8R_0^{CB}e^{4i\Phi_0^{CB}}cos4\theta$}  \quad ;\quad \fcolorbox{red}{white}{$\hat{R}_0e^{4i\hat{\Phi}_0}=-\frac{1}{2}\frac{h^2}{N^2}32R_0^{CB}e^{4i\Phi_0^{CB}}sin4\theta$}
\end{equation*}
\noindent
On en conclut que \textbf{A} et \textbf{B} sont donc au moins à symétrie du carré. \\

\noindent
Qu'en est - il de D ?
\smallbreak
\noindent
Calcul de $d_k$ pour n = 8 couches, avec $d_k = 12k(k-n-1)+4+3n(n+2)$ :
\begin{equation*}
    d_1 = 148 ; \quad
    d_2 = 76 ; \quad
    d_3 = 28 ; \quad
    d_4 = 4 ; \quad
    d_5 = 28 ; \quad
    d_7 = 76 ; \quad
    d_8 = 148
\end{equation*}
\noindent
Le calcul de $\sum_{k=1}^{N}d_kcos{2\delta_k}$ pour la séquence de \textit{Wincker}, nous donne
\begin{equation*}
\sum_{k=1}^{N}d_kcos{2\delta_k} = 96cos2\theta 
\end{equation*}
\smallbreak
\noindent
Le calcul de $\sum_{k=1}^{N}d_ksin{2\delta_k}$ pour la séquence de \textit{Winckler}, nous donne
\begin{equation*}
\sum_{k=1}^{N}d_ksin{2\delta_k} = 0  \quad \forall \theta
\end{equation*}

\noindent
Le calcul de $\sum_{k=1}^{N}d_kcos{2\delta_k}$ pour la séquence de \textit{Wincker}, nous donne
\begin{equation*}
\sum_{k=1}^{N}d_kcos{4\delta_k} = 512cos4\theta 
\end{equation*}
\smallbreak
\noindent
Le calcul de $\sum_{k=1}^{N}d_ksin{2\delta_k}$ pour la séquence de \textit{Winckler}, nous donne
\begin{equation*}
\sum_{k=1}^{N}d_ksin{4\delta_k} = 0  \quad \forall \theta
\end{equation*}

\noindent
or les termes en $R_0$ et $R_1$ de la forme polaire de \textbf{D} s'écrivent respectivement comme suit,

\begin{equation*}
  \tilde{R}_0e^{4i\tilde{\Phi}_0} = \frac{1}{12}\frac{h^3}{N^3}R_0^{CB}e^{4i\Phi_0^{CB}}\sum_{k=1}^{N}d_ke^{4i\delta_k} \quad ;\quad  
  \tilde{R}_1e^{2i\tilde{\Phi}_1} = \frac{1}{12}\frac{h^3}{N^3}R_0^{CB}e^{2i\Phi_1^{CB}}\sum_{k=1}^{N}d_ke^{2i\delta_k}
\end{equation*}
\noindent
$\Rightarrow$
\begin{equation*}
   \fcolorbox{red}{white}{$\tilde{R}_0e^{4i\tilde{\Phi}_0} = \frac{1}{12}\frac{h^3}{N^3}R_0^{CB}e^{4i\Phi_0^{CB}}512cos4\theta$} \quad ;\quad \fcolorbox{red}{white}{$\tilde{R}_1e^{2i\tilde{\Phi}_1} = \frac{1}{12}\frac{h^3}{N^3}R_1^{CB}e^{2i\Phi_1^{CB}}96cos2\theta$}
\end{equation*}
\noindent
On en conclut que \textbf{D} est de propriété orthotrope ordinaire.\\
\smallbreak
\noindent
6. \textit{Ecrire la forme cartésienne du tenseur  \textbf{B}. Commenter sur la forme de ce tenseur et le type de couplages qu’il représente.}
\smallbreak
\noindent
Pour cela on utilise les expressions de passage polaire-cartésien de la page 4 de ce projet, que l'on particularise au cas où $\hat{T}_0=\hat{T}_1=\hat{R}_1=0$ et $\hat{R}_0 \neq 0$, le tenseur \textbf{B} s'écrit alors:
\begin{equation*}
\textbf{B} =
\begin{pmatrix}
\hat{R}_0cos4\hat{\Phi}_0 & -\hat{R}_0cos4\hat{\Phi}_0 & \hat{R}_0sin4\hat{\Phi}_0  \\
-\hat{R}_0cos4\hat{\Phi}_0  & \hat{R}_0cos4\hat{\Phi}_0 & -\hat{R}_0sin4\hat{\Phi}_0 \\
\hat{R}_0sin4\hat{\Phi}_0  & -\hat{R}_0sin4\hat{\Phi}_0 & -\hat{R}_0cos4\hat{\Phi}_0 
\end{pmatrix}
\end{equation*}
\noindent
Or nous avons vu dans la première partie que pour obtenir un couplage traction-torsion on peut imposer $sin4\hat{\Phi}_0 = \pm 1$ avec ${\Phi}^{CB}_0 = \frac{\pi}{8}$ ce qui implique automatiquement la correspondance $cos4\hat{\Phi}_0 = 0$.
% Dans le cas général de matériau orthotrope $\Phi_1 = 0^o$ dans la couche de base et que $\Phi_1 = 0^o$-$\Phi_1 = k\.rac{\pi}{4}$, on prend ici $k=0$. Soit, $sin4\hat{\Phi}_0 = \pm 1$ avec ce qui implique automatiquement la correspondance $cos4\hat{\Phi}_0 = 0$ 
\begin{equation*}
\textbf{B}=
\begin{bmatrix}
0 & 0 &  \hat{R}_0sin4\hat{\Phi}_0 \\
0 & 0 & - \hat{R}_0sin4\hat{\Phi}_0  \\
\hat{R}_0sin4\hat{\Phi}_0 & - \hat{R}_0sin4\hat{\Phi}_0 & 0
\end{bmatrix}
\end{equation*}
\noindent
ainsi cette composante du tenseur vaut $\hat{R}_0sin4\hat{\Phi}_0  = - \frac{1}{2}\frac{h^2}{N^2}32R_0^{CB}sin4\theta$, fixons maitenant le paramètre $\hat{p}_0$ tel que :
\begin{equation*}
\hat{p}_0=-\frac{1}{2}\frac{h^2}{N^2}32R_0^{CB}
\end{equation*}
il vient,
\begin{equation*}
\textbf{B}=
\begin{bmatrix}
0 & 0 &  \hat{p}_0sin4\theta  \\
0 & 0 & - \hat{p}_0sin4\theta \\
\hat{p}_0sin4\theta& - \hat{p}_0sin4\theta & 0
\end{bmatrix}
\end{equation*}
\noindent
On obtient alors un tenseur de couplage de type \textbf{Extension-Twisting and Shearing-Bending} (ETSB), il permet alors d'avoir un couplage élastique du stratifié en "traction-torsion" et/ou cisaillement-flexion, comme vu en première partie, mais pas seulement...\smallbreak
% ou alors si on veut annuler les entièrement les effets de traction-torsion, on tourne la séquence des plis dans le stratifié de telle manière à obtenir $cos4\hat{\Phi}_0 = \pm 1$ ce qui implique que $sin4\hat{\Phi}_0 = 0$.
% \begin{equation*}
% \textbf{B}=
% \begin{bmatrix}
% \pm sin4\theta\hat{R}_0  & \mp sin4\theta\hat{R}_0  & 0  \\
% \mp sin4\theta\hat{R}_0  & \pm sin4\theta\hat{R}_0  & 0  \\
% 0 & 0 & \mp sin4\theta\hat{R}_0 
% \end{bmatrix}
% \end{equation*}
\noindent
En effet, par rapport au cas traité en première partie, en plus d'avoir un stratifié thermo-elastiquement stable (\textbf{V} = 0 et \textbf{U} sphérique) la séquences de \textit{Winckler} permet de fournir en plus \textbf{la possibilité de moduler l'intensité du couplage} élastique \textbf{B} à l'aide du paramètre géométrique $\theta$ désiré.
\smallbreak
\noindent
7. \textit{Comment varient les propriétés élastiques des stratifiés de Winckler en fonction de l’angle $\theta$? Dire pour quelle valeur de l’angle $\theta$ on obtient le couplage B maximum. Pour cette valeur  $\theta_{max}$ écrire aussi les expressions des tenseurs \textbf{A} et \textbf{D}.Tracer les courbes représentant les composantes A11 et D11 en fonction de l’angle du référentiel.}
\smallbreak
\noindent
Les propriétés élastiques des stratifiés de Winckler en fonction de l’angle $\theta$ varie comme 4 fois l'angle, on parle alors d'harmonique ou de \textbf{periodicité d'ordre 4}. En effet, $sin\theta$ va de $0$ à $1$ entre $[0;\frac{\pi}{2}]$ alors qu'ici pour cette intervalle on a,
\begin{align*}
&\theta = 0 \Rightarrow  sin4\theta = 0\\
&\theta = \frac{\pi}{8} \Rightarrow  sin4\theta = 1\\
&\theta = \frac{\pi}{4} \Rightarrow  sin4\theta = 0\\
&\theta = \frac{3\pi}{8} \Rightarrow  sin4\theta = -1\\
&\theta = \frac{\pi}{2} \Rightarrow  sin4\theta = 0
\end{align*}
\noindent
On en conclut que la valeur de $\theta$ qui donne un couplage \textbf{B} maximum est,
\begin{equation*}
\theta_{max} = \frac{\pi}{8} \qquad [\frac{2\pi}{8}]
\end{equation*}
\noindent
le signe $\pm$ ne change que le sens de la déformation géométriquement.
%(et \textbf{B} est à symétrie du carré, on peut se restreindre à l'intervalle $[0;\frac{\pi}{4}]$)
Écrivons maintenant les expressions des tenseurs \textbf{A} et \textbf{D} pour la valeur  $\theta$ = $\theta_{max}$ : \smallbreak
\noindent
Pour \textbf{A}, nous avons vu que $\bar{R}_1 = 0$, et nous avons aussi $\bar{T}_0 = hT_0^{CB}$ et $\bar{T}_1 = hT_1^{CB}$.

\begin{equation*}
\bar{R}_0e^{4i\bar{\Phi}_0} = \underbrace{\frac{h}{N}8R_0^{CB}}_{\bar{p}_0}cos4\theta(cos4\Phi_0^{CB} + isin4\Phi_0^{CB})    
\end{equation*}
\noindent
le tenseur \textbf{A} s'écrit alors:
\begin{equation*}
\textbf{A} =
\begin{pmatrix}
\bar{T}_0 + 2\bar{T}_1 + \bar{R}_0cos4\bar{\Phi}_0 & -\bar{T}_0 + 2\bar{T}_1 -\bar{R}_0cos4\bar{\Phi}_0 & \bar{R}_0sin4\bar{\Phi}_0    \\
 -\bar{T}_0 + 2\bar{T}_1 - \bar{R}_0cos4\bar{\Phi}_0 & \bar{T}_0 + 2\bar{T}_1 + \bar{R}_0cos4\bar{\Phi}_0  & -\bar{R}_0sin4\bar{\Phi}_0  \\
 \bar{R}_0sin4\bar{\Phi}_0 & -\bar{R}_0sin4\bar{\Phi}_0  & \bar{T}_0 -\bar{R}_0sin4\bar{\Phi}_0 
\end{pmatrix}
\end{equation*}
\noindent
Traçons alors la composante $A_{11}$ pour $\Phi_0^{CB}=0$ : $A_{11} = h(T_0^{CB}+ 2T_1^{CB}) +  \bar{p}_0cos4\theta$ 
\begin{figure}[H]
\hspace{-0.5cm}
    \includegraphics[width=160mm]{./images/A11.png}
    \vspace{-1cm}
    \caption{Représentation matérielle graphique}
\end{figure}
\vspace{-0.5cm}
\noindent
On en conclut que \textbf{$A_{11}$} est dans ce cas à \textbf{symétrie du carré}, étant donné que la rigidité menbranaire est identique dans les deux directions principales. Or, dans le cadre d'un couplage traction-torsion nous avons imposé lors du calcul de B que $\Phi_0^{CB}=\frac{\pi}{8}$.\\
Traçons alors la composante $A_{11}$ dans ce cas $\rightarrow$ $A_{11} = h(T_0^{CB}+ 2T_1^{CB})$
\begin{figure}[H]
\hspace{-0.5cm}
    \includegraphics[width=160mm]{./images/isotropos.png}
    \vspace{-1cm}
    \caption{Représentation matérielle graphique}
\end{figure}
\vspace{-0.5cm}
\noindent
On en conclut que \textbf{$A_{11}$} est dans ce cas à \textbf{isotrope}, étant donné que la rigidité menbranaire est identique dans toutes les directions. \\
Donc dans le cadre du couplage traction-torsion \textbf{A} s'écrit:
\begin{equation*}
\textbf{A} =
\begin{pmatrix}
\bar{T}_0 + 2\bar{T}_1 & -\bar{T}_0 + 2\bar{T}_1 & \bar{R}_0sin4\bar{\Phi}_0    \\
 -\bar{T}_0 + 2\bar{T}_1 & \bar{T}_0 + 2\bar{T}_1  & -\bar{R}_0sin4\bar{\Phi}_0  \\
 \bar{R}_0sin4\bar{\Phi}_0 & -\bar{R}_0sin4\bar{\Phi}_0  & \bar{T}_0
\end{pmatrix}
\end{equation*}
\noindent
soit,
\begin{equation*}
\textbf{A} =
\begin{pmatrix}
h(T_0^{CB}+ 2T_1^{CB}) & h(2T_1^{CB}-T_0^{CB})  &  \bar{p}_0cos4\theta \\
h(2T_1^{CB}-T_0^{CB})  & h(T_0^{CB}+ 2T_1^{CB}) & -\bar{p}_0cos4\theta  \\
\bar{p}_0cos4\theta & -\bar{p}_0cos4\theta  & hT_0^{CB} 
\end{pmatrix}
\end{equation*}
\noindent
\noindent
Et pour $\theta_{max} = \frac{\pi}{8}$ on a $cos4\theta = 0$, donc $\bar{R}_0= 0$,
le tenseur \textbf{A} est alors totalement isotrope et s'écrit,
\begin{equation*}
\textbf{A}=
\begin{pmatrix}
h(T_0^{CB}+ 2T_1^{CB}) &  h(2T_1^{CB}-T_0^{CB})  & 0  \\
 h(2T_1^{CB}-T_0^{CB})  & h(T_0^{CB} + 2T_1^{CB}) & 0 \\
0 & 0 & hT_0^{CB}
\end{pmatrix}
\end{equation*}

\smallbreak
\noindent
Pour \textbf{D}, nous avons $\tilde{T}_0 = \frac{h^3}{3}\frac{T_0^{CB}}{4}$ et $\tilde{T}_1 = \frac{h^3}{3}\frac{T_1^{CB}}{4}$ et $\Phi_0=\frac{\pi}{8}$ et $\Phi_1 = 0$
\begin{equation*}
   \tilde{R}_0e^{4i\tilde{\Phi}_0} = \underbrace{\frac{1}{12}\frac{h^3}{N^3}R_0^{CB}512}_{\tilde{p}_0}e^{4i\Phi_0^{CB}}cos4\theta \quad ;\quad \tilde{R}_1e^{2i\tilde{\Phi}_1} = \underbrace{\frac{1}{12}\frac{h^3}{N^3}R_1^{CB}96}_{\tilde{p}_1}e^{2i\Phi_1^{CB}}cos2\theta  
\end{equation*}

\begin{align*}
&D_{11} = \tilde{T}_0 + 2\tilde{T}_1 + \tilde{R}_0cos4\tilde{\Phi}_0 + 4\tilde{R}_1cos2\tilde{\Phi}_1  \\
&D_{22} = \tilde{T}_0 + 2\tilde{T}_1 + \tilde{R}_0cos4\tilde{\Phi}_0 - 4\tilde{R}_1cos2\tilde{\Phi}_1\\
&D_{12} = -\tilde{T}_0 + 2\tilde{T}_1 - \tilde{R}_0cos4\tilde{\Phi}_0 \\
&D_{16} = \tilde{R}_0sin4\tilde{\Phi}_0 + 4\tilde{R}_1sin2\tilde{\Phi}_1 \\
&D_{26} = -\tilde{R}_0sin4\tilde{\Phi}_0 + 4\tilde{R}_1sin2\tilde{\Phi}_1  \\
&D_{66} =  \tilde{T}_0 - \tilde{R}_0cos4\tilde{\Phi}_0 
\end{align*}

\noindent
Traçons alors la composante $D_{11}$ pour $\Phi_0^{CB}=\Phi_1^{CB}=0$ : $D_{11} = \frac{h^3}{12}(T_0^{CB} + 2T_1^{CB}) +  \tilde{p}_0cos4\theta + 4\tilde{p}_1cos2\theta$ :
\begin{figure}[H]
    \includegraphics[width=160mm]{./images/D11.png}
    \vspace{-.5cm}
    \caption{Représentation matérielle graphique}
\end{figure}
\vspace{-0.5cm}
\noindent
On en conclut que \textbf{$D_{11}$} est bien \textbf{orthotrope ordinaire}, étant donné que la rigidité en flexion est plus grande dans la direction 1.\\
Or, dans le cadre d'un couplage traction-torsion nous avons imposé lors du calcul de B que $\Phi_0^{CB}=\frac{\pi}{8}$.\\
Traçons alors la composante $D_{11}$ dans ce cas $\rightarrow$ $D_{11} = \frac{h^3}{12}(T_0^{CB} + 2T_1^{CB}) + 4\tilde{p}_1cos2\theta$
\begin{figure}[H]
    \includegraphics[width=160mm]{./images/ovale.png}
    \caption{Représentation matérielle graphique}
\end{figure}
\noindent
On en conclut que \textbf{$D_{11}$} est bien \textbf{orthotrope ordinaire}.\\
Donc, dans le cadre du couplage traction-torsion on a:
\begin{equation*}
   \tilde{R}_0e^{4i\tilde{\Phi}_0} = \underbrace{\frac{1}{12}\frac{h^3}{N^3}R_0^{CB}512}_{\tilde{p}_0}cos4\theta \quad ;\quad \tilde{R}_1e^{2i\tilde{\Phi}_1} = \underbrace{\frac{1}{12}\frac{h^3}{N^3}R_1^{CB}196}_{\tilde{p}_1}cos2\theta  
\end{equation*}
et \textbf{D} s'écrit:
\begin{equation*}
\textbf{D} =
\begin{pmatrix}
\tilde{T}_0 + 2\tilde{T}_1 + 4\tilde{R}_1cos2\tilde{\Phi}_1 & -\tilde{T}_0 + 2\tilde{T}_1 & \tilde{R}_0sin4\tilde{\Phi}_0 \\
-\tilde{T}_0 + 2\tilde{T}_1 & \tilde{T}_0 + 2\tilde{T}_1 - 4\tilde{R}_1cos2\tilde{\Phi}_1  & -\tilde{R}_0sin4\tilde{\Phi}_0 \\
\tilde{R}_0sin4\tilde{\Phi}_0  & -\tilde{R}_0sin4\tilde{\Phi}_0 & \tilde{T}_0
\end{pmatrix}
\end{equation*}
\noindent
soit, 
\begin{equation*}
\textbf{D} =
\begin{pmatrix}
\frac{h^3}{12}(T_0^{CB} + 2T_1^{CB}) + 4\tilde{p}_1cos2\theta  & -\frac{h^3}{12}(T_0^{CB} - 2T_1^{CB}) & \tilde{p}_0cos4\theta \\
-\frac{h^3}{12}(T_0^{CB} - 2T_1^{CB})  & \frac{h^3}{12}(T_0^{CB} + 2T_1^{CB}) - 4\tilde{p}_1cos2\theta  & - \tilde{p}_0cos4\theta\\
\tilde{p}_0cos4\theta & - \tilde{p}_0cos4\theta & \frac{h^3}{12}T_0^{CB} 
\end{pmatrix}
\end{equation*}
\noindent
Pour $\theta_{max} = \frac{\pi}{8}$ on a $cos4\theta = 0$,
le tenseur \textbf{D} est de propriété d'horthotropie $R_0$ (car $R_0 = 0$), et s'écrit :\\
\textbf{D} =
\begin{equation*}
\begin{pmatrix}
\frac{h^3}{12}(T_0^{CB} + 2T_1^{CB}) + 4\tilde{p_1}\frac{\sqrt{2}}{2} & -\frac{h^3}{12}(T_0^{CB} - 2T_1^{CB}) & 0  \\
-\frac{h^3}{12}(T_0^{CB} - 2T_1^{CB}) & \frac{h^3}{12}(T_0^{CB} + 2T_1^{CB}) + 4\tilde{p_1}\frac{\sqrt{2}}{2} & 0 \\
0 & 0 & \frac{h^3}{12}T_0^{CB}  
\end{pmatrix}
\end{equation*}
\smallbreak
\noindent
 Remarque : pour tracer les courbes à l'aide de \textbf{Matlab} (cf. \textit{Application.m}) nous avons considéré un stratifié composé de huit couche unidirectionnelles en carbone/epoxyde d'épaisseur h, et modules polaires donnée dans le tableau ci-dessous:\\
\noindent
 \begin{tabular}{|c|c|c|c|c|c|c|}
  \hline
  Modules polaires & $T_0$ (GPa)  & $T_1$ (GPa) & $R_0$ (GPa) & $R_1$ (GPa) & $\Phi_0,\Phi_1=0$ & h (mm) \\
  \hline
  & 26,88 & 24,74 & 19,71 & 21,43 & 0 & 1 \\
  \hline
\end{tabular}
\smallbreak
\noindent
8. \textit{En réalité, pour comprendre la réponse d’un stratifié à un chargement mécanique, il faut étudier les tenseurs de souplesse : inverser la loi de comportement thermo-élastique, puis particulariser les équations dans le cas d’un chargement purement mécanique.}

\smallbreak
\noindent
Je ne veux pas mettre ici toutes les étapes du calcul. Mais la démarche est la suivante on part de, 

\begin{equation*}\left\{
    \begin{array}{ll}
        \textbf{N} = \textbf{A}\boldsymbol{\mathcal{E}}^0 + \textbf{B}\boldsymbol{\mathcal{K}} - T_0\textbf{U} \\
        \textbf{M} = \textbf{B}\boldsymbol{\mathcal{E}}^0 + \textbf{D}\boldsymbol{\mathcal{K}}
    \end{array}
\right.
\end{equation*}
\noindent
en prenant la première équation, on isole $\mathcal{E}^0$ en pré-multipliant par $\textbf{A}^{-1}$. On obtient alors:

\begin{equation*}
\boldsymbol{\mathcal{K}} = \textbf{d}\textbf{M} + \textbf{b}_2(\textbf{N}-T_0\textbf{U})
\end{equation*}
\noindent
avec,
\begin{align*}
&\textbf{d} = (\textbf{D}-\textbf{B}\textbf{A}^{-1}\textbf{B})^{-1}\\
&\textbf{b}_2 = -\textbf{d}\textbf{B}\textbf{A}^{-1}
\end{align*}
\noindent
de même pour la deuxième équation, on a
\begin{equation*}
\boldsymbol{\mathcal{E}}^0 = \textbf{a}(\textbf{N}-T_0\textbf{U}) + \textbf{b}_1\textbf{M} 
\end{equation*}
\noindent
avec,
\begin{align*}
&\textbf{a} = (\textbf{A}-\textbf{B}\textbf{D}^{-1}\textbf{B})^{-1}\\
&\textbf{b}_1 = -\textbf{a}\textbf{B}\textbf{D}^{-1}
\end{align*}
On regroupe alors sous forme matricielle et on peut aussi montrer que $\textbf{b}_1 = \textbf{b}^T_2$ et posé $\textbf{b}_1 = \textbf{b}$

\begin{equation*}
\begin{bmatrix}
\boldsymbol{\mathcal{E}}^0 \\
\boldsymbol{\mathcal{K}} 
\end{bmatrix} =
\begin{bmatrix}
\textbf{a} & \textbf{b} \\
\textbf{b}^T & \textbf{d} 
\end{bmatrix}
\begin{bmatrix}
\textbf{N}-T_0\textbf{U}\\
\textbf{M}
\end{bmatrix}
\end{equation*}

\begin{equation*}
\Leftrightarrow \; \; \left\{
    \begin{array}{ll}
        \boldsymbol{\mathcal{E}}^0 = \textbf{a}(\textbf{N}-T_0\textbf{U}) + \textbf{b}\textbf{M} \\
        \boldsymbol{\mathcal{K}} = \textbf{d}\textbf{M} + \textbf{b}^T(\textbf{N}-T_0\textbf{U})
    \end{array}
\right.
\end{equation*}
Particularisons à un chargement purement mécanique,

\begin{equation*}
\Rightarrow \; \; \left\{
    \begin{array}{ll}
        \boldsymbol{\mathcal{E}}^0 = \textbf{a}\textbf{N} + \textbf{b}\textbf{M} \\
        \boldsymbol{\mathcal{K}} = \textbf{d}\textbf{M} + \textbf{b}^T\textbf{N}
    \end{array}
\right.
\end{equation*}
\noindent
% Nous analysons le terme non nul $\textbf{b}\textbf{M}$ qui couple les effets d'effort de membrane et le moment.
% \begin{equation*}
% \boldsymbol{\mathcal{E}}^0_{coupl} = \textbf{b}\textbf{M} =
% \begin{bmatrix}
% \mathcal{E}_{1}\\
% \mathcal{E}_{2}\\
% \mathcal{E}_{6}
% \end{bmatrix}=
% \begin{bmatrix}
% b_{11}M_1 + b_{12}M_2 + b_{16}M_6\\
% b_{12}M_1 + b_{22}M_2 + b_{26}M_6\\
% b_{16}M_1 + b_{26}M_2 + b_{66}M_6
% \end{bmatrix}
% \end{equation*}

% \noindent
% Si on veut obtenir un couplage en "traction-torsion" il faut alors activer les composantes $b_{16}$ et $b_{26}$ qui lie justement le moment de torsion $M_6$ avec les déformations de membranes longitudinal $\mathcal{E}^0_1$ et $\mathcal{E}^0_2$, soit
% Pour cela on utilise les expressions de passage polaire-cartésien de la page 4 de ce projet, que l'on particularise au cas où $\hat{T}_0=\hat{T}_1=\hat{R}_1=0$ et $\hat{R}_0 \neq 0$, le tenseur \textbf{B} s'écrit alors:
% \begin{equation*}
% \textbf{B} =
% \begin{pmatrix}
% \hat{R}_0cos4\hat{\Phi}_0 & -\hat{R}_0cos4\hat{\Phi}_0 & \hat{R}_0sin4\hat{\Phi}_0  \\
% -\hat{R}_0cos4\hat{\Phi}_0  & \hat{R}_0cos4\hat{\Phi}_0 & -\hat{R}_0sin4\hat{\Phi}_0 \\
% \hat{R}_0sin4\hat{\Phi}_0  & -\hat{R}_0sin4\hat{\Phi}_0 & -\hat{R}_0cos4\hat{\Phi}_0 
% \end{pmatrix}
% \end{equation*}
% \noindent
% Or nous avons vu à la question précédente que le terme qui ne s'annule pas dans  \textbf{B} est la partie réel en sinus de la composante polaire $\hat{R}_0e^{4i\hat{\Phi}_0}$.
% % Dans le cas général de matériau orthotrope $\Phi_1 = 0^o$ dans la couche de base et que $\Phi_1 = 0^o$-$\Phi_1 = k\.rac{\pi}{4}$, on prend ici $k=0$. Soit, $sin4\hat{\Phi}_0 = \pm 1$ avec ce qui implique automatiquement la correspondance $cos4\hat{\Phi}_0 = 0$ 
% \begin{equation*}
% \textbf{B}=
% \begin{bmatrix}
% 0 & 0 & \pm \hat{R}_0sin4\hat{\Phi}_0 \\
% 0 & 0 & \mp \hat{R}_0sin4\hat{\Phi}_0  \\
% \pm \hat{R}_0sin4\hat{\Phi}_0 & \hat{R}_0sin4\hat{\Phi}_0 & 0
% \end{bmatrix}
% \end{equation*}
% \noindent

% % \begin{figure}[H]
% %     \centering
% %     \includegraphics[width=40mm]{./images/brouillon.jpg}
% %     \caption{Brouillon du calcul}
% % \end{figure}

\smallbreak
\noindent
9. \textit{Calculer les comportements de souplesse pour les stratifiés de type Winckler en fonction de l’angle $\theta$. Quelle est la forme de ces comportements ?
Comment répond le stratifié si sollicité en traction pure selon l’axe x ? Est-ce que le maximum de cette réponse est obtenu pour la même valeur de $\theta = \theta_{max}$ ?}
\smallbreak
\noindent
Pour calculer le comportements de souplesse en fonction de l'angle $\theta$ on utilise \textbf{Matlab} pour pouvoir faire les calculs de a,b et déterminé à la question précédentes. En effet les calculs de a, b et d en fonction de $\theta$ sont grand, nous donnons alors ici leurs formes (voir le fichier joint \textbf{\textit{Analytique.m}} pour les expressions explicite de ces tenseurs);

% \begin{framed}
% \FSource{./code/inversion.m}
% \end{framed}
\noindent
 

\begin{equation*}a,d = \begin{bmatrix}
$\textbullet$ & $\textbullet$ &  0\\
$\textbullet$  & $\textbullet$ & 0\\
0 & 0 & $\textbullet$
\end{bmatrix} \quad et \quad b = \begin{bmatrix}
0 & 0 & $\textbullet$\\
0 & 0 & $\textbullet$\\
$\textbullet$ & $\textbullet$ & 0
\end{bmatrix}
\end{equation*}
\noindent
avec,
\begin{align*}
&a_{11} = a_{22} \quad et \quad a_{12}=a_{21}.\\ 
&d_{12}=d_{21}\\
&b_{16} = b_{26}  \quad et \quad  b_{61} \neq b_{62}
\end{align*}
\noindent
Étudions alors la réponse à une traction simple:

\begin{framed}
\FSource{./code/abd.m}
\end{framed}

\begin{equation*}
\begin{bmatrix}
\mathcal{E}_1 \\
\mathcal{E}_2 \\
\mathcal{E}_6 \\
\mathcal{K}_1 \\
\mathcal{K}_2 \\
\mathcal{K}_6 
\end{bmatrix} = \begin{bmatrix}
a_{11} & a_{12} & 0       & 0 & 0           & b_{16}\\
a_{12} & a_{11} & 0       & 0 & 0           & b_{26}\\
     0 & 0      & a_{66}  & b_{61} & b_{62} & 0\\
     0 & 0      & b_{61}  & d_{11} & d_{12} & 0\\
     0 & 0      & b_{62}  & d_{12} & d_{22} & 0 \\
b_{16} & b_{26} & 0       & 0      & 0      & d_{66}
\end{bmatrix}\begin{bmatrix}
\mathcal{N}_x \\
0 \\
0 \\
0 \\
0 \\
0 
\end{bmatrix}
\end{equation*}
\noindent
$\Rightarrow$
\begin{equation*}
\begin{bmatrix}
\mathcal{E}_1 \\
\mathcal{E}_2 \\
\mathcal{E}_6 \\
\mathcal{K}_1 \\
\mathcal{K}_2 \\
\mathcal{K}_6 
\end{bmatrix} =\begin{bmatrix}
a_{11}\mathcal{N}_x \\
a_{12}\mathcal{N}_x \\
0\\
0 \\
0 \\
a_{16}\mathcal{N}_x  
\end{bmatrix}
\end{equation*}
\noindent
On obtient pour une traction uniaxial des effets de déformations menbranaire, mais aussi un effet de courbure, due au terme $b_{16}$.
Ce couplage varie en fonction de $\theta$ et puisque cette expression comporte des termes en $cos4\theta$ et $sin4\theta$ on peut imaginer une valeur de $\theta_{max}$ proche de $\frac{\pi}{16}$, en effet d'après la courbe de suivante de $b_{16}$ en fonction de $\theta$:
\begin{figure}[H]
    \includegraphics[width=160mm]{./images/max.png}
    \vspace{-.5cm}
    \caption{Variation de la composante de couplage $b_{16}$}
\end{figure}
\noindent
on vérifie que $\theta_{max} = \frac{\pi}{18}$ (cf. \textit{Application.m})

\noindent
10. \textit{On considère une pale de largeur a et longueur L fabriquée en composite stratifié selon la séquence (1) avec $\theta = \theta_{max}$. La pale est disposée dans le plan $O_{xy}$ et en rotation autour d’un axe vertical $O_z$ à une vitesse angulaire constante $\omega$ : soit L $>>$ a, de telle manière qu'on puisse confondre la direction axiale et celle radiale. Dans ce cas, on peut dire donc que chaque point de la pale est soumis à un effort de traction $N_x = \rho h a \omega 2x$ (x dans [0 , L]). Calculer la réponse de la pale à la sollicitation centrifuge (seule la composante $N_x$ est non nulle) et en particulier l’angle de torsion en bout de pale. Confirmer ce résultat par un calcul éléments finis sous Abaqus.}

\noindent
Pour effectuer la modélisation d'une pale en traction-compression nous avons considéré une pale d'une lonqueur de $6m$, de l'argeur $50cm$ et d'une épaisseur de $8mm$, de masse volumique $1600 Kg/m^3$ et d'une vitesse de rotation de $157rad/s$. Cette pale est constitué d'un stratifié composé de huit couches unidirectionnelles en carbone/epoxyde d'épaisseur h. Détaillons les étapes de la modélisation (cf. \textit{Winckler.cae}):\\
\begin{itemize}
    \item \texttt{Part : 3D Shell Planar}
    \item \texttt{Geometrie : }
    \smallbreak
    \includegraphics[width=80mm]{./images/geometrie.png}
    \smallbreak
    \item \texttt{Property : Elastic - Type : Lamina}\smallbreak
    \begin{tabular}{|c|c|c|c|c|}
  \hline
  Modules élastique & $E_1$ (GPa)  & $E_2$ (GPa) & $G_{12}$ (GPa) & $\nu_{12}$  \\
  \hline
    & 181 &  10,3 & 7,17 & 0,28 \\
  \hline
\end{tabular}
  \smallbreak
 \item \texttt{Composit Layup : Conventional Shell}
   \smallbreak
  \includegraphics[width=140mm]{./images/layup.png}
   \item \texttt{Step : Static, General}
   \item \texttt{Load : Body force, Analytical-Field:\\
   Nx = X*94651.6, Component 1}
     \item \texttt{Boundary condition : ENCASTRE }
     
    \includegraphics[width=100mm]{./images/condition.png}
    
    \item \texttt{Mesh : Family Shell, S4R : finite menbrane strains}
    
    \includegraphics[width=100mm]{./images/mesh.png}
      
    \item \texttt{Job : UR1}
    
    \includegraphics[width=100mm]{./images/job.png}
    \includegraphics[width=20mm]{./images/jauge.png}
    
    \item \texttt{Path dans la lonqueur de la pale:}\\
    \includegraphics[width=120mm]{./images/path1.png}
    
    
    \item \texttt{Path dans l'épaisseur en bout de pale :}\\
    \center{\includegraphics[width=70mm]{./images/RX2.png}}\\
    \flushleft
    \includegraphics[width=120mm]{./images/path2.png}
    
\end{itemize}
\noindent
L'angle de tortion $\alpha$ obtenue en bout de pale est de \textbf{13.94 degré}.\\
\noindent
Théoriquement, pour des plaques minces la théorie cinématique de Love-Kirchoff permet d'écrire: 
\begin{equation*}
\mathcal{K}_6 = 2\frac{\partial^2 W_0(x,y)}{\partial x \partial y}
\end{equation*}
$\Rightarrow$
\noindent
\begin{equation*}
    W_0(x,y) = \frac{1}{2}\int_0^L\int_0^l \mathcal{K}_6 dxdy = \frac{1}{2}K_6*L*l 
\end{equation*}
Donc le déplacement selon $z$  vaut $W_0 = 0.8421 \, m$ et le coté adjacent vaut $adj = \frac{l}{2} = 0.25\,m$ et donc $\alpha=atan(\frac{W_0}{adj}) = 73.47\, \text{degré}$.\\
Cette valeur ne correspond pas à la valeur trouvé numériquement, elle est 5 fois trop grandes, et je ne trouve pas où mon raisonnement est erroné (cf \textit{Application.m})

