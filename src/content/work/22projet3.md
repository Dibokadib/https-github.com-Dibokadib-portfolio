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

### I. Description du problème

<p style="text-align: justify;">
Dans la conception de structures composites il peut être utile de réaliser des stratifiés couplés élastiquement afin de produire des changements de forme lors d’un chargement mécanique : par exemple, le $\textbf{couplage traction-torsion}$ peut être exploité pour construire des aubes ou pales à pas adaptatif, la torsion intervenant lors de la rotation des pales/aubes par effet de la sollicitation produite par effet centrifuge (plus la vitesse de rotation est élevée et plus important sera l’angle de torsion).
</P>

<div style="text-align: center;">
  <figure style="display: inline-block;">
    <img src="/assets/projet3/p3f1.png" alt="Couplage traction/torsion" width="400"/>
    <figcaption>Figure 1 : Couplage traction/torsion</figcaption>
  </figure>
</div>

<p style="text-align: justify;">
\(   \) 
</P>

<p style="text-align: justify;">
Néanmoins, lors de la fabrication, les pièces composites subissent des variations de température et les couplages peuvent produire des déformations/courbures d’origine thermo-élastique, ce qui n’est pas souhaitable pour garantir la stabilité de la géométrie de ces pièces.
</P>

<p style="text-align: justify;">
Il est donc nécessaire de concevoir des stratifiés ayant un couplage thermo-élastique nul $\textbf{V = O}$. On parle dans ce cas de $\textit{stratifiés thermo-élastiquement stables.}$
</P>

<p style="text-align: justify;">
On se propose d’étudier ici comment concevoir des stratifiés de ce type et d’analyser leurs propriétés. Ensuite on considérera une famille particulière de stratifiés thermo-élastiquement stables et on en étudiera le comportement.
</P>

<p style="text-align: justify;">
On considère ici le cas de stratifiés constitués de couches identiques (même matériau et même épaisseur) : les paramètres de conception seront donc le nombre de couches n et la séquence d’angles d’orientation $[\delta k]$ $(k = 1,..., n)$.
</P>

<p style="text-align: justify;">
On utilisera la représentation polaire pour analyser le comportement de stratifiés, partout où il est possible. On notera « CB » les paramètres qui décrivent les propriétés de la couche de base et par des symboles $\,$ $\bar{}$ $\,$, $\,$ $\hat{}$ $\,$ et $\,$ $\tilde{}$ $\,$ les termes en lien respectivement à $\textbf{A}$,  $\textbf{B}$ et  $\textbf{D}$.
</P>

### II. Propriétés générales

<p style="text-align: justify;">
1. Expliquer pourquoi si $V = O$ pour un stratifié à couches identiques, alors le tenseur de couplage élastique B est à symétrie du carré.
</P>

<p style="text-align: justify;">
L'objecif est de montrer que le tenseur $\textbf{B}$ d'ordre 4 de couplage élastique est à symétrie du carré en considérant le cas de figure de stratifié à couches identiques et d'un couplage thermo-élastique nul i.e $\textbf{V=O}$.
</P>

<p style="text-align: justify;">
$\underline{Remarque}$ : le tenseur $\textbf{V}$ permet de coupler les effets de contributions de courbures due à une élévation uniforme de température et les effets de contributions menbranaires due à un gradient thermique dans l'épaisseur du stratifié.
</P>

<p style="text-align: justify;">
Rappelons d'abord la signification d'un stratifié à couches identiques, il correspond à des plis constitués du même matériau et de même épaisseur.
</P>

On part alors de l'expression de $\textbf{V}$:

\begin{equation*}
    \textbf{V}=\frac{1}{2}\sum_{k=1}^{N}{\boldsymbol{\gamma}}(\delta_k)(z^2_k-z^2_{k-1})
\end{equation*}

<p style="text-align: justify;">
avec, $\boldsymbol{\gamma}(\delta_k)$ le tenseur de contrainte thermique par unité de température $(Mpa.^oC^{-1})$,
que l'ont particularise pour des couches identiques :
</P>

\begin{equation*}
    \textbf{V}=\frac{1}{2}\frac{h^2}{N^2}\sum_{k=1}^{N}b_k{\boldsymbol{\gamma}}(\delta_k)
\end{equation*}

<p style="text-align: justify;">
avec le coefficient $b_k = 2k-N-1$
N et h nombre et épaisseur de couche total respectivement.
</P>

<p style="text-align: justify;">
Nous pouvons alors exprimer ce tenseur sous forme cartesienne ou polaire. Pour des raisons de simplifications on écrit le tenseur $\textbf{V}$ sous sa forme polaire à l'aide de la représentation des cercles de Mohr (deux paramètres invariants $\textbf{T}$, $\textbf{R}$ partie sphérique et déviatorique et un angle $\Phi$ de direction principal de contrainte qui varie avec le repère). 
</P>

<p style="text-align: justify;">
La décomposition polaire la plus générale de $\textbf{V}$ nous donne:
</P>

\begin{equation*}
    T_v=\frac{1}{2}\big(\sum_{k=1}^{N}(z^2_k-z^2_{k-1})T_{\gamma}\big)
\end{equation*}

\begin{equation*}
R_v e^{2i\Phi_v} = \frac{1}{2}\big(\sum_{k=1}^{N}(z^2_k-z^2_{k-1})R_{\gamma}e^{2i(\Phi_{\gamma}+\delta_k)}\big)
\end{equation*}

<p style="text-align: justify;">
On particularise pour des couches identiques :
$T_{\gamma}$, $R_{\gamma}$ et $e^{2i\Phi_{\gamma}}$ sont les mêmes pour toutes les couches (constantes).
</P>

$$($z^2_k-z^2_{k-1}$) \rightarrow \frac{h^2}{N^2}b_k $$ 

et 

$$b_{ii} = 0$$

Soit,

\begin{equation*}
    T_v= \frac{1}{2}\frac{h^2}{N^2}T_{\gamma}\big(\sum_{k=1}^{N}b_k\big) = 0
\end{equation*}

\begin{equation*}
 	R_v e^{2i\Phi_v} = \frac{1}{2}\frac{h^2}{N^2}R_{\gamma}e^{2i\Phi_{\gamma}}\big(\sum_{k=1}^{N}b_k e^{2i\delta_k}\big)
\end{equation*}

<p style="text-align: justify;">
Ainsi, un stratifié à couche identique est thermo-élastiquement stable, i.e $\textbf{V = 0}$, si et seulement si 
\begin{equation}
    R_v e^{2i\Phi_v} = 0 \;\;\; \Leftrightarrow \;\;\; \fcolorbox{red}{ }{$\sum_{k=1}^{N}b_ke^{2i\delta_k} = 0$} \;\; ou \;\; \fcolorbox{red}{ }{$R_{\gamma} = 0$} 
    \label{condition}
\end{equation}
</P>

<p style="text-align: justify;">
$\underline{Conclusion}$ : pour avoir $\textbf{V=0}$ dans un stratifié à couche identique il faut satisfaire au moins l'une des deux conditions [\ref{condition}] ci-dessus, soit:
</P><br>
 
> <p style="text-align: justify;"> 1) La première équation est une condition géométrique, elle dépend de la combinaison des angles (séquences) donné au plis et du coefficient $b_k$. Cette condition complexe donne en réalité deux conditions, une venant de la partie réelle et une autre de la partie imaginaire.
</P>
<br>

> <p style="text-align: justify;">  2) La deuxième est une condition matérielle, $R_{\gamma} = 0$ signifie que $\boldsymbol{\gamma}$ est purement sphérique. Pour obtenir cela, il faut supposer que les propriétés du matériau de la couche de base est à symétrie du carré (ou cas trivial i.e isotrope), ainsi $R_1^{CB} = 0$ et impliquant $\hat{R}_1 = \tilde{R}_1 = \bar{R}_1 = 0$, ce qui nous ramène à la propriété de symétrie du carré pour $\textbf{A}$,$\textbf{B}$ et $\textbf{D}$.
</P>
<br>

<p style="text-align: justify;">
Considérons dans la suite la première condition afin de satisfaire la propriété de symétrie du carré seulement pour $\textbf{B}$.
</P>

<p style="text-align: justify;">
Le tenseur $\textbf{B}$ d'ordre 4 possède par ailleurs 4 composante polaire. On effectue le même déroulement que précédemment.
</P>

<p style="text-align: justify;">
On part alors de l'expression de $\textbf{B}$:
</P>

<p style="text-align: justify;">
\begin{equation*}
    \textbf{B}=\frac{1}{2}\sum_{k=1}^{N}\textbf{Q}(\delta_k)(z^2_k-z^2_{k-1})
\end{equation*}
</P>

<p style="text-align: justify;">
avec $\textbf{Q}(\delta_k)$ le tenseur d'ordre $4$ d'élasticité en $2D$ $(Mpa)$, que l'ont particularise pour des couches identiques :
</P>

<p style="text-align: justify;">
\begin{equation*}
    \textbf{B}=\frac{1}{2}\frac{h^2}{N^2}\sum_{k=1}^{N}b_k\textbf{Q}(\delta_k)
\end{equation*}
</P>

<p style="text-align: justify;">
La décomposition polaire de $\textbf{B}$ est donné par $\hat{T}_0$, $\hat{T}_1$, $\hat{R}_0e^{4i\hat{\Phi}_0}$ et $\hat{R}_1e^{2i\hat{\Phi}_1}$ que l'ont écrit sous forme particularisé pour des couches identiques. 
</P>

<p style="text-align: justify;">
Comme pour le calcul de $T_v$ on a $\hat{T}_0$ = $\hat{T}_1 = 0$ car $\sum_{k=1}^{N}b_k = 0$.
</P>

<p style="text-align: justify;">
Pour $\hat{R}_0e^{4i\hat{\Phi}_0}= \frac{1}{2}\big(\sum_{k=1}^{N}R_{0k}e^{4i(\Phi_{0k}+\delta_k)}(z^2_k-z^2_{k-1})\big)$  on a,
</P>

<p style="text-align: justify;">
$$
\hat{R}_0e^{4i\hat{\Phi}_0} =
$$
</P>

<p style="text-align: justify;">
$$ 
	\frac{1}{2}\frac{h^2}{N^2}\big(\sum_{k=1}^{N}b_k{R}_0^{CB}e^{4i(\Phi_0^{CB}+\delta_k)}\big)
$$
</P>

<p style="text-align: justify;">
$$
   =\frac{1}{2}\frac{h^2}{N^2}{R}_0^{CB}e^{4i\Phi_0^{CB}}\big(\sum_{k=1}^{N}b_ke^{4i\delta_k}\big)
$$
</P>

<p style="text-align: justify;">
et pour $\hat{R}_1e^{2i\hat{\Phi}_1}= \frac{1}{2}\big(\sum_{k=1}^{N}R_{1k}e^{2i(\Phi_{1k}+\delta_k)}(z^2_k-z^2_{k-1})\big)$  on a,
</P>

<p style="text-align: justify;">
$$
   \hat{R}_1e^{2i\hat{\Phi}_1} = \frac{1}{2}\frac{h^2}{N^2}\big(\sum_{k=1}^{N}b_k{R}_1^{CB}e^{2i(\Phi_1^{CB}+\delta_k)}\big) 
$$
</P>

<p style="text-align: justify;">
$$
= \frac{1}{2}\frac{h^2}{N^2}{R}_1^{CB}e^{2i\Phi_1^{CB}}\big(\sum_{k=1}^{N}b_ke^{2i\delta_k}\big) 
$$
</P>

<p style="text-align: justify;">
Or d'après la condition [\ref{condition}] pour que $\textbf{V} = 0$, il vient:
</P>

<p style="text-align: justify;">
$$
     \hat{R}_1e^{2i\hat{\Phi}_1}=0 \Leftrightarrow  \fcolorbox{red}{ }{$\hat{R}_1=0$}
$$
</P>

<p style="text-align: justify;">
Le paramètre polaire $\hat{R}_1$ étant nul, nous avons bien la symétrie du carré seulement pour le tenseur $\textbf{B}$, en donnant une orientation adéquate de la séquences d'angles. Cela est satisfait par exemple dans la configuration de bi-couche à $0^o$ et $90^o$ (cross-ply).
</P>

<p style="text-align: justify;">
$\textit{Quelles sont donc les composantes non nulles de $\textbf{B}$ ?}$
</P>

<p style="text-align: justify;">
D'après la question précédente, la composante polaire non nulle de $\textbf{B}$ est
</P>

<p style="text-align: justify;">
$$
\hat{R}_0e^{4i\hat{\Phi}_0}   
$$
</P>

donc,

<p style="text-align: justify;">
\begin{equation*}
\fcolorbox{red}{ }{$\hat{R}_0cos{4\hat{\Phi}_0}$} \;\; et \;\; \fcolorbox{red}{ }{$\hat{R}_0sin{4\hat{\Phi}_0}$}  
\end{equation*}
</P>

<p style="text-align: justify;">
Ecrire dans ce cas la représentation cartésienne de $\textbf{B}$ et analyser les divers types de couplages élastiques qui peuvent se manifester.
</P>

<p style="text-align: justify;">
Dans le cas générique les composantes polaire de $\textit{Verchery}$ pour un tenseur d'ordre $4$ s'écrivent:
</P>

<p style="text-align: justify;">
\begin{align}
&L_{1111} = \quad T_0 + 2T_1 + R_0cos4\Phi_0 + 4R_1cos2\Phi_1\\
&L_{1112} = \quad \quad  \quad  \quad \quad \quad R_0sin4\Phi_0 + 2R_1sin2\Phi_1\\
&L_{1122} =-T_0 + 2T_1 - R_0cos4\Phi_0\\
&L_{1212} = \quad T_0  \quad \quad \quad - R_0cos4\Phi_0\\
&L_{2212} = \quad \quad \quad  \quad \quad \;- R_0sin4\Phi_0 + 2R_1sin2\Phi_1\\
&L_{2222} = \quad T_0 + 2T_1 + R_0cos4\Phi_0 - 4R_1cos2\Phi_1\\
\end{align}
</p>

<p style="text-align: justify;">
et dans notre cas où $\hat{T}_0=\hat{T}_1=\hat{R}_1=0$ le tenseur $\textbf{B}$ s'écrit alors:
</p>


<p style="text-align: justify;">
\begin{equation}
\textbf{B} =
\begin{pmatrix}
\hat{R}_0cos4\hat{\Phi}_0 & -\hat{R}_0cos4\hat{\Phi}_0 & \hat{R}_0sin4\hat{\Phi}_0  \\
-\hat{R}_0cos4\hat{\Phi}_0  & \hat{R}_0cos4\hat{\Phi}_0 & -\hat{R}_0sin4\hat{\Phi}_0 \\
\hat{R}_0sin4\hat{\Phi}_0  & -\hat{R}_0sin4\hat{\Phi}_0 & -\hat{R}_0cos4\hat{\Phi}_0 
\end{pmatrix}
\end{equation}
<p style="text-align: justify;">

<p style="text-align: justify;">
Pour analyser les divers couplages élastiques, écrivons la loi de comportement thermo-élastique (sans gradient thermique et $\textbf{V}=0$):
</P>

<p style="text-align: justify;">
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
</P>

<p style="text-align: justify;">
$$
\Rightarrow \; \; \left\{
    \begin{array}{ll}
        \textbf{N} = \textbf{A}\boldsymbol{\mathcal{E}}^0 + \textbf{B}\boldsymbol{\mathcal{K}} - T_0\textbf{U} \\
        \textbf{M} = \textbf{B}\boldsymbol{\mathcal{E}}^0 + \textbf{D}\boldsymbol{\mathcal{K}}
    \end{array}
\right.
$$
</P>

<p style="text-align: justify;">
Nous analysons le terme non nul $\textbf{B}\boldsymbol{\mathcal{E}}^0$ qui couple les effets d'efforts de membranes et de moments.
</P>

<p style="text-align: justify;">
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
</P>

<p style="text-align: justify;">
Si on veut obtenir un couplage en "traction-torsion" il faut activer les composantes $B_{16}$ et $B_{26}$ qui lie justement le moment de torsion $M_6$ avec les déformations de membranes longitudinal $\mathcal{E}^0_1$ et $\mathcal{E}^0_2$ (en réalité il faut d'abord inverser la loi pour parler d'un tel couplage cf. question 8).
</P>

<p style="text-align: justify;">
A titre d'exemple, les constructeurs de pales d'hélicoptères exploite cette effet de traction-torsion. En effet, lors de la rotation des pales la force centrifuge exercé de manière radial dans le plan des pales va déclencher une torsion de ces pales, cette effet est d'un grand intérêt aérodynamique.
</P>

<p style="text-align: justify;">
Pour cela on peut imposer $sin4\hat{\Phi}_0 = \pm 1$ avec ${\Phi}_0 = \frac{\pi}{8}$ ce qui implique automatiquement la correspondance $cos4\hat{\Phi}_0 = 0$ 
% Dans le cas général de matériau orthotrope $\Phi_1 = 0^o$ dans la couche de base et que $\Phi_1 = 0^o$-$\Phi_1 = k\frac{\pi}{4}$, on prend ici $k=0$. Soit, $sin4\hat{\Phi}_0 = \pm 1$ avec ce qui implique automatiquement la correspondance $cos4\hat{\Phi}_0 = 0$
</P>

<p style="text-align: justify;">
\begin{equation*}
\textbf{B}=
\begin{bmatrix}
0 & 0 & \pm\hat{R}_0  \\
0 & 0 & \mp\hat{R}_0  \\
\pm\hat{R}_0 & \mp\hat{R}_0 & 0
\end{bmatrix}
\end{equation*}
$\Rightarrow \textbf{Extension-Twisting and Shearing-Bending (ETSB) Matrix}$
</P>

<p style="text-align: justify;">
ou alors si on veut annuler entièrement les effets de traction-torsion, il suffit de tourner la séquence des plis dans le stratifié de telle manière à obtenir $cos4\hat{\Phi}_0 = \pm 1$ ce qui implique que $sin4\hat{\Phi}_0 = 0$, en prenant par exemple ${\Phi}_0 = {\Phi}_1 = 0$ (orthotropie ordinaire)
</P>

<p style="text-align: justify;">
\begin{equation*}
\textbf{B}=
\begin{bmatrix}
\pm\hat{R}_0  & \mp\hat{R}_0  & 0  \\
\mp\hat{R}_0  & \pm\hat{R}_0  & 0  \\
0 & 0 & \mp\hat{R}_0 
\end{bmatrix}
\end{equation*}
$\Rightarrow \textbf{Extension-Bending and Shearing-Twisting (EBST) Matrix}$
</P>

<p style="text-align: justify;">
on alors ici un couplage en flexion-traction et cisaillement-torsion.
</P>

<p style="text-align: justify;">
2. Un autre aspect important de la stabilité thermo-élastique est lié à la forme du tenseur $\textbf{U}$ : celui-ci mesure les efforts thermiques dans le stratifié induits par une élévation uniforme de température. Est-ce que le tenseur $\textbf{U}$ peut être nul ?
</P>

<p style="text-align: justify;">
On écrit pour cela le tenseur \textbf{U} sous sa forme polaire la plus générale avec $T_u$ et $R_u e^{2i\Phi_u}$ :
</P>

<p style="text-align: justify;">
\begin{equation*}
    T_u=\big(\sum_{k=1}^{N}(z_k-z_{k-1})T_{\gamma}\big) R_u e^{2i\Phi_u} = \big(\sum_{k=1}^{N}(z_k-z_{k-1})R_{\gamma}e^{2i(\Phi_{\gamma}+\delta_k)}\big)
\end{equation*}
</P>

<p style="text-align: justify;">
on particularise pour des couches identiques :
</P>

$$
T_u= \frac{h}{N}\big(\sum_{k=1}^{N}T_{\gamma}\big) = hT_{\gamma}
$$

$$
R_u e^{2i\Phi_u} = \frac{h}{N}\big(\sum_{k=1}^{N}R_{\gamma}e^{2i(\Phi_{\gamma}+\delta_k}\big)=\frac{h}{N}R_{\gamma}e^{2i\Phi_{\gamma}}\big(\sum_{k=1}^{N}e^{2i\delta_k}\big)
$$

<p style="text-align: justify;">
Ainsi, un stratifié à couche identique donne $\textbf{U = 0}$, si et seulement si,
</P>

<p style="text-align: justify;">
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
</P>


<p style="text-align: justify;">
Comme déjà dit, $R_{\gamma} = 0$ signifie que $\boldsymbol{\gamma}$ est purement sphérique, c'est à dire que les propriétés du matériau de la couche de base du stratifié est à symétrie du carré ($R_1^{CB} = 0$) ou cas trivial i.e isotrope.
</P>


<p style="text-align: justify;">
Pour voir comment annuler $T_{\gamma}$ et $R_{\gamma}$, donnons avant tout leurs expressions à l'aide de la représentation polaire sur les cercles de Mohr.
</P>

<p style="text-align: justify;">
En effet $T_{\gamma}$ est l'invariant sphérique qui représente l'abscisse du centre du cercle de Mohr et $R_{\gamma}$ l'invariant déviatorique qui représente le rayon du cercle de Mohr, soit,
</P>

<p style="text-align: justify;">
\begin{equation*}
T_{\gamma} = \frac{\gamma_1 + \gamma_2}{2}
\end{equation*}
\begin{equation*}
R_{\gamma} = \sqrt{\Big(\frac{\gamma_1 - \gamma_2}{2}\Big)^2}
\end{equation*}
</P>

<p style="text-align: justify;">
Calculons ensuite les composante du tenseur $\boldsymbol{\gamma}$ tel que $\boldsymbol{\gamma}=\textbf{Q}\boldsymbol{\alpha}$, en prenant un comportement au moins de type orthotrope pour le tenseur de dilatation thermique $\boldsymbol{\alpha}$ et le tenseur de rigidité $\textbf{Q}$ exprimé dans la couche de base (i.e dans ces axes d'orthotropie $\rightarrow$ $Q_{16}$ = Q_{26} = 0$):
</P>

<p style="text-align: justify;">
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
</P>

ainsi en remplaçant, 

<p style="text-align: justify;">
\begin{equation*}
T_{\gamma} = \frac{1}{2}\big[\alpha_1(Q_{11} + Q_{12}) + \alpha_2(Q_{12} + Q_{22})\big]
\end{equation*}
</P>

<p style="text-align: justify;">
\begin{equation*}
R_{\gamma} = \frac{1}{2}[\alpha_1(Q_{11} - Q_{12}) + \alpha_2(Q_{12} - Q_{22})] 
\end{equation*}
</P>

<p style="text-align: justify;">
Dans le cas particulier où la couche de base est à symétrie du carré les différentes composantes s'identifient comme suit,
</P>

$$ \alpha_1 = \alpha_2 $$
$$ Q_{11} = Q_{22} $$

et on obtient,

\begin{equation*}
\fcolorbox{red}{ }{$R_{\gamma} = 0$} 
\end{equation*}

<p style="text-align: justify;">
On alors $\textbf{V}=0$ (condition matérielle à satisfaire cf. [\ref{condition}]) et aussi $\textbf{U}$ $\text{sphérique}$ 
</P>

<p style="text-align: justify;">
Et, pour obtenir $\textbf{U}$ nul pour les deux solutions possible il faut que $T_\gamma$ soit nul.
</P>

<p style="text-align: justify;">
Est-il possible d'obtenir  $T_{\gamma} = 0$ ?
En utilisant les expressions de la page 4 donnant les composantes polaire pour le tenseur $\textbf{B}$, il vient,
</P>

<p style="text-align: justify;">
\begin{align*}
Q_{11} + Q_{12} &= 4T_1 + 4R_1cos2\Phi_1\\ 
Q_{12} + Q_{22} &= 4T_1 - 4R_1cos2\Phi_1 
\end{align*}
Or $\Phi_1$ = 0 pour la couche de base, 
\begin{align*}
Q_{11} + Q_{12} &= 4T_1 + 4R_1\\
Q_{12} + Q_{22} &= 4T_1 - 4R_1
\end{align*}
</P>

<p style="text-align: justify;">
ce qui conduit à,

\begin{equation*}
T_{\gamma} = \frac{1}{2}\big[4\alpha_1(T_1 + R_1) + 4\alpha_2(T_1 - R_1)\big] = 0
\end{equation*}
</P>

<p style="text-align: justify;">
Étant donné que $T_1$, $\alpha_1$, $\alpha_2$ sont strictement positif dans le cas de matériau d'usage commun. De plus la dilatation thermique $\alpha_2$ défini dans le sens de la matrice est supérieur à $\alpha_1$ dans la direction des fibres. Donc la condition à satisfaire est telle que :
</P>

\begin{equation*}
\fcolorbox{red}{ }{$R_1 = T_1\frac{(\alpha_1+\alpha_2)}{(\alpha_2R_1 - \alpha_1)} > 0$}
\end{equation*}

<p style="text-align: justify;">
Les matériaux d'usage commun ne satisfont pas cette condition. Ainsi $T_{\gamma}$ n'est pas nul en générale, même si théoriquement cela est possible, et donc faisable seulement pour des matériaux $\textbf{ad-hoc}$. Donc $\textbf{U}$ est en générale non nul.
</P>

<p style="text-align: justify;">
Si on imagine une plaque stratifiée de forme rectangulaire dans un référentiel d’axes $x-y$, quelle forme du tenseur $\textbf{U}$ garantit que la plaque reste de forme rectangulaire ?
</P>

<p style="text-align: justify;">
Le tenseur des efforts linéique thermique par unité de température uniforme $\textbf{U}$ mesure l'élévation ou la diminution uniforme de température. Si la composante $\textbf{U}_6$ = $\textbf{U}_{xy}$ est nul alors les directions de $\textbf{U}$ sont alignés dans les directions des axes principales d'orthotropie de la couche.
</P>

<p style="text-align: justify;">
\begin{equation*}
 \textbf{U}=
   \begin{bmatrix}
U_{xx} \\
U_{yy}\\
0 
\end{bmatrix} 
\Rightarrow \textbf{Déformation dans le plan uniquement}
\end{equation*}
</P>

Ainsi, la plaque rectangulaire reste rectangulaire.

<p style="text-align: justify;">
Quel est le lien avec le tenseur A de rigidité en membrane ?
</P>

<p style="text-align: justify;">
Le tenseur $\textbf{A}$ d'ordre 4 possède 4 composantes polaires. L'expression de $\textbf{A}$ sous forme polaire est donc donné par $\hat{T}_0$, $\hat{T}_1$, $\hat{R}_0e^{4i\hat{\Phi}_0}$ et $\hat{R}_1e^{2i\hat{\Phi}_1}$ que l'ont écrit sous forme particularisé pour des couches identiques, on a pour $\hat{R}_1e^{2i\hat{\Phi}_1}$:
</P>

<p style="text-align: justify;">
\begin{equation*}
\hat{R}_1e^{2i\hat{\Phi}_1} =
  \frac{1}{2}\frac{h}{N}{R}_1^{CB}e^{2i\Phi_1^{CB}}\big(\sum_{k=1}^{N}e^{2i\delta_k}\big) 
\end{equation*}
</P>

<p style="text-align: justify;">
Le lien entre $\textbf{A}$ et $\textbf{U}$ est qu'ils dépendent tout les deux du terme $\fcolorbox{red}{ }{$\sum_{k=1}^{N}e^{2i\delta_k}$}.$
</P>

<p style="text-align: justify;">
Si on veut que la stabilité thermo-élastique soit étendue au comportement de dilatation thermique dans le plan, quelle condition faut-il imposer ?
</P>

<p style="text-align: justify;">
Si maintenant on impose de plus que $U_{xx}$=$U_{yy}$ c'est à dire $\textbf{U}$ sphérique,
</P>

<p style="text-align: justify;">
\begin{equation*}
 \textbf{U}=
   \begin{bmatrix}
U_{xx} \\
U_{xx}\\
0 
\end{bmatrix} 
\end{equation*}
</P>

<p style="text-align: justify;">
alors la déformation de la plaque rectangulaire est une plaque rectangulaire avec les mêmes proportion de dimension, i.e $\textbf{transformation homothétique}.$
</P>

<p style="text-align: justify;">
Quelle conséquence sur le comportement de rigidité en membrane (tenseur $\textbf{A}$) ?
</P>

<p style="text-align: justify;">
Puisque que cette condition de sphéricité pour $\textbf{U}$ est obtenue dans le cas le plus courant comme vue précédemment, lorsque
</P>

<p style="text-align: justify;">
\begin{equation*}
\sum_{k=1}^{N}e^{2i\delta_k} = 0
\end{equation*}

et d'après l'analyse précèdent le tenseur de rigidité en membrane $\textbf{A}$ $\textbf{sera à symétrie du carré}$, car $\bar{R}_1 = 0$ dans ce cas.
</P>

<p style="text-align: justify;">
A partir des expressions de la $CLPT$ écrites en polaire, écrire l’ensemble des conditions qui doivent être satisfaites afin d’obtenir la stabilité thermo-élastique (à la fois $\textbf{V} = O$ et $\textbf{U}$ sphérique) en fonction des angles de la stratification.
</P>

<li> 
<p style="text-align: justify;"> Pour avoir un stratifié a comportement thermo-élastiquement stable il suffit de prendre une couche de base à $R_1^{CB} = 0$, i.e à symétrie du carré, au quel cas $\textbf{V}=0$ et $\textbf{U}$ est sphérique quelque soit la séquence, mais $\textbf{A}$, $\textbf{B}$ et $\textbf{D}$ sont nécessairement aussi à symétrie du carré (limitant) et les tissus sont équilibrés (rare et coûteux par rapport aux UD).
</P>

<li> 
<p style="text-align: justify;"> L'autre solution est de prendre des couches UD, dans ce cas le stratifié a un comportement thermo-élastiquement stable s'il vérifie la séquence [$\delta_k$] tel que, $\textbf{V}=0$ $\Leftrightarrow$ $\sum_{k=1}^{N}b_ke^{2i\delta_k} = 0$ et vérifie U sphérique si $\sum_{k=1}^{N}e^{2i\delta_k}=0$.
</P><br>

$\textbf{Résumé}$ :

<p style="text-align: justify;">

<div style="overflow-x: auto; border-radius: 10px; background-color: #080a10; padding: 10px;">
  <table style="border-collapse: collapse; width: 100%; overflow: hidden; border-radius: 10px;">
    <thead style="background-color: #333; color: white; border-radius: 10px;">
      <tr>
        <th style="border: 1px solid white; padding: 10px ; border-top-left-radius: 10px;">&nbsp;</th>
        <th style="border: 1px solid white; padding: 10px;">$1^{ere}$ cas : propriété matérielle</th>
        <th style="border: 1px solid white; padding: 10px;">$2^{eme}$ cas : propriété séquentielle</th>
      </tr>
    </thead>
    <tbody style="color: white;">
      <tr>
        <td style="border: 1px solid white; padding: 10px; border-bottom-left-radius: 0px;"><strong>V</strong> = 0</td>
        <td style="border: 1px solid white; padding: 10px;">$R_1^{CB} = 0$</td>
        <td style="border: 1px solid white; padding: 10px;">$\sum_{k=1}^{N}b_ke^{2i\delta_k} = 0$</td>
      </tr>
      <tr>
        <td style="border: 1px solid white; padding: 10px;"><strong>U</strong> sphérique</td>
        <td style="border: 1px solid white; padding: 10px;">$R_1^{CB} = 0$</td>
        <td style="border: 1px solid white; padding: 10px;">$\sum_{k=1}^{N}e^{2i\delta_k} = 0$</td>
      </tr>
      <tr>
        <td style="border: 1px solid white; padding: 10px;"><strong>U</strong> nul</td>
        <td style="border: 1px solid white; padding: 10px;">ad-hoc</td>
        <td style="border: 1px solid white; padding: 10px;">$\varnothing$</td>
      </tr>
      <tr>
        <td style="border: 1px solid white; padding: 10px; border-top-right-radius: 20px;">conséquence</td>
        <td style="border: 1px solid white; padding: 10px;">$A, B, D$ à symétrie du carré</td>
        <td style="border: 1px solid white; padding: 10px; border-bottom-right-radius: 20px;">$B, A$ à symétrie du carré</td>
      </tr>
    </tbody>
  </table>
</div>


</P>

<p style="text-align: justify;">
Quel est le nombre minimal de couches nécessaire pour satisfaire l’ensemble de ces conditions ?
</P>

<p style="text-align: justify;">
Dans le cadre du $2^{eme}$ cas du tableau nous avons en réalité 4 conditions à satisfaire:
</P>

<p style="text-align: justify;">
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
</P>

<p style="text-align: justify;">
Soit 4 équations à résoudre donc il faut au minimum 4 ou 5 couches pour satisfaire l'ensemble de ces conditions.
</P><br>

<p style="text-align: justify;">
Comment résoudre le problème pour un nombre de couches plus élevé ?
</P><br>

<p style="text-align: justify;">
Pour résoudre un problème par exemple à 6 couches, on peut fixer une orientation de l'une des couches puis résoudre le problème pour 5 couches en prenant une couche parmi les cinq définissant le repère. Mais pour beaucoup plus de couches il faut résoudre un problème d'optimisation.
</P><br>

<p style="text-align: justify;">
Écrire aussi les conditions à satisfaire pour obtenir le découplage élastique et comparer avec les relations établies pour la stabilité thermo-élastique. Commenter
</P><br>

<p style="text-align: justify;">
Pour obtenir le découplage élastique il faut que $\textbf{B} = \textbf{0}$, cette condition est satisfaite si et seulement si :
</P><br>

<p style="text-align: justify;">
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
        T_0^{CB}\big(\sum_{k=1}^{N}b_k\big) = 0 \quad (\top) \\
        T_1^{CB}\big(\sum_{k=1}^{N}b_k\big) = 0 \quad  (\top)\\
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
</P>

( $\top$ : toujours vrai )  
<br>


<p style="text-align: justify;">
On s'aperçoit alors que deux des conditions sur quatre à satisfaire pour obtenir $\textbf{B} = \textbf{0}$ correspondent aux deux conditions à satisfaire pour obtenir $\textbf{V} = \textbf{0}$
</P>
<br>

<p style="text-align: justify;">
Une solution bien connue pour la stabilité thermo-élastique a été obtenue de manière empirique par $\textit{Winckler}$ (1985). Il s’agit de stratifiés à $n = 8$ couches caractérisés par la séquence :
</P>

$$
[\theta/(\theta - 90)_2/\theta/-\theta/(90 - \theta)_2/-\theta]
$$

### III. Stratifiés de $\textit{Winckler}$

<p style="text-align: justify;">
Vérifier que les stratifiés de $\textit{Winckler}$ sont thermo-élastiquement stables quelle que soit la valeur de l’angle $\theta$.
</P>

<p style="text-align: justify;">
Calcul de $b_k$ pour $n = 8$ couches, avec $b_k = 2k-n-1$ :
</P>

<p style="text-align: justify;">
\begin{align*}
    &b_1 = -7 ; \\
    &b_2 = -5 ; \\
    &b_3 = -3 ; \\
    &b_4 = -1 ; \\
    &b_5 = 1 ; \\
    &b_6 = 3 ; \\
    &b_7 = 5 ; \\
    &b_8 = 7
\end{align*}
</P>

<p style="text-align: justify;">
Calcul de $\sum_{k=1}^{N}b_kcos{2\delta_k}$ pour la séquence de $\textit{Wincker}$,
</P>

<p style="text-align: justify;">
\begin{align*}
\sum_{k=1}^{N}b_kcos{2\delta_k} &= -7cos2\theta - 5cos2(\theta - 90) - 3cos2(\theta - 90) - cos2\theta\\
&+ cos2(-\theta) + 3cos2(90 - \theta) + 5cos2(90 - \theta) + 7cos2(-\theta) 
\end{align*}
or $cos(x) = cos(-x)$
\begin{align*}
\sum_{k=1}^{N}b_kcos{2\delta_k} &= -7cos2\theta - 5cos2(\theta - 90) - 3cos2(\theta - 90) - cos2\theta\\
&+ cos2(\theta) + 3cos2(\theta - 90) + 5cos2(\theta - 90 ) + 7cos2(\theta) = 0
\end{align*}
</P>

donc,

<p style="text-align: justify;">
\begin{equation*}
\fcolorbox{red}{ }{$\sum_{k=1}^{N}b_kcos{2\delta_k} = 0 \quad \forall \theta$}
\end{equation*}
</P>

<p style="text-align: justify;">
Calcul de $\sum_{k=1}^{N}b_ksin{2\delta_k}$ pour la séquence de \textit{Winckler},
</P>

<p style="text-align: justify;">
\begin{align*}
\sum_{k=1}^{N}b_ksin{2\delta_k} &= -7sin2\theta - 5sin2(\theta - 90) - 3sin2(\theta - 90) - sin2\theta\\
&+ sin2(-\theta) + 3sin2(90 - \theta) + 5sin2(90 - \theta) + 7sin2(-\theta)
\end{align*}
or $sin(-x) = -sin(x)$, $sin(x-90) = -cos(x)$ et $sin(2x) = 2sin(x)cos(x)$, donc
\begin{align*}
\sum_{k=1}^{N}b_ksin{2\delta_k} &= -16sin2\theta - 16sin2(\theta - 90)\\
&= -32sin\theta cos\theta + 32sin\theta cos\theta = 0 
\end{align*}
</P>

donc, 

<p style="text-align: justify;">
\begin{equation*}
 \fcolorbox{red}{ }{$\sum_{k=1}^{N}b_ksin{2\delta_k} = 0 \quad \forall \theta$}
\end{equation*}
</P>

<p style="text-align: justify;">
Ainsi les deux conditions à satisfaire pour avoir $\textbf{V} = \textbf{0}$ sont vérifiées.
</P>

<p style="text-align: justify;">
De la même façon en calculant  $\sum_{k=1}^{N}cos{2\delta_k}$ et  $\sum_{k=1}^{N}sin{2\delta_k}$ pour la séquence de $\textit{Winckler}$ ont trouve que ces sommes s'annule 
</P>

<p style="text-align: justify;">
\begin{equation*}
  \fcolorbox{red}{ }{$\sum_{k=1}^{N}cos{2\delta_k} = 0$} \quad ; \quad  \fcolorbox{red}{ }{$\sum_{k=1}^{N}sin{2\delta_k} = 0$}
\end{equation*}
</P>

<p style="text-align: justify;">
et satisfont donc les conditions pour que $\textbf{U}$ soit sphérique.
</P><br>

<p style="text-align: justify;">
$\textbf{Conclusion}$ : La séquence de $\textit{Winckler}$ vérifie les conditions pour que le stratifié soit thermo-elastiquement stable ($\textbf{V} = \textbf{0}$) et pour que les dilatations thermique s'effectue sphériquement dans le plan ($\textbf{U}$ sphérique).
</P><br>

<p style="text-align: justify;">
5. Calculer les propriétés de rigidité en membrane et flexion, $\textbf{A}$ et $\textbf{D}$, ainsi que le couplage élastique $\textbf{B}$ pour ces stratifiés. Commenter sur la forme et les propriétés de ces comportements.
</P><br>

On a vu que la séquence de $\textit{Winckler}$ permet d'avoir,

<p style="text-align: justify;">
\begin{equation*}
  \sum_{k=1}^{N}e^{2i\delta_k} = 0 \quad \text{et} \quad \sum_{k=1}^{N}b_ke^{2i\delta_k} = 0
\end{equation*}
</P>

<p style="text-align: justify;">
or le terme en $R_1$ de la forme polaire de $\textbf{A}$ et $\textbf{B}$ s'écrivent respectivement comme suit,
</P>

<p style="text-align: justify;">
\begin{equation*}
   \bar{R}_1e^{2i\bar{\Phi}_1} = \frac{h}{N}R_1^{CB}e^{2i\Phi_1^{CB}}\sum_{k=1}^{N}e^{2i\delta_k} \quad ;\quad \hat{R}_1e^{2i\hat{\Phi}_1} = \frac{1}{2}\frac{h^2}{N^2}R_1^{CB}e^{2i\Phi_1^{CB}}\sum_{k=1}^{N}b_ke^{2i\delta_k} 
\end{equation*}
</P>

donc,

<p style="text-align: justify;">
\begin{equation*}
    \fcolorbox{red}{ }{$\bar{R}_1=0$} \quad ;\quad  \fcolorbox{red}{ }{$\hat{R}_1 = 0$}
\end{equation*}
</P>

<p style="text-align: justify;">
On en conclut que $\textbf{A}$ et $\textbf{B}$ ont au moins comme propriété la symétrie du carré.
</P><br>

<p style="text-align: justify;">
Pour s'en assurer on calcul $\sum_{k=1}^{N}cos{4\delta_k}$ pour la séquence de $\textit{Wincker}$ qui nous donne,
</P>

<p style="text-align: justify;">
\begin{equation*}
\sum_{k=1}^{N}cos{4\delta_k} = 8cos4\theta \, ,
\end{equation*}
</P>

<p style="text-align: justify;">
le calcul de $\sum_{k=1}^{N}sin{4\delta_k}$ pour la séquence de $\textit{Winckler}$, nous donne
</P>

<p style="text-align: justify;">
\begin{equation*}
\sum_{k=1}^{N}sin{4\delta_k} = 0 \quad \forall \theta \, ,
\end{equation*}
</P>

calcul de $\sum_{k=1}^{N}b_kcos{4\delta_k}$ pour la séquence de $\textit{Winckler}$, nous donne

<p style="text-align: justify;">
\begin{equation*}
\sum_{k=1}^{N}b_kcos{4\delta_k} = 0 \quad \forall \theta 
\end{equation*}
</P>

<p style="text-align: justify;">
et calcul de $\sum_{k=1}^{N}b_ksin{4\delta_k}$ pour la séquence de \textit{Winckler}, nous donne
</P>

<p style="text-align: justify;">
\begin{equation*}
\sum_{k=1}^{N}b_ksin{4\delta_k} = -32sin4\theta
\end{equation*}
</P>

<p style="text-align: justify;">
or le terme en $R_0$ de la forme polaire de $\textbf{A}$ et 
$\textbf{B}$ s'écrivent respectivement comme suit,
</P>

<p style="text-align: justify;">
\begin{equation*}
   \bar{R}_0e^{4i\bar{\Phi}_0} = \frac{h}{N}R_0^{CB}e^{4i\Phi_0^{CB}}\sum_{k=1}^{N}e^{4i\delta_k} \quad ;\quad \hat{R}_0e^{4i\hat{\Phi}_0} = \frac{1}{2}\frac{h^2}{N^2}R_0^{CB}e^{4i\Phi_0^{CB}}\sum_{k=1}^{N}b_ke^{4i\delta_k} 
\end{equation*}
</P>
donc,
<p style="text-align: justify;">
\begin{equation*}
   \fcolorbox{red}{ }{$\bar{R}_0e^{4i\bar{\Phi}_0} = \frac{h}{N}8R_0^{CB}e^{4i\Phi_0^{CB}}cos4\theta$}  \quad ;\quad \fcolorbox{red}{ }{$\hat{R}_0e^{4i\hat{\Phi}_0}=-\frac{1}{2}\frac{h^2}{N^2}32R_0^{CB}e^{4i\Phi_0^{CB}}sin4\theta$}
\end{equation*}
</P>

<p style="text-align: justify;">
On en conclut que $\textbf{A}$ et $\textbf{B}$ sont donc au moins à symétrie du carré.
</P>

Qu'en est - il de $\textbf{D}$ ?

<p style="text-align: justify;">
Calcul de $d_k$ pour n = 8 couches, avec $d_k = 12k(k-n-1)+4+3n(n+2)$ :
</P>


<p style="text-align: justify;">
\begin{align*}
    &d_1 = 148 \\
    &d_2 = 76  \\
    &d_3 = 28  \\
    &d_4 = 4 \\
    &d_5 = 28  \\
    &d_7 = 76  \\
    &d_8 = 148
\end{align*}
</P>

<p style="text-align: justify;">
Le calcul de $\sum_{k=1}^{N}d_kcos{2\delta_k}$ pour la séquence de $\textit{Wincker}$, nous donne
</P>

<p style="text-align: justify;">
\begin{equation*}
\sum_{k=1}^{N}d_kcos{2\delta_k} = 96cos2\theta 
\end{equation*}
</P>

<p style="text-align: justify;">
Le calcul de $\sum_{k=1}^{N}d_ksin{2\delta_k}$ pour la séquence de $\textit{Winckler}$, nous donne
</P>

<p style="text-align: justify;">
\begin{equation*}
\sum_{k=1}^{N}d_ksin{2\delta_k} = 0  \quad \forall \theta
\end{equation*}
</P>

<p style="text-align: justify;">
Le calcul de $\sum_{k=1}^{N}d_kcos{2\delta_k}$ pour la séquence de $\textit{Wincker}$, nous donne
</P>

<p style="text-align: justify;">
\begin{equation*}
\sum_{k=1}^{N}d_kcos{4\delta_k} = 512cos4\theta 
\end{equation*}
</P>

<p style="text-align: justify;">
Le calcul de $\sum_{k=1}^{N}d_ksin{2\delta_k}$ pour la séquence de \textit{Winckler}, nous donne
\begin{equation*}
\sum_{k=1}^{N}d_ksin{4\delta_k} = 0  \quad \forall \theta
\end{equation*}
</P>

<p style="text-align: justify;">
or les termes en $R_0$ et $R_1$ de la forme polaire de $\textbf{D}$ s'écrivent respectivement comme suit,
</P>

<p style="text-align: justify;">
\begin{equation*}
  \tilde{R}_0e^{4i\tilde{\Phi}_0} = \frac{1}{12}\frac{h^3}{N^3}R_0^{CB}e^{4i\Phi_0^{CB}}\sum_{k=1}^{N}d_ke^{4i\delta_k} \quad ;\quad  
  \tilde{R}_1e^{2i\tilde{\Phi}_1} = \frac{1}{12}\frac{h^3}{N^3}R_0^{CB}e^{2i\Phi_1^{CB}}\sum_{k=1}^{N}d_ke^{2i\delta_k}
\end{equation*}
</P>

<p style="text-align: justify;">
$\Rightarrow$
\begin{equation*}
   \fcolorbox{red}{ }{$\tilde{R}_0e^{4i\tilde{\Phi}_0} = \frac{1}{12}\frac{h^3}{N^3}R_0^{CB}e^{4i\Phi_0^{CB}}512cos4\theta$} \quad ;\quad \fcolorbox{red}{ }{$\tilde{R}_1e^{2i\tilde{\Phi}_1} = \frac{1}{12}\frac{h^3}{N^3}R_1^{CB}e^{2i\Phi_1^{CB}}96cos2\theta$}
\end{equation*}
</P>

<p style="text-align: justify;">
On en conclut que $\textbf{D}$ est de propriété orthotrope ordinaire.
</P><br>

<p style="text-align: justify;">
6. Ecrire la forme cartésienne du tenseur  $\textbf{B}$. Commenter sur la forme de ce tenseur et le type de couplages qu’il représente.
</P><br>

<p style="text-align: justify;">
Pour cela on utilise les expressions de passage polaire-cartésien de la page 4 de ce projet, que l'on particularise au cas où $\hat{T}_0=\hat{T}_1=\hat{R}_1=0$ et $\hat{R}_0 \neq 0$, le tenseur $\textbf{B}$ s'écrit alors:
</P><br>

<p style="text-align: justify;">
\begin{equation*}
\textbf{B} =
\begin{pmatrix}
\hat{R}_0cos4\hat{\Phi}_0 & -\hat{R}_0cos4\hat{\Phi}_0 & \hat{R}_0sin4\hat{\Phi}_0  \\
-\hat{R}_0cos4\hat{\Phi}_0  & \hat{R}_0cos4\hat{\Phi}_0 & -\hat{R}_0sin4\hat{\Phi}_0 \\
\hat{R}_0sin4\hat{\Phi}_0  & -\hat{R}_0sin4\hat{\Phi}_0 & -\hat{R}_0cos4\hat{\Phi}_0 
\end{pmatrix}
\end{equation*}
</P>

<p style="text-align: justify;">
Or nous avons vu dans la première partie que pour obtenir un couplage traction-torsion on peut imposer $sin4\hat{\Phi}_0 = \pm 1$ avec ${\Phi}^{CB}_0 = \frac{\pi}{8}$ ce qui implique automatiquement la correspondance $cos4\hat{\Phi}_0 = 0$.
% Dans le cas général de matériau orthotrope $\Phi_1 = 0^o$ dans la couche de base et que $\Phi_1 = 0^o$-$\Phi_1 = k\frac{\pi}{4}$, on prend ici $k=0$. Soit, $sin4\hat{\Phi}_0 = \pm 1$ avec ce qui implique automatiquement la correspondance $cos4\hat{\Phi}_0 = 0$ 
</P>

<p style="text-align: justify;">
\begin{equation*}
\textbf{B}=
\begin{bmatrix}
0 & 0 &  \hat{R}_0sin4\hat{\Phi}_0 \\
0 & 0 & - \hat{R}_0sin4\hat{\Phi}_0  \\
\hat{R}_0sin4\hat{\Phi}_0 & - \hat{R}_0sin4\hat{\Phi}_0 & 0
\end{bmatrix}
\end{equation*}
</P>

<p style="text-align: justify;">
ainsi cette composante du tenseur vaut $\hat{R}_0sin4\hat{\Phi}_0  = - \frac{1}{2}\frac{h^2}{N^2}32R_0^{CB}sin4\theta$, fixons maitenant le paramètre $\hat{p}_0$ tel que :
</P>

<p style="text-align: justify;">
\begin{equation*}
\hat{p}_0=-\frac{1}{2}\frac{h^2}{N^2}32R_0^{CB}
\end{equation*}
</P>

il vient,

<p style="text-align: justify;">
\begin{equation*}
\textbf{B}=
\begin{bmatrix}
0 & 0 &  \hat{p}_0sin4\theta  \\
0 & 0 & - \hat{p}_0sin4\theta \\
\hat{p}_0sin4\theta& - \hat{p}_0sin4\theta & 0
\end{bmatrix}
\end{equation*}
</P>

<p style="text-align: justify;">
On obtient alors un tenseur de couplage de type $\textbf{Extension-Twisting and Shearing-Bending}$ (ETSB), il permet alors d'avoir un couplage élastique du stratifié en "traction-torsion" et/ou cisaillement-flexion, comme vu en première partie, mais pas seulement...
% ou alors si on veut annuler les entièrement les effets de traction-torsion, on tourne la séquence des plis dans le stratifié de telle manière à obtenir $cos4\hat{\Phi}_0 = \pm 1$ ce qui implique que $sin4\hat{\Phi}_0 = 0$.  
</P><br>  

<p style="text-align: justify;">
En effet, par rapport au cas traité en première partie, en plus d'avoir un stratifié thermo-elastiquement stable ($\textbf{V} = 0$ et $\textbf{U}$ sphérique) la séquences de $\textit{Winckler}$ permet de fournir en plus $\textbf{la possibilité de moduler l'intensité du couplage}$ élastique $\textbf{B}$ à l'aide du paramètre géométrique $\theta$ désiré.
</P><br>
  
<p style="text-align: justify;">
7. Comment varient les propriétés élastiques des stratifiés de \textit{Winckler} en fonction de l’angle $\theta$? Dire pour quelle valeur de l’angle $\theta$ on obtient le couplage $\textbf{B}$ maximum. Pour cette valeur  $\theta_{max}$ écrire aussi les expressions des tenseurs $\textbf{A}$ et $\textbf{D}$. Tracer les courbes représentant les composantes $A_11$ et $D_11$ en fonction de l’angle du référentiel.
</P><br>  

<p style="text-align: justify;">
Les propriétés élastiques des stratifiés de Winckler en fonction de l’angle $\theta$ varie comme $4$ fois l'angle, on parle alors d'harmonique ou de $\textbf{periodicité d'ordre 4}$. En effet, $sin\theta$ va de $0$ à $1$ entre $[0;\frac{\pi}{2}]$ alors qu'ici pour cette intervalle on a,
</P><br>  

<p style="text-align: justify;">
\begin{align*}
&\theta = 0 \Rightarrow  sin4\theta = 0\\
&\theta = \frac{\pi}{8} \Rightarrow  sin4\theta = 1\\
&\theta = \frac{\pi}{4} \Rightarrow  sin4\theta = 0\\
&\theta = \frac{3\pi}{8} \Rightarrow  sin4\theta = -1\\
&\theta = \frac{\pi}{2} \Rightarrow  sin4\theta = 0
\end{align*}
</P>

On en conclut que la valeur de $\theta$ qui donne un couplage $\textbf{B}$ maximum est,

<p style="text-align: justify;">
\begin{equation*}
\theta_{max} = \frac{\pi}{8} \qquad [\frac{2\pi}{8}]
\end{equation*}
</P>

<p style="text-align: justify;">
le signe $\pm$ ne change que le sens de la déformation géométriquement.
Écrivons maintenant les expressions des tenseurs $\textbf{A}$ et $\textbf{D}$ pour la valeur  $\theta$ = $\theta_{max}$ : 
</P>

<p style="text-align: justify;">
Pour $\textbf{A}$, nous avons vu que $\bar{R}_1 = 0$, et nous avons aussi $\bar{T}_0 = hT_0^{CB}$ et $\bar{T}_1 = hT_1^{CB}$.
</P>

<p style="text-align: justify;">
\begin{equation*}
\bar{R}_0e^{4i\bar{\Phi}_0} = \underbrace{\frac{h}{N}8R_0^{CB}}_{\bar{p}_0}cos4\theta(cos4\Phi_0^{CB} + isin4\Phi_0^{CB})    
\end{equation*}
</P>

le tenseur $\textbf{A}$ s'écrit alors:

<p style="text-align: justify;">
\begin{equation*}
\textbf{A} =
\begin{pmatrix}
\bar{T}_0 + 2\bar{T}_1 + \bar{R}_0cos4\bar{\Phi}_0 & -\bar{T}_0 + 2\bar{T}_1 -\bar{R}_0cos4\bar{\Phi}_0 & \bar{R}_0sin4\bar{\Phi}_0    \\
 -\bar{T}_0 + 2\bar{T}_1 - \bar{R}_0cos4\bar{\Phi}_0 & \bar{T}_0 + 2\bar{T}_1 + \bar{R}_0cos4\bar{\Phi}_0  & -\bar{R}_0sin4\bar{\Phi}_0  \\
 \bar{R}_0sin4\bar{\Phi}_0 & -\bar{R}_0sin4\bar{\Phi}_0  & \bar{T}_0 -\bar{R}_0sin4\bar{\Phi}_0 
\end{pmatrix}
\end{equation*}
</P>

<p style="text-align: justify;">
Traçons alors la composante $A_{11}$ pour $\Phi_0^{CB}=0$ : $A_{11} = h(T_0^{CB}+ 2T_1^{CB}) +  \bar{p}_0cos4\theta$ 
</P><br>  

<div style="text-align: center;">
  <figure style="display: inline-block;">
    <img src="/assets/projet3/p3f2.png" alt="Représentation matérielle graphique" width="600"/>
    <figcaption>Figure 2 : Représentation matérielle graphique</figcaption>
  </figure>
</div><br>

<p style="text-align: justify;">
On en conclut que $A_{11}$ est dans ce cas à $\textbf{symétrie du carré}$, étant donné que la rigidité menbranaire est identique dans les deux directions principales. Or, dans le cadre d'un couplage traction-torsion nous avons imposé lors du calcul de $\textbf{B}$ que $\Phi_0^{CB}=\frac{\pi}{8}$.

Traçons alors la composante $A_{11}$ dans ce cas $\rightarrow$ $A_{11} = h(T_0^{CB}+ 2T_1^{CB})$
</P><br>  

<div style="text-align: center;">
  <figure style="display: inline-block;">
    <img src="/assets/projet3/p3f3.png" alt="Représentation matérielle graphique" width="600"/>
    <figcaption>Figure 3 : Représentation matérielle graphique</figcaption>
  </figure>
</div><br>

<p style="text-align: justify;">
On en conclut que $A_{11}$ est dans ce cas à $\textbf{isotrope}$, étant donné que la rigidité menbranaire est identique dans toutes les directions.<br>
Donc dans le cadre du couplage traction-torsion $\textbf{A}$ s'écrit:
</P><br>


<p style="text-align: justify;">
\begin{equation*}
\textbf{A} =
\begin{pmatrix}
\bar{T}_0 + 2\bar{T}_1 & -\bar{T}_0 + 2\bar{T}_1 & \bar{R}_0sin4\bar{\Phi}_0    \\
 -\bar{T}_0 + 2\bar{T}_1 & \bar{T}_0 + 2\bar{T}_1  & -\bar{R}_0sin4\bar{\Phi}_0  \\
 \bar{R}_0sin4\bar{\Phi}_0 & -\bar{R}_0sin4\bar{\Phi}_0  & \bar{T}_0
\end{pmatrix}
\end{equation*}
</div>

soit,

<p style="text-align: justify;">
\begin{equation*}
\textbf{A} =
\begin{pmatrix}
h(T_0^{CB}+ 2T_1^{CB}) & h(2T_1^{CB}-T_0^{CB})  &  \bar{p}_0cos4\theta \\
h(2T_1^{CB}-T_0^{CB})  & h(T_0^{CB}+ 2T_1^{CB}) & -\bar{p}_0cos4\theta  \\
\bar{p}_0cos4\theta & -\bar{p}_0cos4\theta  & hT_0^{CB} 
\end{pmatrix}
\end{equation*}
</div>

<p style="text-align: justify;">
Et pour $\theta_{max} = \frac{\pi}{8}$ on a $cos4\theta = 0$, donc $\bar{R}_0= 0$, le tenseur $\textbf{A}$ est alors totalement isotrope et s'écrit,
</P><br>

<p style="text-align: justify;">
\begin{equation*}
\textbf{A}=
\begin{pmatrix}
h(T_0^{CB}+ 2T_1^{CB}) &  h(2T_1^{CB}-T_0^{CB})  & 0  \\
 h(2T_1^{CB}-T_0^{CB})  & h(T_0^{CB} + 2T_1^{CB}) & 0 \\
0 & 0 & hT_0^{CB}
\end{pmatrix}
\end{equation*}
</div>

<p style="text-align: justify;">
Pour $\textbf{D}$, nous avons $\tilde{T}_0 = \frac{h^3}{3}\frac{T_0^{CB}}{4}$ et $\tilde{T}_1 = \frac{h^3}{3}\frac{T_1^{CB}}{4}$ et $\Phi_0=\frac{\pi}{8}$ et $\Phi_1 = 0$
</P>

<p style="text-align: justify;">
\begin{equation*}
   \tilde{R}_0e^{4i\tilde{\Phi}_0} = \underbrace{\frac{1}{12}\frac{h^3}{N^3}R_0^{CB}512}_{\tilde{p}_0}e^{4i\Phi_0^{CB}}cos4\theta \quad ;\quad \tilde{R}_1e^{2i\tilde{\Phi}_1} = \underbrace{\frac{1}{12}\frac{h^3}{N^3}R_1^{CB}96}_{\tilde{p}_1}e^{2i\Phi_1^{CB}}cos2\theta  
\end{equation*}
</P>

<p style="text-align: justify;">
\begin{align*}
&D_{11} = \tilde{T}_0 + 2\tilde{T}_1 + \tilde{R}_0cos4\tilde{\Phi}_0 + 4\tilde{R}_1cos2\tilde{\Phi}_1  \\
&D_{22} = \tilde{T}_0 + 2\tilde{T}_1 + \tilde{R}_0cos4\tilde{\Phi}_0 - 4\tilde{R}_1cos2\tilde{\Phi}_1\\
&D_{12} = -\tilde{T}_0 + 2\tilde{T}_1 - \tilde{R}_0cos4\tilde{\Phi}_0 \\
&D_{16} = \tilde{R}_0sin4\tilde{\Phi}_0 + 4\tilde{R}_1sin2\tilde{\Phi}_1 \\
&D_{26} = -\tilde{R}_0sin4\tilde{\Phi}_0 + 4\tilde{R}_1sin2\tilde{\Phi}_1  \\
&D_{66} =  \tilde{T}_0 - \tilde{R}_0cos4\tilde{\Phi}_0 
\end{align*}
</P>

<p style="text-align: justify;">
Traçons alors la composante $D_{11}$ pour $\Phi_0^{CB}=\Phi_1^{CB}=0$ : $D_{11} = \frac{h^3}{12}(T_0^{CB} + 2T_1^{CB}) +  \tilde{p}_0cos4\theta + 4\tilde{p}_1cos2\theta$ :
</P><br>  

<div style="text-align: center;">
  <figure style="display: inline-block;">
    <img src="/assets/projet3/p3f4.png" alt="Représentation matérielle graphique" width="600"/>
    <figcaption>Figure 4 : Représentation matérielle graphique</figcaption>
  </figure>
</div><br>

<p style="text-align: justify;">
On en conclut que $D_{11}$ est bien $\textbf{orthotrope ordinaire}$, étant donné que la rigidité en flexion est plus grande dans la direction 1.<br>
Or, dans le cadre d'un couplage traction-torsion nous avons imposé lors du calcul de $\textbf{B}$ que $\Phi_0^{CB}=\frac{\pi}{8}$.<br>  

Traçons alors la composante $D_{11}$ dans ce cas $\rightarrow$ $D_{11} = \frac{h^3}{12}(T_0^{CB} + 2T_1^{CB}) + 4\tilde{p}_1cos2\theta$
</P><br>  

<div style="text-align: center;">
  <figure style="display: inline-block;">
    <img src="/assets/projet3/p3f5.png" alt="Représentation matérielle graphique" width="600"/>
    <figcaption>Figure 5 : Représentation matérielle graphique</figcaption>
  </figure>
</div><br>  

<p style="text-align: justify;">
On en conclut que $D_{11}$ est bien $\textbf{orthotrope ordinaire}$.
</P>

<p style="text-align: justify;">
Donc, dans le cadre du couplage traction-torsion on a:
</P>

<p style="text-align: justify;">
\begin{equation*}
   \tilde{R}_0e^{4i\tilde{\Phi}_0} = \underbrace{\frac{1}{12}\frac{h^3}{N^3}R_0^{CB}512}_{\tilde{p}_0}cos4\theta \quad ;\quad \tilde{R}_1e^{2i\tilde{\Phi}_1} = \underbrace{\frac{1}{12}\frac{h^3}{N^3}R_1^{CB}196}_{\tilde{p}_1}cos2\theta  
\end{equation*}
</P>

et $\textbf{D}$ s'écrit:

<p style="text-align: justify;">
\begin{equation*}
\textbf{D} =
\begin{pmatrix}
\tilde{T}_0 + 2\tilde{T}_1 + 4\tilde{R}_1cos2\tilde{\Phi}_1 & -\tilde{T}_0 + 2\tilde{T}_1 & \tilde{R}_0sin4\tilde{\Phi}_0 \\
-\tilde{T}_0 + 2\tilde{T}_1 & \tilde{T}_0 + 2\tilde{T}_1 - 4\tilde{R}_1cos2\tilde{\Phi}_1  & -\tilde{R}_0sin4\tilde{\Phi}_0 \\
\tilde{R}_0sin4\tilde{\Phi}_0  & -\tilde{R}_0sin4\tilde{\Phi}_0 & \tilde{T}_0
\end{pmatrix}
\end{equation*}
</P>

soit, 

<p style="text-align: justify;">
\begin{equation*}
\textbf{D} =
\begin{pmatrix}
\frac{h^3}{12}(T_0^{CB} + 2T_1^{CB}) + 4\tilde{p}_1cos2\theta  & -\frac{h^3}{12}(T_0^{CB} - 2T_1^{CB}) & \tilde{p}_0cos4\theta \\
-\frac{h^3}{12}(T_0^{CB} - 2T_1^{CB})  & \frac{h^3}{12}(T_0^{CB} + 2T_1^{CB}) - 4\tilde{p}_1cos2\theta  & - \tilde{p}_0cos4\theta\\
\tilde{p}_0cos4\theta & - \tilde{p}_0cos4\theta & \frac{h^3}{12}T_0^{CB} 
\end{pmatrix}
\end{equation*}
</P>

<p style="text-align: justify;">
Pour $\theta_{max} = \frac{\pi}{8}$ on a $cos4\theta = 0$,
le tenseur $\textbf{D}$ est de propriété d'horthotropie $R_0$ (car $R_0 = 0$), et s'écrit :
</P>

<p style="text-align: justify;">
\begin{equation*}
\textbf{D} =
\begin{pmatrix}
\frac{h^3}{12}(T_0^{CB} + 2T_1^{CB}) + 4\tilde{p_1}\frac{\sqrt{2}}{2} & -\frac{h^3}{12}(T_0^{CB} - 2T_1^{CB}) & 0  \\
-\frac{h^3}{12}(T_0^{CB} - 2T_1^{CB}) & \frac{h^3}{12}(T_0^{CB} + 2T_1^{CB}) + 4\tilde{p_1}\frac{\sqrt{2}}{2} & 0 \\
0 & 0 & \frac{h^3}{12}T_0^{CB}  
\end{pmatrix}
\end{equation*}
</P>

<p style="text-align: justify;">
 Remarque : pour tracer les courbes à l'aide de $\textbf{Matlab}$ (cf. $\textit{Application.m}$) nous avons considéré un stratifié composé de huit couche unidirectionnelles en carbone/epoxyde d'épaisseur h, et modules polaires donnée dans le tableau ci-dessous:
</P>

<p style="text-align: justify;">
<div style="overflow-x: auto; border-radius: 10px; background-color: #080a10; padding: 10px;">
<table style="border-collapse: collapse; width: 100%; color: white;">
  <thead>
    <tr style="background-color: #333;">
      <th style="border: 1px solid white; padding: 10px; border-top-left-radius: 10px;">Modules polaires</th>
      <th style="border: 1px solid white; padding: 10px;">$T_0$ (GPa)</th>
      <th style="border: 1px solid white; padding: 10px;">$T_1$ (GPa)</th>
      <th style="border: 1px solid white; padding: 10px;">$R_0$ (GPa)</th>
      <th style="border: 1px solid white; padding: 10px;">$R_1$ (GPa)</th>
      <th style="border: 1px solid white; padding: 10px;">$\Phi_0,\Phi_1=0$</th>
      <th style="border: 1px solid white; padding: 10px; border-top-right-radius: 10px;">h (mm)</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td style="border: 1px solid white; padding: 10px; border-bottom-left-radius: 10px;"></td>
      <td style="border: 1px solid white; padding: 10px;">26,88</td>
      <td style="border: 1px solid white; padding: 10px;">24,74</td>
      <td style="border: 1px solid white; padding: 10px;">19,71</td>
      <td style="border: 1px solid white; padding: 10px;">21,43</td>
      <td style="border: 1px solid white; padding: 10px;">0</td>
      <td style="border: 1px solid white; padding: 10px; border-bottom-right-radius: 10px;">1</td>
    </tr>
  </tbody>
</table>
</div>
</p><br> 

<p style="text-align: justify;">
8. En réalité, pour comprendre la réponse d’un stratifié à un chargement mécanique, il faut étudier les tenseurs de souplesse : inverser la loi de comportement thermo-élastique, puis particulariser les équations dans le cas d’un chargement purement mécanique.<br> 
</P>

<p style="text-align: justify;">
Je ne veux pas mettre ici toutes les étapes du calcul. Mais la démarche est la suivante on part de, 
</P>

<p style="text-align: justify;">
\begin{equation*}\left\{
    \begin{array}{ll}
        \textbf{N} = \textbf{A}\boldsymbol{\mathcal{E}}^0 + \textbf{B}\boldsymbol{\mathcal{K}} - T_0\textbf{U} \\
        \textbf{M} = \textbf{B}\boldsymbol{\mathcal{E}}^0 + \textbf{D}\boldsymbol{\mathcal{K}}
    \end{array}
\right.
\end{equation*}
</P>

<p style="text-align: justify;">
en prenant la première équation, on isole $\mathcal{E}^0$ en pré-multipliant par $\textbf{A}^{-1}$. On obtient alors:
</P>

<p style="text-align: justify;">
\begin{equation*}
\boldsymbol{\mathcal{K}} = \textbf{d}\textbf{M} + \textbf{b}_2(\textbf{N}-T_0\textbf{U})
\end{equation*}
</P>

avec,

<p style="text-align: justify;">
\begin{align*}
&\textbf{d} = (\textbf{D}-\textbf{B}\textbf{A}^{-1}\textbf{B})^{-1}\\
&\textbf{b}_2 = -\textbf{d}\textbf{B}\textbf{A}^{-1}
\end{align*}
\noindent
de même pour la deuxième équation, on a
\begin{equation*}
\boldsymbol{\mathcal{E}}^0 = \textbf{a}(\textbf{N}-T_0\textbf{U}) + \textbf{b}_1\textbf{M} 
\end{equation*}
</P>

avec,

<p style="text-align: justify;">
\begin{align*}
&\textbf{a} = (\textbf{A}-\textbf{B}\textbf{D}^{-1}\textbf{B})^{-1}\\
&\textbf{b}_1 = -\textbf{a}\textbf{B}\textbf{D}^{-1}
\end{align*}
</P>

<p style="text-align: justify;">
On regroupe alors sous forme matricielle et on peut aussi montrer que $\textbf{b}_1 = \textbf{b}^T_2$ et posé $\textbf{b}_1 = \textbf{b}$
</P>

<p style="text-align: justify;">
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
</P>

<p style="text-align: justify;">
\begin{equation*}
\Leftrightarrow \; \; \left\{
    \begin{array}{ll}
        \boldsymbol{\mathcal{E}}^0 = \textbf{a}(\textbf{N}-T_0\textbf{U}) + \textbf{b}\textbf{M} \\
        \boldsymbol{\mathcal{K}} = \textbf{d}\textbf{M} + \textbf{b}^T(\textbf{N}-T_0\textbf{U})
    \end{array}
\right.
\end{equation*}
</P>

<p style="text-align: justify;">
Particularisons à un chargement purement mécanique,
</P>

<p style="text-align: justify;">
\begin{equation*}
\Rightarrow \; \; \left\{
    \begin{array}{ll}
        \boldsymbol{\mathcal{E}}^0 = \textbf{a}\textbf{N} + \textbf{b}\textbf{M} \\
        \boldsymbol{\mathcal{K}} = \textbf{d}\textbf{M} + \textbf{b}^T\textbf{N}
    \end{array}
\right.
\end{equation*}
</P>

<p style="text-align: justify;">
9. Calculer les comportements de souplesse pour les stratifiés de type Winckler en fonction de l’angle $\theta$. Quelle est la forme de ces comportements ?
</P><br>  

<p style="text-align: justify;">
Comment répond le stratifié si sollicité en traction pure selon l’axe x ? Est-ce que le maximum de cette réponse est obtenu pour la même valeur de $\theta = \theta_{max}$ ?
</P><br>  

<p style="text-align: justify;">
Pour calculer le comportements de souplesse en fonction de l'angle $\theta$ on utilise $\textbf{Matlab}$ pour pouvoir faire les calculs de a,b et déterminé à la question précédentes. En effet les calculs de a, b et d en fonction de $\theta$ sont grand, nous donnons alors ici leurs formes (voir le fichier joint $\textbf{Analytique.m}$ pour les expressions explicite de ces tenseurs);
</P><br>  
 
<p style="text-align: justify;">
\begin{equation*}a,d = 
\begin{bmatrix}
⚪ & ⚪ &  0\\
⚪  & ⚪ & 0\\
0 & 0 & ⚪
\end{bmatrix} \quad et \quad b = \begin{bmatrix}
0 & 0 & ⚪\\
0 & 0 & ⚪\\
⚪ & ⚪ & 0
\end{bmatrix}
\end{equation*}
</P>

avec,

<p style="text-align: justify;">
\begin{align*}
&a_{11} = a_{22} \quad et \quad a_{12}=a_{21}.\\ 
&d_{12}=d_{21}\\
&b_{16} = b_{26}  \quad et \quad  b_{61} \neq b_{62}
\end{align*}
</P>

Étudions alors la réponse à une traction simple:

<div style=" padding: 10px; border-radius: 5px; background-color: #080a10; color: #fff;">

```matlab
% Calcul de a b et d

A1=inv(A); D1=inv(D);

a=inv(A-(B*(D1*B)));
d=inv(D-(B*(A1*B)));
b1=-a*(B*D1);
b2=-d*(B*A1);

% Souplesse
S = [a(1,1), a(1,2), a(1,3), b1(1,1), b1(1,2) , b1(1,3)
     a(2,1), a(2,2), a(2,3), b1(2,1), b1(2,2) , b1(2,3)
     a(3,1), a(3,2), a(3,3), b1(3,1), b1(3,2) , b1(3,3)
     b2(1,1), b2(2,1), b2(3,1), d(1,1), d(1,2) , d(1,3)
     b2(1,2), b2(2,2), b2(3,2), d(2,1), d(2,2) , d(2,3)
     b2(1,3), b2(2,3), b2(3,3), d(3,1), d(3,2) , d(3,3) ];

% Application traction simple
F = [Nx, 0, 0, 0, 0, 0];
F=F';
EK = S*F;
</div>
```

<p style="text-align: justify;">
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
</P><br>  

<p style="text-align: justify;">
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
</P>

<p style="text-align: justify;">
On obtient pour une traction uniaxial des effets de déformations menbranaire, mais aussi un effet de courbure, due au terme $b_{16}$.
Ce couplage varie en fonction de $\theta$ et puisque cette expression comporte des termes en $cos4\theta$ et $sin4\theta$ on peut imaginer une valeur de $\theta_{max}$ proche de $\frac{\pi}{16}$, en effet d'après la courbe de suivante de $b_{16}$ en fonction de $\theta$ :
</P><br>  

<div style="text-align: center;">

  <figure style="display: inline-block;">
    <img src="/assets/projet3/p3f6.png" alt="Variation de la composante de couplage " width="600"/>
    <figcaption>Figure 6 : Variation de la composante de couplage $b_{16}$</figcaption>
  </figure>
</div><br>

<p style="text-align: justify;">
on vérifie que $\theta_{max} = \frac{\pi}{18}$ (cf. \textit{Application.m})
</P>

<p style="text-align: justify;">
$10.$ On considère une pale de largeur a et longueur $L$ fabriquée en composite stratifié selon la séquence $(1)$ avec $\theta = \theta_{max}$. 
</P><br>

<p style="text-align: justify;">
La pale est disposée dans le plan $O_{xy}$ et en rotation autour d’un axe vertical $O_z$ à une vitesse angulaire constante $\omega$ : soit $L >> a$, de telle manière qu'on puisse confondre la direction axiale et celle radiale. Dans ce cas, on peut dire donc que chaque point de la pale est soumis à un effort de traction $N_x = \rho . h . a . \omega . 2x$ ($x$ dans $[0 , L]$). Calculer la réponse de la pale à la sollicitation centrifuge (seule la composante $N_x$ est non nulle) et en particulier l’angle de torsion en bout de pale. Confirmer ce résultat par un calcul éléments finis sous Abaqus.
</P><br>

<p style="text-align: justify;">
Pour effectuer la modélisation d'une pale en traction-compression nous avons considéré une pale d'une lonqueur de $6m$, de l'argeur $50cm$ et d'une épaisseur de $8mm$, de masse volumique $1600 Kg/m^3$ et d'une vitesse de rotation de $157rad/s$. Cette pale est constitué d'un stratifié composé de huit couches unidirectionnelles en carbone/epoxyde d'épaisseur h. Détaillons les étapes de la modélisation (cf. $\textit{Winckler.cae}$):
</P><br>  

<p style="text-align: justify;">

<li> $\texttt{Part : 3D Shell Planar}$
<li> $\texttt{Geometrie : }$
    <div style="text-align: center;">
    <figure style="display: inline-block;">
    <img src="/assets/projet3/p3f7.png" alt=" geometrie " width="500"/>
    </figure>
    </div>
<li> $\texttt{Property : Elastic - Type : Lamina}$

<div style="overflow-x: auto; border-radius: 10px; background-color: #080a10; padding: 10px;">
<table style="border-collapse: collapse; width: 100%; color: white;">
  <thead>
    <tr style="background-color: #333;">
      <th style="border: 1px solid white; padding: 10px; border-top-left-radius: 10px;">Modules élastique</th>
      <th style="border: 1px solid white; padding: 10px;">$E_1$ (GPa)</th>
      <th style="border: 1px solid white; padding: 10px;">$E_2$ (GPa)</th>
      <th style="border: 1px solid white; padding: 10px;">$G_{12}$ (GPa)</th>
      <th style="border: 1px solid white; padding: 10px; border-top-right-radius: 10px;">$\nu_{12}$</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td style="border: 1px solid white; padding: 10px; border-bottom-left-radius: 10px;"></td>
      <td style="border: 1px solid white; padding: 10px;">181</td>
      <td style="border: 1px solid white; padding: 10px;">10,3</td>
      <td style="border: 1px solid white; padding: 10px;">7,17</td>
      <td style="border: 1px solid white; padding: 10px; border-bottom-right-radius: 10px;">0,28</td>
    </tr>
  </tbody>
</table>
</div>


<li> $\texttt{Composit Layup : Conventional Shell}$
    <div style="text-align: center;">
    <figure style="display: inline-block;">
    <img src="/assets/projet3/p3f8.png" alt=" step and load " width="700"/>
    </figure>
    </div>
<li> $\texttt{Step : Static, General}$
<li> $\texttt{Load : Body force, Analytical-Field:
   Nx = X*94651.6, Component 1}$
<li> $\texttt{Boundary condition : ENCASTRE}$
    <div style="text-align: center;">
    <figure style="display: inline-block;">
    <img src="/assets/projet3/p3f9.png" alt=" condition " width="500"/>
    </figure>
    </div>
    
<li> $\texttt{Mesh : Family Shell, S4R : finite menbrane strains}$
    <div style="text-align: center;">
    <figure style="display: inline-block;">
    <img src="/assets/projet3/p3f10.png" alt=" mesh " width="500"/>
    </figure>
    </div>
      
<li> $\texttt{Job : UR1}$
    <div style="text-align: center;">
    <figure style="display: inline-block;">
    <img src="/assets/projet3/p3f11.png" alt=" job " width="500"/>
    <img src="/assets/projet3/p3f12.png" alt=" jauge " width="100"/>
    </figure>
    </div>
<li> $\texttt{Plot}$
    <div style="text-align: center;">
    <figure style="display: inline-block;">
    <img src="/assets/projet3/p3f13.png" alt=" jauge " width="600"/>
    </figure>
    </div>
    
    
<li> $\texttt{Path dans l'épaisseur en bout de pale :}$
    <div style="text-align: center;">
    <figure style="display: inline-block;">
    <img src="/assets/projet3/p3f14.png" alt=" RX2 " width="500"/>
    </figure>
    </div>

<li> $\texttt{Path dans l'épaisseur en bout de pale :}$
	<div style="text-align: center;">
    <figure style="display: inline-block;">
    <img src="/assets/projet3/p3f15.png" alt=" path2 " width="600"/>
    </figure>
    </div>    
</P><br>

<p style="text-align: justify;">
L'angle de tortion $\alpha$ obtenue en bout de pale est de $\textbf{13.94 degré}$.
</P><br>

<p style="text-align: justify;">
Théoriquement, pour des plaques minces la théorie cinématique de $\textit{Love-Kirchoff}$ permet d'écrire: 
</P>

<p style="text-align: justify;">
\begin{equation*}
\mathcal{K}_6 = 2\frac{\partial^2 W_0(x,y)}{\partial x \partial y}
\end{equation*}
</P><br>

<p style="text-align: justify;">
$\Rightarrow$
\begin{equation*}
    W_0(x,y) = \frac{1}{2}\int_0^L\int_0^l \mathcal{K}_6 dxdy = \frac{1}{2}K_6*L*l 
\end{equation*}
</P><br>

<p style="text-align: justify;">
Donc le déplacement selon $z$  vaut $W_0 = 0.8421 \, m$ et le coté adjacent vaut $adj = \frac{l}{2} = 0.25\,m$ et donc $\alpha=atan(\frac{W_0}{adj}) = 73.47\, \text{degré}$.
</P><br>

<p style="text-align: justify;">
Cette valeur ne correspond pas à la valeur trouvé numériquement, elle est 5 fois trop grandes, et je ne trouve pas où mon raisonnement est erroné (cf $\textit{Application.m}$)
</P>

### IV. Conclusion

<p style="text-align: justify;">
Nous avons pu étudier lors de ce projet la conception de structures stratifiés couplés élastiquement afin de produire des changements de forme lors d’un chargement mécanique. Le couplage traction-torsion peut être exploité pour construire des aubes ou pales à pas adaptatif, la torsion intervenant lors de la rotation des pales/aubes par effet de la sollicitation produite par effet centrifuge. Ce type de composite peut-être obtenu avec de stratifié à tissu équilibré, mais aussi en donnant une orientation bien spécifique à la séquence des couches. Nous avons alors étudié celle proposé par $\textit{Winckler}$ en 1985 :
</P><br>

<div style="text-align: center;">
  <figure style="display: inline-block;">
    <img src="/assets/projet3/p3f16.png" alt="Winckler " width="600"/>
    <figcaption>Figure 16 : Séquence de Winckler</figcaption>
  </figure>
</div><br>

<p style="text-align: justify;">
Cette séquence permet d'obtenir des stratifiés ayant un couplage thermo-élastiquement nul $\textbf{V = O}$ et $\textbf{U}$ sphérique, i.e thermo-élastiquement stables. Mais aussi la séquence de $\textit{Winckler}$ permet de fournir en plus la possibilité de moduler l'intensité du couplage élastique $\textbf{B}$ à l'aide du paramètre géométrique $\theta$ désiré.
</P>
