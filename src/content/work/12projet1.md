---
title: Projet 1
publishDate: 2020-03-04 00:00:00
img: /assets/stock-3.jpg
img_alt: Pearls of silky soft white cotton, bubble up under vibrant lighting
description: |
  Explosion dynamique d'une barre de fer console.
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

### Introduction
<p style="text-align: justify;">
Dans le cadre de problèmes de mécanique vibratoire avancée, l'étude des structures doit se faire en tenant compte des forces inertielles d'entraînement afin d'obtenir les équations dynamiques du problème traité.
</p>

<p style="text-align: justify;">
Dans notre étude, nous nous concentrons sur la propagation d'une onde longitudinale au sein d'une barre déclenchée par une explosion. 
</p>

<p style="text-align: justify;">
L'objectif est de comparer différentes méthodes de schéma d'intégration temporelle pour évaluer leur précision en fonction du pas de temps. Nous comparerons également la procédure PASAPAS de Cast3m avec celle utilisée dans Abaqus.
</p>

### Description du problème
<p style="text-align: justify;">
Considèrons la barre de fer illustrée en Figure 1 de longueur \( L \) et soumise à une pression uniforme \( p(t) \) à son extrémité libre en \( x = 0 \). À l'autre extrémité, en \( x = L \), la barre est encastrée. La forme de la charge en fonction du temps simule l'effet d'une explosion à proximité de l'extrémité libre, tel que, la barre est au repos à l'instant \( t = 0 \) et soumis l'explosion d'une durée \( T_p \) de \( 3{,}88 \times 10^{-5} \) secondes (cf. figure 2.) :
</P>

<p style="text-align: justify;">
\[ p(t) = 
\begin{cases} 
0.1 \, \text{MPa} & \text{pour } &  0 < t < T_p \\
\\
0 \, \text{MPa} & \text{pour } &  t > T_p 
\end{cases} \]
</P>

<p style="text-align: justify;">
Nous supposons dans cette étude, la barre élastique, isotrope et homogène de module de Young \( E \) et de coefficient de Poisson \( \nu \). La masse volumique \( \rho \), la section \( S \) de la barre seront prises comme constantes :
</p>
<p style="text-align: justify;">
<li> \( E = 207 \times 10^9 \, \text{Pa} \)
</p>
<p style="text-align: justify;">
<li> \( \nu = 0.3 \)
</p>
<p style="text-align: justify;">
<li> \( \rho = 7800 \, \text{kg/m}^3 \)
</p>
<p style="text-align: justify;">
<li> \( S = 0.2 \, \text{m} \times 0.2 \, \text{m} \)
</p>
<p style="text-align: justify;">
<li> \( L = 1 \, \text{m} \)
</p>

<p style="text-align: justify;">
\(   \)
</P>

<p style="text-align: justify;">
Pour la résolution de ce problème, nous avons utilisé des méthodes directes d'intégration de la dynamique, en suivant un plan méthodique :</P>

<p style="text-align: justify;">
\(   \)
</P>

<p style="text-align: justify;">
\(1. \)  Formulation du problème
</P> 

<p style="text-align: justify;">
\(2. \)  Analyse préliminaire
</P> 

<p style="text-align: justify;">
\(3. \)  Algorithme de solution et implémentation dans Castem
</P> 

<p style="text-align: justify;">
\(4. \)  Résolution numérique, vérification et analyse
</P> 

### 1. Formulation du problème  

<p style="text-align: justify;">
Nous avons d'abord considéré un modèle d'une barre de traction en 1D, élastique linéaire, pour décrire la structure de la formulation.
</P> 
<p style="text-align: justify;">
\(   \)

</P>
<p style="text-align: justify;">
Plan du sous-problème :
</P> 
<p style="text-align: justify;">
\(   \)
</P> 
<p style="text-align: justify;">
\(\text{a.} \) Conditions initiales et aux limites
</P> 
<p style="text-align: justify;"> 
\(\text{b.} \) Formulation forte
</P>
<p style="text-align: justify;">
\(\text{c.} \) Formulation faible
</P>
<p style="text-align: justify;">
\(\text{d.} \) Formulation semi-discrète
</P> 

#### 1.a Conditions initiales et aux limites  :

<p style="text-align: justify;">
- La barre est encastrée à son extrémité droite, en $x = L$ : : 
</P> 
$$ u(L,t) = 0 \; \text{(H_1)} $$

<p style="text-align: justify;">
- La barre est soumise à un chargement \( p(t) \) à son extrémité gauche en $x = 0$ : :
</P> 
 $$ N(0,t) = p(t)S \) \(\text{(H_2)} \) $$

<p style="text-align: justify;">
- On choisit un état de repos initial pour la barre  
  (H_3) : En déplacement  \( u(x,0) = 0 \)  
  (H_4) : En vitesse : \( \frac{\partial u(x,0)}{\partial t} = 0 \) 
</P> 
#### 1.b Formulation forte

##### (i) Équation d'équilibre

\[
\frac{\partial N(x,t)}{\partial x} = \rho S \frac{\partial^2 u(x,t)}{\partial t^2}
\]

##### (ii) Loi de comportement
\[
N(x,t) = ES \epsilon(x,t)
\]

##### (iii) Équation de compatibilité
\[
\epsilon(x,t) = \frac{\partial u(x,t)}{\partial x} 
\]

En injectant (\(\epsilon(x,t) = \frac{\partial u(x,t)}{\partial x}\)) dans (\(N(x,t) = ES \epsilon(x,t)\)) et ceci dans (\(\frac{\partial N(x,t)}{\partial x} = \rho S \frac{\partial^2 u(x,t)}{\partial t^2}\)) nous obtenons l'équation de d'Alembert qui est une équation différentielle hyperbolique : 

\[
\frac{\partial^2 u(x,t)}{\partial x^2} - c_d^2 \frac{\partial^2 u(x,t)}{\partial t^2} = 0
\]

avec, 
\[
c_d = \sqrt{\frac{E}{\rho}}
\]

la célérité des ondes longitudinales.

#### 1.c Formulation faible

Établissons la formulation faible du problème en introduisant d'abord les ensembles de champs admissibles \(\mathcal{U}_{ad}\) et admissible à zéro \(\mathcal{U}_{ad}^o\) par rapport aux données du problème considéré:

\[
\mathcal{U}_{ad} = \{ u(x) \in H^1(\Omega) \; \vert \; u(L) = 0 \}
\]

L'espace cinématiquement admissible à zéro \(\mathcal{U}_{ad}^o\) est ici équivalent à \(\mathcal{U}_{ad}\):

\[
\mathcal{U}_{ad}^o = \{ \hat u(x) \in H^1(\Omega) \; \vert \; \hat u(L) = 0 \}
\]

#### Remarque

Les fonctions \( u \) et \( \hat{u} \) doivent être suffisamment régulières puisque la formulation locale requiert a priori un champ cinématique admissible deux fois continûment différentiable par morceaux. Ainsi, l'espace de Sobolev \( H^1 \) (\( \Omega \)) dont les fonctions sont de carré intégrable suffit pour vérifier cette condition.

Il s'agit maintenant d'exprimer l'équation du mouvement par une forme intégrale équivalente par dualisation, c'est-à-dire par multiplication par un champ de déplacement virtuel \( \hat{u} \) cinématiquement admissible et par intégration sur le support géométrique \( \Omega=[0,L] \).

\[
\int_{\Omega} \frac{\partial N(x,t)}{\partial x}\hat{u}(x) \,dx = \int_{\Omega}\rho S \frac{\partial^2 u(x,t)}{\partial t^2}\hat{u}(x) \,dx
\]

On utilise ensuite l'Intégration Par Parties, ce qui nous ramène à:

\[
\big[ N(x,t) \hat{u} \big]_0^L - \int_{\Omega} N(x,t) \frac{d \hat{u}(x)}{d x} \,dx = \int_{\Omega}\rho S \frac{\partial^2 u(x,t)}{\partial t^2}\hat{u}(x) \,dx
\]

Puis à l'aide des conditions aux limites (H_1) et (H_2), on a :

\[
-p(t)S\hat{u}(0) - \int_{\Omega} N(x,t) \frac{d \hat{u}(x)}{d x}\,dx = \int_{\Omega}\rho S \frac{\partial^2 u(x,t)}{\partial t^2}\hat{u}(x) \, dx
\]

Et en utilisant la loi de comportement, il vient :

\[
\int_{\Omega} ES\frac{\partial u(x,t)}{\partial x} \frac{d \hat{u}(x)}{d x}\,dx + \int_{\Omega}\rho S \frac{\partial^2 u(x,t)}{\partial t^2}\hat{u}(x) \, dx = -p(t)S\hat{u}(0) 
\]

Finalement, la formulation faible du problème s'écrit :

Trouver \( u \) \(\in\) \(\mathcal{U}_{ad}\) tel que:

\[
a_k(u,\hat{u}) + a_m\left( \frac{\partial^2 u}{\partial t^2},\hat{u}\right)= l(\hat{u}) \qquad \forall \hat{u} \, \in \, \mathcal{U}_{ad}^o
\]

avec

\[
a_k(u,\hat{u}) = ES\int_{\Omega} \frac{\partial u(x,t)}{\partial x} \frac{d \hat{u}(x)}{d x}\,dx
\]

\[
a_m\left( \frac{\partial^2 u}{\partial t^2},\hat{u}\right) = \rho S\int_{\Omega} \frac{\partial^2 u(x,t)}{\partial t^2}\hat{u}(x) \, dx
\]

qui représentent respectivement l'opérateur de raideur modale et l'opérateur de masse modale qui sont tous deux de forme bilinéaire, continue et coercive (définie positive);

et

\[
l(\hat{u}) = -p(t)S\hat{u}(0) 
\]

qui est une forme linéaire et continue.

#### 1.d Formulation semi-discrète

La solution analytique de l'équation est en général inaccessible puisque les ensembles de champs admissibles sont de dimension infinie. On est donc conduit à chercher une solution approchée par une méthode numérique. On commence par séparer nos variables temporelles et spatiales, soit : 

\[
u(x,t) = w(x)d(t) 
\]


\noindent
La méthode de Rayleigh-Ritz suppose que la solution peut être recherchée comme une combinaison linéaire de fonctions de forme correctement choisies dans l'espace cinématiquement admissible. On pose ainsi :
\begin{equation}
w(x)=\sum\limits_{i=1}^nN_i(x)
\end{equation}

\noindent
où les $N_i(x)$ sont les modes propres.\\

\noindent
En injectant dans $(16)$, on obtient la solution semi-discrétisée en espace :
\begin{equation}
u(x,t)=\sum\limits_{i=1}^nN_i(x)d(t) 
\label{u}
\end{equation}

\noindent
Ainsi, on discrétise \textbf{la formulation faible du problème d'ondes continues} en introduisant respectivement un sous-espace des fonctions essai $u$ et test $\hat{u}$ : $\mathcal{U}_{ad}^h$ $\subset$ $\mathcal{U}_{ad}$

\noindent
Et en posant :
\begin{equation}
\hat{u}=N_j(x)  \qquad j=[1,..,n] 
\label{v}
\end{equation}

\noindent
On obtient le nouveau problème discrétisé en injectant (\ref{u}) et (\ref{v}) dans (\ref{faible}) : 

Trouver $\hat{u}$ $\in$ $\mathcal{U}_{ad}$ tel que:
\begin{equation}
a_k(N_i,N_j) + a_m(N_i,N_j)= l(N_j) \qquad \forall \hat{u} \, \in \, \mathcal{U}_{ad} 
\label{faible_N}
\end{equation}

\noindent
avec
\begin{equation}
a_k(N_i,N_j) = d(t)\underbrace{ES\int_{\Omega} \frac{d N_i(x)}{d x} \frac{d N_j(x)}{d x}\,dx}_{K_{ij}} 
\end{equation}

\begin{equation}
a_m(N_i,N_j) =  \ddot{d}(t)\underbrace{\rho S \int_{\Omega} N_i(x) N_j(x) \, dx}_{M_{ij}}
\end{equation}

\noindent
et,
\begin{equation}
l(N_j) = \underbrace{-p(t)N_j(0)}_{F_j(t)} 
\end{equation}

\noindent
Que l'on peut écrire sous la forme suivante :
\begin{equation}
M_{ij}\ddot{d}(t) + K_{ij}d(t) = F_j(t)
\end{equation}

\noindent
%En utilisant l'orthogonalité des modes:
%Pour i$\neq$j : ... à faire pour trouver le problème semi-discrétisé sous forme modale.

\section*{Analyse préliminaire}
\label{analyse}
\noindent Il est judicieux de commencer à calculer les différentes caractéristiques fondamentales du système à l'aide d'une analyse élémentaire de propagation d'ondes afin de prédire les ordres de grandeur de la solution cherchée.

\begin{enumerate}
\item
Calculons la vitesse $c_d$ de propagation  des ondes longitudinales dans la barre d'aprés (\ref{c}) et des donnés suivantes :
\end{enumerate}

\begin{align*} 
&E = 207.10^9 \, Pa \\
&\rho = 7800 \, kg/m^3 \\
&c_d = \sqrt{\frac{E}{\rho}} = 5.15.10^3 \, m.s^{-1} 
\end{align*}
\begin{enumerate}[start = 2]
\item Calculons le temps $T_0$ nécessaire à une onde pour parcourir toute la longueur de la barre :
\end{enumerate}

\begin{align*} 
&L = 1 \, m \\
&T_0 = \frac{L}{c_d} = 1.94.10^{-4} \, s
\end{align*}
\begin{enumerate}[start =  3]
\item L'explosion génère un impact en $x=0$ provoquant une onde de pression très brève entre les instants $t=0$ et $t=T_p$ se propageant dans la barre à la vitesse $c_d$. Calculons la distance $L_p$ parcourue par cette onde à t = $T_p$ :
\end{enumerate}

\begin{align*} 
&T_p = 3.88.10^{-5} \, s\\
&L_p = c_d \, T_p = 0.1999 \, m
\end{align*}
\begin{enumerate}[start = 4]
\item Supposons d'adopter pour le système un modèle de barre 1D et une discrétisation en espace par des éléments finis linéaires.
\end{enumerate}

\begin{enumerate}[label=(\alph*).]
\item Donnons l'expression des matrices élémentaires $M_e$ et $K_e$,
respectivement de masse et de rigidité.
\end{enumerate}

\begin{equation}
([K_e] + [M_e])\Bigg( \begin{matrix}
UX_1  \\
UX_2 
\end{matrix}\Bigg) = \begin{pmatrix}
0 \\
0 
\end{pmatrix}
\end{equation}
\noindent
Le système à résoudre s'écrit alors, avec une matrice de masse consistante:
\begin{equation}
\Bigg(\underbrace{\frac{ES}{l_e}\begin{bmatrix}
1 &-1 \\
-1 & 1 
\end{bmatrix}}_{K_e} \, \, \underbrace{- \, \frac{\omega^2 ml_e}{6}\begin{bmatrix}
2 & 1 \\
1 & 2 
\end{bmatrix}}_{M_e}\Bigg)\Bigg( \begin{matrix}
UX_1  \\
UX_2 
\end{matrix}\Bigg) = \begin{pmatrix}
0 \\
0 
\end{pmatrix}
\end{equation}
\noindent
Dans le cas d'une matrice de masse condensée (\textit{lumped mass}), le système s'écrit :
\begin{equation}
\Bigg(\underbrace{\frac{ES}{l_e}\begin{bmatrix}
1 &-1 \\
-1 & 1 
\end{bmatrix}}_{K_e} \, \, \underbrace{- \, \frac{\omega^2 ml_e}{2}\begin{bmatrix}
1 & 0 \\
0 & 1 
\end{bmatrix}}_{M_e}\Bigg)\Bigg( \begin{matrix}
UX_1  \\
UX_2 
\end{matrix}\Bigg) = \begin{pmatrix}
0 \\
0 
\end{pmatrix}
\end{equation}

\begin{enumerate}[label=(\alph*)., start= 2]
%[label=\Alph*,start = b]
\item Calculons la pulsation la plus élevée d'un élément de longueur $l_e$ lorsque l'on utilise une matrice de masse consistante
\end{enumerate}
\begin{equation}
det\Bigg(\frac{6E}{\rho l_e^2}\begin{bmatrix}
1 & -1 \\
-1 & 1 
\end{bmatrix} - \omega^2\begin{bmatrix}
2 & 1 \\
1 & 2 
\end{bmatrix}\Bigg) = 0
\end{equation}

\begin{equation}
det\Bigg(\begin{matrix}
\frac{6E}{\rho l_e^2} - 2\omega^2 & -\frac{6E}{\rho l_e^2} - \omega^2 \\
-\frac{6E}{\rho l_e^2} - \omega^2 & \frac{6E}{\rho l_e^2} - 2\omega^2
\end{matrix}\Bigg) = \omega^4 -\frac{12E\omega^2}{\rho l_e^2} = 0
\end{equation}

\begin{equation}
\Rightarrow \omega_{1,2} = \pm \sqrt{\frac{12E}{\rho l_e^2}}
\end{equation}
\noindent
Nous gardons les pulsations positives
\begin{equation}
\omega = \frac{2\sqrt{3}}{l_e}\sqrt{\frac{E}{\rho}} = 2\sqrt{3} \, \frac{c_v}{l_e}\\
\end{equation}
\noindent
Calculons la pulsation la plus élevée d'un élément de longueur $l_e$ lorsque l'on utilise une matrice de masse condensée
\begin{equation}
det\Bigg(\frac{2E}{\rho l_e^2}\begin{bmatrix}
1 & -1 \\
-1 & 1 
\end{bmatrix} - \omega^2\begin{bmatrix}
1 & 0 \\
0 & 1 
\end{bmatrix}\Bigg) = 0
\end{equation}

\begin{equation}
det\Bigg(\begin{matrix}
\frac{2E}{\rho l_e^2} - \omega^2 & -\frac{2E}{\rho l_e^2} \\
-\frac{2E}{\rho l_e^2} & \frac{2E}{\rho l_e^2} - \omega^2
\end{matrix}\Bigg) = \omega^4 -\frac{4E\omega^2}{\rho l_e^2} = 0
\end{equation}

\begin{equation}
\Rightarrow \omega_{1,2} = \pm \sqrt{\frac{4E}{\rho l_e^2}}
\end{equation}
\noindent
Nous gardons les pulsations positives

\begin{equation}
\omega = \frac{2}{l_e}\sqrt{\frac{E}{\rho}} = 2 \, \frac{c_v}{l_e}
\end{equation}
\begin{enumerate}[label=(\alph*)., start= 3]

\item Déterminons le pas de temps $\Delta t_0$ pour respecter la condition de stabilité (dite de Courant) lorsque la discrétisation en espace est faite par des éléments de taille uniforme $l_e$, et en temps par la différence finie centrée. Déterminons le $\Delta t_0$ pour les deux cas suivants :
\end{enumerate}

\begin{enumerate}[label=i.]
\item Éléments fins avec matrice de masse consistante
\end{enumerate}

\begin{equation}
\Delta t_0 \leq \frac{\Omega_c}{\omega_h} = \frac{\Omega_c \, l_e}{2 \, c_v}
\end{equation}

Dans ce cas : $\Omega_c = 2$ et $\beta = 0$ , $\gamma = \frac{1}{2}$.
\begin{enumerate}[label=ii.]
\item Éléments finis avec matrice de masse condensée
\end{enumerate}
\begin{equation}
\Delta t_0 \leq \frac{\Omega_c}{\omega_h} = \frac{\Omega_c \, l_e}{2\sqrt{3} \, c_v}
\end{equation}

Dans ce cas : $\Omega_c = 2\sqrt{3}$ et $\beta = \frac{1}{6}$, $\gamma = \frac{1}{2}$.


\section*{Algorithme de solution et implémentation dans Cast3M}

\begin{enumerate}
\item Donnons sur un schéma les éléments de modélisation nécessaires au calcul :\\
Comme évoqué dans la partie description du modèle\ref{description}, le problème étudié est schématisé comme suit,
\end{enumerate}

\begin{figure}[H] \centering 
\includegraphics[width=9cm]{./images/description.png}
\caption{Modélisation de la barre}
\end{figure}

qui en passant à la discrétisation spatiale nous donne le schéma suivant,

\begin{figure}[H] \centering 
\includegraphics[width=7cm]{./images/schema.jpg}
\caption{Modélisation discrète de la barre}
\end{figure}
\noindent
ainsi les paramètres d'entrée nécessaires au calcul sont alors implémentés au code, 
%gibiane en commentaire car consomme bcps de memoir a allouer
%\lstinputlisting[language=gibiane,firstline=40, lastline=53]{./gibiane/tp1.dgibi}
et à l'aide de l'étude préliminaire\ref{analyse} sur l'analyse des grandeurs caractéristiques du problème, on effectue la discrétisation temporelle.
%\lstinputlisting[language=gibiane,firstline=100, lastline=125]{./gibiane/tp1.dgibi}
\begin{enumerate}[start= 2]
\item Détaillons les conditions aux limites et initiales.\\
Comme il à été énoncé dans la description du problème\ref{description}, nous pouvons écrire les conditions aux limites et initiales.
\end{enumerate}
%\lstinputlisting[language=gibiane,firstline=76, lastline=77]{./gibiane/tp1.dgibi}
%\lstinputlisting[language=gibiane,firstline=128, lastline=130]{./gibiane/tp1.dgibi}
\begin{enumerate}[start= 3]
\item Les intégrations temporelles qui peuvent être envisagées sont les suivantes :
\end{enumerate}


\bigbreak
\noindent
\begin{tabular}{|l||c|c|c|c|c|}
\hline Schéma & $\gamma$ &  $\beta$ &  $\Omega_c$  & $\frac{\Delta T}{T}$ & Propriétés\\
\cline{1-6}  Newmark différence centrée & $\frac{1}{2}$ & 0 & 2 & $\frac{-w^2h^2}{24}$ & E $\&$ C-S \\
Newmark accélération moyenne & $\frac{1}{2}$ & $\frac{1}{4}$ & $\infty$ & $\frac{w^2h^2}{12}$  & I-S \\
Newmark accélération moyenne modifiée & $\frac{1}{2}$ + $\alpha$ & $\frac{(1+\alpha)^2}{4}$ &  $\infty$  & ($\frac{1}{12}+\frac{\alpha^2}{4})w^2h^2$ & I-S \\
HHT & $\gamma_0$ & $\beta_1$  & - & - & I-S\\
$\alpha$-généralisé $\rho$ = 0.905 & $\gamma_1$ & $\beta_1$ & - & - & I-S \\
$\alpha$-généralisé $\rho$ = 0.5 & $\gamma_1$ & $\beta_1$ & - & - & I-S \\ 
\hline 
\end{tabular}
\noindent
E : Explicite\\
I-S : Inconditionnellement Stable\\
C-S : Conditionnellement Stable\\
$\gamma_0$ = $\frac{1}{2}$ - $\alpha_f$\\
$\gamma_1$ = $\frac{1}{2}$ + $\alpha_m$ - $\alpha_f$ \\
$\beta_1$ = $\frac{1}{4}$ + $\frac{1}{2}$($\gamma_m$ - $\gamma_f$) 

\medbreak
\noindent
Les trois premières sont des méthodes faisant partie de la famille de Newmark tandis que les trois dernières en sont des méthodes dérivées. En outre, il est à souligner que la méthode HHT est un cas particulier de la méthode $\alpha$-généralisé lorsque $\alpha_m=0$. 
\begin{enumerate}[start = 4]
\item Listons les étapes à suivre dans Cast3M 
\end{enumerate}
\noindent
\texttt{OPTIONS}\\
\quad Tout d'abord on commence par compléter les options du modèle, à savoir la dimension de l'espace, le modèle de calcul, le choix de l'élément à fabriquer et la méthode de calcul des déformations.Puis, le choix de l'algorithme à utiliser, présentés à la question précédente, est laissé libre à l'utilisateur.
\noindent
\texttt{DONNÉES}\\
\quad On entre ensuite les données du problème, autrement dit la géométrie, les propriétés matérielles, le chargement - en particulier la durée de l'explosion $T_p$ - et le nombre d'éléments à discrétiser.
\noindent
\textbf{\texttt{MAILLAGE}}\\
\quad Ceci étant fait, le maillage est ainsi réalisé.
\noindent
\texttt{MODÈLE}\\
\quad Le modèle du problème est alors construit à partir des propriétés géométriques et matérielles définies précédemment. Les matrices de raideur $K$ et de masse $M$ sont calculées auxquelles les conditions aux limites sont appliquées.
\noindent
\texttt{CHARGEMENT}\\
\quad De ce fait, on représente alors le step de pression avant et après l'explosion, ainsi que le chargement impliqué de ladite explosion.
\noindent
\texttt{DICRÉTISATION TEMPORELLE}\\
\quad Nous procédons a posteriori à la discrétisation temporelle en veillant à sélectionner la vitesse de propagation $c_d$, une matrice de masse consistante ou lumpée ainsi que le pas de temps critique et le temps $T_0$ mis par l'onde pour traverser la barre. Puis, on choisit un pas de temps et on écrit les conditions initiales.
\noindent
\texttt{RÉSOLUTION}
\quad À présent, nous pouvons effectuer la résolution de ce problème.La procédure \texttt{DYNAMIC} permet de réaliser un calcul dynamique pas à pas pour les algorithmes vus précédemment. Ceci est fait via les données du problème, à savoir les déplacement et vitesse initiaux, les matrices de rigidité, de masse, et d'amortissement. On va, de plus, sauvegarder les instants de calcul. Suivant le choix du schéma d'intégration, on va entrer les paramètres $\alpha$, $\gamma$, $\alpha$ ou $\rho^\infty$.
\noindent
\texttt{POST-TRAITEMENT}
\quad Finalement, on post-traite les résultats intéressants, que l'on décrira à la question suivante, et on les représente graphiquement.

\begin{enumerate}[start= 5]
\item Les grandeurs pertinentes à post-traiter sont:
\end{enumerate}


\noindent
\ding{226} le déplacement longitudinal\\
\ding{226} le déplacement transversal\\
\ding{226} la vitesse longitudinale\\
\ding{226} l'effort normal\\
\ding{226} le travail des efforts extérieurs\\
\ding{226} le travail des efforts intérieurs\\
\ding{226} le bilan d'énergie\\


\begin{enumerate}[start=6 ]
\item Écrivons le programme gibiane pour Cast3M:\\
(cf. fichier .dgibi joint )
\end{enumerate}

\newpage
\section*{Résolution numérique, vérification et analyse}

\begin{enumerate}
\item Vérifions les résultats obtenus et en particulier l'effet du pas de temps.
\end{enumerate}
\quad
Le programme a été lancé pour les différents algorithmes et en faisant varier le pas de temps tout en veillant à ce qu'il soit maintenu inférieur au pas de temps critique déterminé via la condition de stabilité.\\
Leurs représentations sont jointes dans les fichiers \textit{postscript}.

\noindent
Notre première étude portera sur un pas de temps égal au pas de temps critique.\\
Pour le premier algorithme, le schéma de Newmark différence centrée, les grandeurs tracées présentent de fortes variations ou alors un seul pic tout le long de la barre. Ceci s'explique par le fait que ce schéma est conditionnellement stable, il faut donc un pas de temps inférieur au pas de temps critique pour ce schéma afin de visualiser une solution correcte.\\
Les autres schémas présentent des résultats plus cohérents, et en particulier le schéma de Newmark accélération moyenne modifiée, étant donné qu'ils sont tous inconditionnellement stables.
\noindent
Des pas de temps cinq et dix fois inférieurs ont aussi été traités mais l'effet du pas temps sera particulièrement observable pour le premier schéma dans la mesure où il s'agit du seul schéma conditionnellement stable. En effet, nous constatons que ce schéma présente une solution convenable uniquement dès lors que le pas de temps est près de deux fois inférieur au pas de temps critique.

%à faire ensemble je passe au tp batiment -OK
\begin{enumerate}[start = 2]
\item Analysons la solution du problème en la comparant avec l'analyse de l'exemple du chapitre 9.4 du livre \textit{Getting Started with Abaqus: Interactive Edition}. La solution donne la contrainte longitudinale selon $\b{$e$}_3$ soit (S33) en trois points différents le long de la barre ($0.25 m$, $0.5 m$, et $0.75 m$)
\end{enumerate}
\vspace*{-0.7cm}
 \flushleft
\begin{figure}[H] \centering
\includegraphics[width=9cm]{./images/abaqus.png}
\vspace*{-0.7cm}
\caption{Solution en contraintes sur Abaqus}
\end{figure}
\begin{figure}[H] \centering
\noindent
Sur Cast3m, nous avons pris les mêmes points pour comparer la solution en contraintes.
\includegraphics[width=9cm]{./images/castem.png}
\vspace*{-0.5cm}
\caption{Solution en contraintes sur Cast3m}
\end{figure}

\underline{Conclusion} : Nous obtenons une solution quasi-identique pour le même problème résolu par un logiciel de programmation (Cast3m) et un logiciel de simulation (Abaqus), ce qui permet de valider nos résultats.


\begin{enumerate}[start = 3]

\item Comparons les avantages et les inconvénients des méthodes explicite et implicite.
\end{enumerate}

La résolution d’un problème de dynamique s'interroge tout d'abord sur la nature du problème posé afin de choisir le schéma d’intégration temporelle, qui définit la discrétisation en temps de l’équation d’équilibre globale.

Les schémas sont classés en deux grandes familles: les schémas explicites et les schémas implicites. 

Les premiers utilisent un pas de temps très réduit et ont pour avantage de permettre la résolution de l’équilibre en dynamique rapide, pour des phénomènes très violents (chocs, explosions, ...). L'avantage dudit schéma réside dans l'adaptation du phénomène à observer qui impose déjà un pas de temps très petit (souvent de l’ordre de la microseconde). L’inconvénient du schéma explicite est qu’il existe un pas de temps critique à ne pas dépasser et les temps de calcul sont très longs.
%de n’avoir à inverser qu’une matrice de masse (lumpée). Il n’y a, de plus, pas de boucle de convergence car la résolution est directe. 
%Enfin, la qualité de la solution obtenue n’est pas garantie par un critère de résidu en équilibre. 
D'autre part, les phénomènes rapides induisent de fortes non-linéarités sur un temps très court, le schéma explicite est donc beaucoup plus efficace que le schéma implicite.

En revanche, les schémas implicites sont utiles en dynamique lente (vibrations, ébranlements des structures), on est en présence d’un comportement linéaire des matériaux, pour une durée des phénomènes s’étalant sur plusieurs secondes. Les schémas implicites ont donc cette fois-ci tout leur intérêt car ils permettent d’utiliser des pas de temps très grands et sont alors très rapides en temps de calcul. La solution numérique reste stable et sa qualité est contrôlée par le critère de résidu en équilibre, tout comme en quasi-statique. Cependant, la convergence s’avère parfois très délicate en présence de fortes non-linéarités. 

%Résumé du texte

\begin{enumerate}[start = 4]
\item Réanalysons le problème dans le cas où la barre est libre à l'extrémité droite et non plus encastrée.
\end{enumerate}
Nous pouvons, par exemple, comparer le déplacement selon $x$ pour les deux cas.

\begin{figure}[H]
\begin{minipage}[t]{8cm}
\includegraphics[width=8.5cm]{./images/libre.png}
\caption{Déplacement $u_x$ : bord libre }
\end{minipage}
\begin{minipage}[t]{8cm}   
\caption{Déplacement $u_x$ : bord encastré}
\includegraphics[width=8.5cm]{./images/encastre.png}
\end{minipage}
\end{figure}
Nous remarquons que l'onde en $x = 0$ va s'annuler puis devenir négative dans le cas d'un encastrement; et s'additionner dans le cas d'un bord libre, car l'onde est dans ce cas non inversé comme schématisé ci-dessous:
\vspace{-1cm}
\begin{figure}[H] \centering 
\includegraphics[width=8.5cm]{./images/onde.png}
\vspace{-2.8cm}
\caption{Modélisation du comportement de l'onde}
\end{figure}

