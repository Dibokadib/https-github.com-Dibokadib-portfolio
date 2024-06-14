---
title: Projet 1
publishDate: 2020-03-04 00:00:00
img: /assets/stock-1.jpg
img_alt: Pearls of silky soft white cotton, bubble up under vibrant lighting
description: |
  Réponse d'une barre à une explosion.
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

### I. Introduction
<p style="text-align: justify;">
Dans le cadre de problèmes de mécanique vibratoire avancée, l'étude des structures doit se faire en tenant compte des forces inertielles d'entraînement afin d'obtenir les équations dynamiques du problème traité.
</p>

<p style="text-align: justify;">
Dans notre étude, nous nous concentrons sur la propagation d'une onde longitudinale au sein d'une barre déclenchée par une explosion. 
</p>

<p style="text-align: justify;">
L'objectif est de comparer différentes méthodes de schéma d'intégration temporelle pour évaluer leur précision en fonction du pas de temps. Nous comparerons également la procédure PASAPAS de Cast3m avec celle utilisée dans Abaqus.
</p>

<a id="sect2"></a>
### II. Description du problème
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
III. Formulation du problème
</P> 

<p style="text-align: justify;">
IV.  Analyse préliminaire
</P> 

<p style="text-align: justify;">
V. Algorithme de solution et implémentation dans Castem
</P> 

<p style="text-align: justify;">
VI. Résolution numérique, vérification et analyse
</P> 

### III. Formulation du problème  

<p style="text-align: justify;">
Nous avons d'abord considéré un modèle d'une barre de traction en 1D, que nous appelerons dès à present une tige, élastique linéaire, pour décrire la structure de la formulation.
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
III.1 Conditions initiales et aux limites
</P> 
<p style="text-align: justify;"> 
III.2 Formulation forte
</P>
<p style="text-align: justify;">
III.3 Formulation faible
</P>
<p style="text-align: justify;">
III.4 Formulation semi-discrète
</P> 

#### 1. Conditions initiales et aux limites  :

<p style="text-align: justify;">
- La barre est encastrée à son extrémité droite, en $x = L$ :
</P>

<a id="eq1"></a>
$$ u(L,t) = 0  \tag{$C_1$}$$

<p style="text-align: justify;">
- La barre est soumise à un chargement \( p(t) \) à son extrémité gauche, en $x = 0$ :
</P> 

<a id="eq2"></a>
 $$ N(0,t) = p(t)S \tag{$C_2$}$$

<p style="text-align: justify;">
Supposons l'état de repos initial pour la barre, 
</P> 

en déplacement : 
<a id="eq3"></a>
$$ u(x,0) = 0, \tag{$C_3$} $$
et en vitesse :
<p style="text-align: justify;">
$$\frac{\partial u(x,0)}{\partial t} = 0. \tag{$C_4$} $$
</P> 

#### 2. Formulation forte

<li> Équation d'équilibre : 

<a id="eq4"></a>
$$\frac{\partial N(x,t)}{\partial x} = \rho S \frac{\partial^2 u(x,t)}{\partial t^2} \tag{$i$} $$

<li> Loi de comportement : 
<a id="eq5"></a>
$$ N(x,t) = ES \varepsilon(x,t) \tag{$ii$} $$

<li> Équation de compatibilité : 

<a id="eq6"></a>
$$ \epsilon(x,t) = \frac{\partial u(x,t)}{\partial x} \tag{$iii$} $$

<p style="text-align: justify;">
En injectant <a href="#eq6">$(i)$</a> dans <a href="#eq5">$(ii)$</a> et puis dans <a href="#eq4">$(iii)$</a> nous obtenons l'équation $\partial$'Alembert qui est une équation différentielle hyperbolique : 
</P>

<a id="eq7"></a>
$$ \frac{\partial^2 u(x,t)}{\partial x^2} - c_d^2 \frac{\partial^2 u(x,t)}{\partial t^2} = 0 \tag{$1$} $$

avec, 

<a id="eq8"></a>
$$ c_d = \sqrt{\frac{E}{\rho}} \, \tag{$2$} $$

la célérité des ondes longitudinales.

#### 3. Formulation faible

<p style="text-align: justify;">
Établissons la formulation faible du problème en introduisant d'abord les ensembles de champs admissibles $\mathcal{U}_{ad}^o$ et admissible à zéro $\mathcal{U}_{ad}^o$ par rapport aux données du problème considéré :
</P> 
<a id="eq9"></a>
$$ \mathcal{U}_{ad} = u(x) \in H^1(\Omega) \quad \vert \quad u(L) = 0 . \tag{$3$}  $$

<p style="text-align: justify;">
L'espace cinématiquement admissible à zéro \(\mathcal{U}_{ad}^o\) est ici équivalent à \(\mathcal{U}_{ad}\) :
</P> 
<a id="eq9"></a>
$$ \mathcal{U}_{ad}^o = \hat u(x) \in H^1(\Omega) \quad \vert \quad \hat u(L) = 0 \tag{$4$} $$

##### Remarque

<p style="text-align: justify;">
Les fonctions $u$ et $\hat{u}$ doivent être suffisamment régulières puisque la formulation locale requiert a priori un champ cinématique admissible deux fois continûment différentiable par morceaux. Ainsi, l'espace de Sobolev $H^1(\Omega)$ dont les fonctions sont de carré intégrable suffit pour vérifier cette condition.
</P>

<p style="text-align: justify;">
\(   \) 
</P>

<p style="text-align: justify;">
Il s'agit maintenant d'exprimer l'équation du mouvement <a href="#eq7">$(1)$</a> par une forme intégrale équivalente par dualisation, c'est-à-dire par multiplication par un champ de déplacement virtuel $\hat{u}$ cinématiquement admissible et par intégration sur le support géométrique $\Omega=[0,L]$, il vient alors :
</P>

<a id="eq9"></a>
$$ \int_{\Omega} \frac{\partial N(x,t)}{\partial x}\hat{u}(x) dx = \int_{\Omega}\rho S \frac{\partial^2 u(x,t)}{\partial t^2}\hat{u}(x) dx \tag{$5$} $$

On utilise ensuite l'Intégration Par Parties, ce qui nous ramène à :

<p style="text-align: justify;">
<a id="eq10"></a>
$$ \big[ N(x,t) \hat{u} \big]_0^L - \int_{\Omega} N(x,t) \frac{d \hat{u}(x)}{d x} dx = \int_{\Omega}\rho S \frac{\partial^2 u(x,t)}{\partial t^2}\hat{u}(x) dx \tag{$6$} $$
</P>

Puis à l'aide des conditions aux limites <a href="#eq1">$(H_1)$</a> et <a href="#eq2">$(H_2)$</a>, on a :
<a id="eq11"></a>
$$ -p(t)S\hat{u}(0) - \int_{\Omega} N(x,t) \frac{d \hat{u}(x)}{d x} dx = \int_{\Omega}\rho S \frac{\partial^2 u(x,t)}{\partial t^2}\hat{u}(x) dx \tag{$7$} $$

Et en utilisant la loi de comportement, il vient :
<a id="eq12"></a>
$$ \int_{\Omega} ES\frac{\partial u(x,t)}{\partial x} \frac{d \hat{u}(x)}{d x} dx + \int_{\Omega}\rho S \frac{\partial^2 u(x,t)}{\partial t^2}\hat{u}(x)  dx = -p(t)S\hat{u}(0) \tag{$8$} $$

Finalement, la formulation faible du problème s'écrit :

Trouver $u$ $\in$ $\mathcal{U}_{ad}$ tel que:

<a id="eq13"></a>
$$ a_k(u,\hat{u}) + a_m\left( \frac{\partial^2 u}{\partial t^2},\hat{u}\right)= l(\hat{u}) \qquad \forall \hat{u} \, \in \, \mathcal{U}_{ad}^o \tag{$9$} $$

avec,

<a id="eq14"></a>
$$ a_k(u,\hat{u}) = ES\int_{\Omega} \frac{\partial u(x,t)}{\partial x} \frac{d \hat{u}(x)}{d x}dx \tag{$10$} $$
<a id="eq15"></a>
$$ a_m\left( \frac{\partial^2 u}{\partial t^2},\hat{u}\right) = \rho S\int_{\Omega} \frac{\partial^2 u(x,t)}{\partial t^2}\hat{u}(x)  dx \tag{$11$} $$
<p style="text-align: justify;">
qui représentent respectivement l'opérateur de raideur modale et l'opérateur de masse modale qui sont tous deux de forme bilinéaire, continue et coercive (définie positive);
</P>
et,
<a id="eq16"></a>
$$ l(\hat{u}) = -p(t)S\hat{u}(0) \tag{$12$}$$

qui est une forme linéaire et continue.

#### 4. Formulation semi-discrète
<p style="text-align: justify;">
La solution analytique de l'équation <a href="#eq13">$(9)$</a> est en général inaccessible puisque les ensembles de champs admissibles sont de dimension infinie. On est donc conduit à chercher une solution approchée par une méthode numérique. On commence par séparer nos variables temporelles et spatiales, soit : 
</P>
<a id="eq17"></a>
$$ u(x,t) = w(x)d(t) \tag{$13$}$$

<p style="text-align: justify;">
La méthode de Rayleigh-Ritz suppose que la solution peut être recherchée comme une combinaison linéaire de fonctions de forme correctement choisies dans l'espace cinématiquement admissible. On pose ainsi :
</P>
<a id="eq18"></a>
\begin{equation}
w(x)=\sum\limits_{i=1}^nN_i(x) \tag{$14$}
\end{equation}

où les $N_i(x)$ sont les modes propres.
<p style="text-align: justify;">
\(   \) 
</P>

En injectant dans <a href="#eq17">$(13)$</a>, on obtient la solution semi-discrétisée en espace :

<a id="eq19"></a>
\begin{equation}
u(x,t)=\sum\limits_{i=1}^nN_i(x)d(t)  \tag{$15$}
\label{u}
\end{equation}

<p style="text-align: justify;">
Ainsi, on \( \textit{discrétise} \) la formulation faible du problème d'ondes \( \textit{continue} \) en introduisant respectivement un sous-espace des fonctions essai $u$ et test $\hat{u}$, $\mathcal{U}_{ad}^h$ $\subset$ $\mathcal{U}_{ad}$, et en posant :
</P>

<a id="eq20"></a>
\begin{equation}
\hat{u}=N_j(x)  \qquad j=[1,..,n]  \tag{$16$}
\end{equation}

On obtient le nouveau problème discrétisé en injectant <a href="#eq19">$(15)$</a>, et <a href="#eq20">$(16)$</a> dans <a href="#eq13">$(9)$</a> : 

<p style="text-align: justify;">
\(   \) 
</P>

Trouver $\hat{u}$ $\in$ $\mathcal{U}_{ad}$, tel que:

\begin{equation}
a_k(N_i,N_j) + a_m(N_i,N_j)= l(N_j) \qquad \forall \hat{u} \, \in \, \mathcal{U}_{ad} \tag{$17$}
\label{faible_N}
\end{equation}

<p style="text-align: justify;">
avec,
</P>

<p style="text-align: justify;">
\begin{equation}
a_k(N_i,N_j) = d(t)\underbrace{ES\int_{\Omega} \frac{d N_i(x)}{d x} \frac{d N_j(x)}{d x}dx}_{K_{ij}} \tag{$18$}
\end{equation}
</P>

<p style="text-align: justify;">
\begin{equation}
a_m(N_i,N_j) =  \ddot{d}(t)\underbrace{\rho S \int_{\Omega} N_i(x) N_j(x)dx}_{M_{ij}} \tag{$19$}
\end{equation}
</P>

et,
\begin{equation}
l(N_j) = \underbrace{-p(t)N_j(0)}_{F_j(t)}  \tag{$20$}
\end{equation}

Que l'on peut écrire sous la forme suivante :
\begin{equation}
M_{ij}\ddot{d}(t) + K_{ij}d(t) = F_j(t) \tag{$21$}
\end{equation}

<a id="sect4"></a>
### IV. Analyse préliminaire

<p style="text-align: justify;">
\(   \) 
</P>

<p style="text-align: justify;">
Il est judicieux de commencer à calculer les différentes caractéristiques fondamentales du système à l'aide d'une analyse élémentaire de propagation d'ondes afin de prédire les ordres de grandeur de la solution cherchée.  
</P>
<p style="text-align: justify;">
\(   \) 
</P>
<p style="text-align: justify;">
<li> Calculons la vitesse $c_d$ de propagation des ondes longitudinales dans la barre à l'aide du résultat suivante:
</P>
<p style="text-align: justify;">

$$c_d = \sqrt{\frac{E}{\rho}} = 5.15.10^3 \, m.s^{-1} \tag{$22$} $$

<p style="text-align: justify;">
<li> Calculons le temps $T_0$ nécessaire à une onde pour parcourir toute la longueur de la barre :
</P>

<p style="text-align: justify;">
\begin{align*} 
&L = 1 \, m \\ 
&T_0 = \frac{L}{c_d} = 1.94.10^{-4} \, s \tag{$23$}
\end{align*}
</P>

<p style="text-align: justify;">
L'explosion génère un impact en $x=0$ provoquant une onde de pression très brève entre les instants $t=0$ et $t=T_p$ se propageant dans la barre à la vitesse $c_d$.
</P>
<li>  Calculons la distance $L_p$ parcourue par cette onde à $t = T_p$ :

<p style="text-align: justify;">
\begin{align*} 
&T_p = 3.88.10^{-5} s\\
&L_p = c_d \, T_p = 0.1999 m \tag{$24$}
\end{align*}
</P>

<p style="text-align: justify;">
Supposons maintenant d'adopter pour le système un modèle de barre 1D et une discrétisation en espace par des éléments finis linéaires.
</P>

<p style="text-align: justify;">
<li> Donnons l'expression des matrices élémentaires $M_e$ et $K_e$,
respectivement de masse et de rigidité.
</P>

<p style="text-align: justify;">
\begin{equation}
([K_e] + [M_e])\Bigg( \begin{matrix}
UX_1  \\
UX_2 \tag{$25$}
\end{matrix}\Bigg) = \begin{pmatrix}
0 \\
0 
\end{pmatrix}
\end{equation}
</P>


Le système à résoudre s'écrit alors, avec une matrice de masse consistante:

<p style="text-align: justify;">
\begin{equation}
\Bigg(\underbrace{\frac{ES}{l_e}\begin{bmatrix}
1 &-1 \\
-1 & 1 
\end{bmatrix}}_{K_e} \, \, \underbrace{- \, \frac{\omega^2 ml_e}{6}\begin{bmatrix} \tag{$26$}
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
</P>

Dans le cas d'une matrice de masse condensée (\textit{lumped mass}), le système s'écrit :

<p style="text-align: justify;">
\begin{equation}
\Bigg(\underbrace{\frac{ES}{l_e}\begin{bmatrix}
1 &-1 \\
-1 & 1 
\end{bmatrix}}_{K_e} \, \, \underbrace{- \, \frac{\omega^2 ml_e}{2}\begin{bmatrix} \tag{$27$}
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
</P>

<p style="text-align: justify;">
<li> Calculons la pulsation la plus élevée d'un élément de longueur $l_e$ lorsque l'on utilise une matrice de masse consistante
</P>

<p style="text-align: justify;">
\begin{equation}
det\Bigg(\frac{6E}{\rho l_e^2}\begin{bmatrix}
1 & -1 \\
-1 & 1 
\end{bmatrix} - \omega^2\begin{bmatrix} \tag{$28$}
2 & 1 \\
1 & 2 
\end{bmatrix}\Bigg) = 0
\end{equation}
</P>

<p style="text-align: justify;">
\begin{equation}
det\Bigg(\begin{matrix}
\frac{6E}{\rho l_e^2} - 2\omega^2 & -\frac{6E}{\rho l_e^2} - \omega^2 \\
-\frac{6E}{\rho l_e^2} - \omega^2 & \frac{6E}{\rho l_e^2} - 2\omega^2 \tag{$29$}
\end{matrix}\Bigg) = \omega^4 -\frac{12E\omega^2}{\rho l_e^2} = 0
\end{equation}
</P>

<p style="text-align: justify;">
\begin{equation}
\Rightarrow \omega_{1,2} = \pm \sqrt{\frac{12E}{\rho l_e^2}} \tag{$30$}
\end{equation}
Nous gardons les pulsations positives
\begin{equation}
\omega = \frac{2\sqrt{3}}{l_e}\sqrt{\frac{E}{\rho}} = 2\sqrt{3} \, \frac{c_v}{l_e}\\ \tag{$31$}
\end{equation}
</P>


<li> Calculons la pulsation la plus élevée d'un élément de longueur $l_e$ lorsque l'on utilise une matrice de masse condensée

<p style="text-align: justify;">
\begin{equation}
det\Bigg(\frac{2E}{\rho l_e^2}\begin{bmatrix}
1 & -1 \\
-1 & 1 
\end{bmatrix} - \omega^2\begin{bmatrix} \tag{$32$}
1 & 0 \\
0 & 1 
\end{bmatrix}\Bigg) = 0
\end{equation}
</P>

<p style="text-align: justify;">
\begin{equation}
det\Bigg(\begin{matrix}
\frac{2E}{\rho l_e^2} - \omega^2 & -\frac{2E}{\rho l_e^2} \\
-\frac{2E}{\rho l_e^2} & \frac{2E}{\rho l_e^2} - \omega^2     \tag{$33$}
\end{matrix}\Bigg) = \omega^4 -\frac{4E\omega^2}{\rho l_e^2} = 0
\end{equation}
</P>

<p style="text-align: justify;">
\begin{equation}
\Rightarrow \omega_{1,2} = \pm \sqrt{\frac{4E}{\rho l_e^2}} \tag{$34$}
\end{equation}
</P>

Nous gardons les pulsations positives

<p style="text-align: justify;">
\begin{equation}
\omega = \frac{2}{l_e}\sqrt{\frac{E}{\rho}} = 2 \, \frac{c_v}{l_e}  \tag{$35$}
\end{equation}
</P>

<p style="text-align: justify;">
Déterminons le pas de temps $\Delta t_0$ pour respecter la condition de stabilité (dite de Courant) lorsque la discrétisation en espace est faite par des éléments de taille uniforme $l_e$, et en temps par la différence finie centrée. Déterminons le $\Delta t_0$ pour les deux cas suivants :
</P>

<p style="text-align: justify;">
\(   \) 
</P>

<p style="text-align: justify;">
\( (i) \) Éléments fins avec matrice de masse consistante :
</P>

<p style="text-align: justify;">
\begin{equation}
\Delta t_0 \leq \frac{\Omega_c}{\omega_h} = \frac{\Omega_c \, l_e}{2 \, c_v} \tag{$36$}
\end{equation}
</P>

Dans ce cas, $\Omega_c = 2$ et $\beta = 0$ , $\gamma = \frac{1}{2}$.
<p style="text-align: justify;">
\(   \) 
</P>
<p style="text-align: justify;">
\( (ii) \) Éléments finis avec matrice de masse condensée : 
</P>

<p style="text-align: justify;">
\begin{equation}
\Delta t_0 \leq \frac{\Omega_c}{\omega_h} = \frac{\Omega_c \, l_e}{2\sqrt{3} \, c_v} \tag{$37$}
\end{equation}
</P>

Dans ce cas, $\Omega_c = 2\sqrt{3}$ et $\beta = \frac{1}{6}$, $\gamma = \frac{1}{2}$.


### V. Algorithme de solution et implémentation dans Cast3M

<p style="text-align: justify;">
1. Donnons sur un schéma les éléments de modélisation nécessaires au calcul :
Comme évoqué dans la partie description du modèle <a href="#sect2">(section II.)</a> , le problème étudié est schématisé comme suit,
</P>
<p style="text-align: justify;">
\(   \) 
</P>

<div style="text-align: center;">
  <figure style="display: inline-block;">
    <img src="/assets/projet1/schema1.png" alt="Modélisation de la barre" width="600"/>
    <figcaption>Figure 1 : Modélisation de la barre</figcaption>
  </figure>
</div>

<p style="text-align: justify;">
\(   \) 
</P>
qui en passant à la discrétisation spatiale nous donne le schéma suivant,
<p style="text-align: justify;">
\(   \) 
</P>

<div style="text-align: center;">
  <figure style="display: inline-block;">
    <img src="/assets/projet1/schema2.png" alt="Modélisation de la barre" width="600"/>
    <figcaption>Figure 2 : Modélisation discrète de la barre</figcaption>
  </figure>
</div>

<p style="text-align: justify;">
\(   \) 
</P>
<p style="text-align: justify;">
ainsi les paramètres d'entrée nécessaires au calcul sont alors implémentés au code, et à l'aide de l'étude préliminaire <a href="#sect4">(section IV.)</a> sur l'analyse des grandeurs caractéristiques du problème, on effectue la discrétisation temporelle.
</P>
<p style="text-align: justify;">
\(   \) 
</P>
<p style="text-align: justify;">
2. Détaillons les conditions aux limites et initiales.
Comme il à été énoncé dans la description du problème <a href="#sect2">(section II.)</a> , nous pouvons écrire les conditions aux limites et initiales.
</P>
<p style="text-align: justify;">
\(   \) 
</P>
<p style="text-align: justify;">
3. Les intégrations temporelles qui peuvent être envisagées sont les suivantes :
</P>
<p style="text-align: justify;">
\(   \) 
</P>

<p style="text-align: justify;">

<div>
\[
\begin{array}{|l||c|c|c|c|c|}
\hline Schéma & \gamma &  \beta &  \Omega_c  & \frac{\Delta T}{T} & Propriétés\\
\hline \text{Newmark 1} & \frac{1}{2} & 0 & 2 & \frac{-w^2h^2}{24} & E \& C-S \\
\text{Newmark 2} & \frac{1}{2} & \frac{1}{4} & \infty & \frac{w^2h^2}{12}  & I-S \\
\text{Newmark 3} & \frac{1}{2} + \alpha & \frac{(1+\alpha)^2}{4} &  \infty  & \left(\frac{1}{12}+\frac{\alpha^2}{4}\right)w^2h^2 & I-S \\
HHT & \gamma_0 & \beta_1  & - & - & I-S\\
\alpha-\text{généralisé} \, {1} & \gamma_1 & \beta_1 & - & - & I-S \\
\alpha-\text{généralisé} \, {2} & \gamma_1 & \beta_1 & - & - & I-S \\ 
\hline 
\end{array}
\]
</div>
</P>
<p style="text-align: justify;">
&#x25E6; Newmark 1 : différence centrée <br>
&#x25E6; Newmark 2 : accélération moyenne <br>
&#x25E6; Newmark 3 : accélération moyenne modifiée <br>
&#x25E6; $\alpha$-généralisé 1 : $\rho = 0.905$ <br>
&#x25E6; $\alpha$-généralisé 2 : $\rho = 0.5$ <br>
&#x25E6; E : Explicite <br>
&#x25E6; I-S : Inconditionnellement Stable <br>
&#x25E6; C-S : Conditionnellement Stable <br>
&#x25E6; $\gamma_0$ = $\frac{1}{2}$ - $\alpha_f$ <br>
&#x25E6; $\gamma_1$ = $\frac{1}{2}$ + $\alpha_m$ - $\alpha_f$ <br>
&#x25E6; $\beta_1$ = $\frac{1}{4}$ + $\frac{1}{2}$($\gamma_m$ - $\gamma_f$) 
</p>

<p style="text-align: justify;">
\(   \) 
</P>
Les trois premières sont des méthodes faisant partie de la famille de Newmark tandis que les trois dernières en sont des méthodes dérivées. En outre, il est à souligner que la méthode HHT est un cas particulier de la méthode $\alpha$-généralisé lorsque $\alpha_m=0$. 
<p style="text-align: justify;">
\(   \) 
</P>
<p style="text-align: justify;">
4. Listons les étapes à suivre dans Cast3M 
<p style="text-align: justify;">
\(   \) 
</P>

<div style="text-align: center;">
  <figure style="display: inline-block;">
    <img src="/assets/projet1/schema3.jpg" alt=" logo castem " width="250"/>
    <figcaption> </figcaption>
  </figure>
</div>

$\texttt{OPTIONS}$<br>
<p style="text-align: justify;">
Tout d'abord on commence par compléter les options du modèle, à savoir la dimension de l'espace, le modèle de calcul, le choix de l'élément à fabriquer et la méthode de calcul des déformations.Puis, le choix de l'algorithme à utiliser, présentés à la question précédente, est laissé libre à l'utilisateur.<br>
<p style="text-align: justify;">
\(   \) 
</P>

$\texttt{DONNÉES}$<br>
<p style="text-align: justify;">
On entre ensuite les données du problème, autrement dit la géométrie, les propriétés matérielles, le chargement - en particulier la durée de l'explosion $T_p$ - et le nombre d'éléments à discrétiser.<br>
</P>
<p style="text-align: justify;">
\(   \) 
</P>

$\texttt{MAILLAGE}$<br>
<p style="text-align: justify;">
Ceci étant fait, le maillage est ainsi réalisé.<br>
<p style="text-align: justify;">
\(   \) 
</P>

$\texttt{MODÈLE}$<br>
<p style="text-align: justify;">
Le modèle du problème est alors construit à partir des propriétés géométriques et matérielles définies précédemment. Les matrices de raideur $K$ et de masse $M$ sont calculées auxquelles les conditions aux limites sont appliquées.<br>
</P>
<p style="text-align: justify;">
\(   \) 
</P>

$\texttt{CHARGEMENT}$<br>
<p style="text-align: justify;">
De ce fait, on représente alors le step de pression avant et après l'explosion, ainsi que le chargement impliqué de ladite explosion.<br>
</P>
<p style="text-align: justify;">
\(   \) 
</P>

$\texttt{DICRÉTISATION TEMPORELLE}$<br>
<p style="text-align: justify;">
Nous procédons a posteriori à la discrétisation temporelle en veillant à sélectionner la vitesse de propagation $c_d$, une matrice de masse consistante ou lumpée ainsi que le pas de temps critique et le temps $T_0$ mis par l'onde pour traverser la barre. Puis, on choisit un pas de temps et on écrit les conditions initiales.<br>
</P>
<p style="text-align: justify;">
\(   \) 
</P>

$\texttt{RÉSOLUTION}$<br>
<p style="text-align: justify;">
À présent, nous pouvons effectuer la résolution de ce problème.La procédure \texttt{DYNAMIC} permet de réaliser un calcul dynamique pas à pas pour les algorithmes vus précédemment. Ceci est fait via les données du problème, à savoir les déplacement et vitesse initiaux, les matrices de rigidité, de masse, et d'amortissement. On va, de plus, sauvegarder les instants de calcul. Suivant le choix du schéma d'intégration, on va entrer les paramètres $\alpha$, $\gamma$, $\alpha$ ou $\rho^\infty$.<br>
</P>
<p style="text-align: justify;">
\(   \) 
</P>

$\texttt{POST-TRAITEMENT}$<br>
<p style="text-align: justify;">
Finalement, on post-traite les résultats intéressants, que l'on décrira à la question suivante, et on les représente graphiquement.
</P>

<p style="text-align: justify;">
5. Les grandeurs pertinentes à post-traiter sont:
</P>
<p style="text-align: justify;">
\(   \) 
</P>
<p style="text-align: justify;">
&#9744; le déplacement longitudinal<br>
&#9744; le déplacement transversal<br> 
&#9744; la vitesse longitudinale<br> 
&#9744; l'effort normal<br>
&#9744; le travail des efforts extérieurs<br>
&#9744; le travail des efforts intérieurs<br>
&#9744; le bilan d'énergie<br>
</P>
<p style="text-align: justify;">
\(   \) 
</P>
<p style="text-align: justify;">
6. Écrivons le programme gibiane pour Cast3M (cf. fichier .dgibi)
</P>

### VI. Résolution numérique, vérification et analyse

<p style="text-align: justify;">
1. Vérifions les résultats obtenus et en particulier l'effet du pas de temps.
</P>
<p style="text-align: justify;">
\(   \) 
</P>
<p style="text-align: justify;">
Le programme a été lancé pour les différents algorithmes et en faisant varier le pas de temps tout en veillant à ce qu'il soit maintenu inférieur au pas de temps critique déterminé via la condition de stabilité.
Leurs représentations sont jointes dans les fichiers \( \textit{postscript} \).
</P>
<p style="text-align: justify;">
\(   \) 
</P>
<p style="text-align: justify;">
Notre première étude portera sur un pas de temps égal au pas de temps critique.
Pour le premier algorithme, le schéma de Newmark différence centrée, les grandeurs tracées présentent de fortes variations ou alors un seul pic tout le long de la barre. Ceci s'explique par le fait que ce schéma est conditionnellement stable, il faut donc un pas de temps inférieur au pas de temps critique pour ce schéma afin de visualiser une solution correcte.
Les autres schémas présentent des résultats plus cohérents, et en particulier le schéma de Newmark accélération moyenne modifiée, étant donné qu'ils sont tous inconditionnellement stables.
</P>
<p style="text-align: justify;">
\(   \) 
</P>
<p style="text-align: justify;">
Des pas de temps cinq et dix fois inférieurs ont aussi été traités mais l'effet du pas temps sera particulièrement observable pour le premier schéma dans la mesure où il s'agit du seul schéma conditionnellement stable. En effet, nous constatons que ce schéma présente une solution convenable uniquement dès lors que le pas de temps est près de deux fois inférieur au pas de temps critique.
</P>
<p style="text-align: justify;">
\(   \) 
</P>
<p style="text-align: justify;">
2.  Analysons la solution du problème en la comparant avec l'analyse de l'exemple du chapitre 9.4 du livre \( \textit{Getting Started with Abaqus: Interactive Edition} \). La solution donne la contrainte longitudinale selon $\textbf{e}_3$ soit (S33) en trois points différents le long de la barre ($0.25 m$, $0.5 m$, et $0.75 m$)
</P>
<p style="text-align: justify;">
\(   \) 
</P>
<p style="text-align: justify;">

<div style="text-align: center;">
  <figure style="display: inline-block;">
    <img src="/assets/projet1/schema4.png" alt="Figure 3 : Solution en contrainte sur Abaqus" width="500"/>
    <figcaption>Figure 3 : Solution en contrainte sur Abaqus</figcaption>
  </figure>
</div>

<div style="text-align: center;">
  <figure style="display: inline-block;">
    <img src="/assets/projet1/schema5.png" alt="Figure 3 : Solution en contrainte sur castem" width="500"/>
    <figcaption>Figure 4 : Solution en contrainte sur Cast3m</figcaption>
  </figure>
</div>


</P>
<p style="text-align: justify;">
\(   \) 
</P>
<p style="text-align: justify;">
\underline{Conclusion} : Nous obtenons une solution quasi-identique pour le même problème résolu par un logiciel de programmation (Cast3m) et un logiciel de simulation (Abaqus), ce qui permet de valider nos résultats.
</P>
<p style="text-align: justify;">
\(   \) 
</P>
<p style="text-align: justify;">
3. Comparons les avantages et les inconvénients des méthodes explicite et implicite.
</P>
<p style="text-align: justify;">
\(   \) 
</P>
<p style="text-align: justify;">
La résolution d’un problème de dynamique s'interroge tout d'abord sur la nature du problème posé afin de choisir le schéma d’intégration temporelle, qui définit la discrétisation en temps de l’équation d’équilibre globale.
</P>
<p style="text-align: justify;">
\(   \) 
</P>
<p style="text-align: justify;">
</P>
Les schémas sont classés en deux grandes familles: les schémas explicites et les schémas implicites. 
<p style="text-align: justify;">
\(   \) 
</P>
<p style="text-align: justify;">
Les premiers utilisent un pas de temps très réduit et ont pour avantage de permettre la résolution de l’équilibre en dynamique rapide, pour des phénomènes très violents (chocs, explosions, ...). L'avantage dudit schéma réside dans l'adaptation du phénomène à observer qui impose déjà un pas de temps très petit (souvent de l’ordre de la microseconde). L’inconvénient du schéma explicite est qu’il existe un pas de temps critique à ne pas dépasser et les temps de calcul sont très longs.
</P>
<p style="text-align: justify;">
\(   \) 
</P>
<p style="text-align: justify;">
D'autre part, les phénomènes rapides induisent de fortes non-linéarités sur un temps très court, le schéma explicite est donc beaucoup plus efficace que le schéma implicite.
</P>
<p style="text-align: justify;">
\(   \) 
</P>
<p style="text-align: justify;">
En revanche, les schémas implicites sont utiles en dynamique lente (vibrations, ébranlements des structures), on est en présence d’un comportement linéaire des matériaux, pour une durée des phénomènes s’étalant sur plusieurs secondes. Les schémas implicites ont donc cette fois-ci tout leur intérêt car ils permettent d’utiliser des pas de temps très grands et sont alors très rapides en temps de calcul. La solution numérique reste stable et sa qualité est contrôlée par le critère de résidu en équilibre, tout comme en quasi-statique. Cependant, la convergence s’avère parfois très délicate en présence de fortes non-linéarités. 
</P>
<p style="text-align: justify;">
\(   \) 
</P>
<p style="text-align: justify;">
4. Réanalysons le problème dans le cas où la barre est libre à l'extrémité droite et non plus encastrée.
</P>
<p style="text-align: justify;">
\(   \) 
</P>
Nous pouvons, par exemple, comparer le déplacement selon $x$ pour les deux cas.

<div style="text-align: center;">
  <figure style="display: inline-block;">
    <img src="/assets/projet1/libre.png" alt="Déplacement  bord libre" width="500"/>
    <figcaption>Figure 5 : Déplacement $u_x$ : bord libre</figcaption>
  </figure>
</div>

<div style="text-align: center;">
  <figure style="display: inline-block;">
    <img src="/assets/projet1/encastre.png" alt="Déplacement bord encastré" width="500"/>
    <figcaption>Figure 6 : Déplacement $u_x$ : bord encastré</figcaption>
  </figure>
</div>

<p style="text-align: justify;">
\(   \) 
</P>
<p style="text-align: justify;">
Nous remarquons que l'onde en $x = 0$ va s'annuler puis devenir négative dans le cas d'un encastrement; et s'additionner dans le cas d'un bord libre, car l'onde est dans ce cas non inversé comme schématisé ci-dessous:
</P>

<div style="text-align: center;">
  <figure style="display: inline-block;">
    <img src="/assets/projet1/onde.png" alt="Comportement de l'onde" width="500"/>
    <figcaption>Figure 7 : Comportement de l'onde</figcaption>
  </figure>
</div>

### VII. Conclusion

<p style="text-align: justify;">
Nous avons, au cours de ces travaux, étudié la dynamique d'une structure de type barre soumise à une explosion à une extrémité et encastrée à l'autre extrémité. Nous avons tout d'abord mis en place la formulation faible du problème et après avoir fait une analyse préliminaire afin de connaître les grandeurs caractéristiques du problème posé, nous avons établi l'expression du problème d'équilibre à résoudre à l'aide de la procédure $\texttt{DYNAMIC}$ de Cast3m. Ce travail nous à permis d'utiliser les méthodes d'intégration temporelle issues des schémas explicites. Ces schémas sont bien adaptés pour l'application à des problèmes de dynamique rapide. Ainsi après avoir complété le code gibiane, nous avons pu mettre en évidence l'impact du pas de temps suivant différents algorithmes d'intégration. Les résultats ont aussi été validé par une comparaison du même problème traité sur Abaqus. Nous avons aussi calculé la réponse en changeant les conditions aux limites, pour une extrémité libre, les ondes augmentent en amplitude, le risque de fracture est donc plus important.
</P>