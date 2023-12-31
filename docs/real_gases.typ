#set par(leading: 0.55em, first-line-indent: 1.8em, justify: true)
#set text(font: "New Computer Modern")
#show raw: set text(font: "New Computer Modern Mono")

= Ideal Gas Equation

#link("https://demonstrations.wolfram.com/NewtonsLawOfCooling/")[
Wolfram's Article on Newton's law of Cooling]

$ P V  = n R T (kelvin) $

$ P tilde(V) = R T (kelvin) $

= Real Gases 

Real Gases do not behave according to ideal gas law. We often 
use more complicated equations of 
state @felder2020elementary @kontogeorgis2019taking 
@coleman2023distillation @adewumi2023phaserelations.

Key assumptions in ideal gas are @adewumi2023phaserelations:

- gas molecules or atoms take no space
- no attractive/repulsive forces between gas molecules or atoms


== Van Der Waals
To make the ideal gas more realistic, Van Der Waals equation 
accounts for space between gas molecules/atoms and the forces 
between them as well @adewumi2023phaserelations. The Equation 
of State (EOS) looks like @felder2020elementary:

$ (P + a/ (tilde(V)^2))(tilde(V) - b)  = R T (kelvin) $

$ a = (27 R^2 T_c^2) / (64 P_c ) $
$ b = (R T_c)/ (8 P_c) $


== SRK EOS

Gas/vapour molecules are not uniformly shaped. The force between
gas molecules depends on the shape. Redlich Kwong and Soave Redlich 
Kwong (SRK) equation of state helps account for this 
shape problem @adewumi2023phaserelations. This is accounted for 
in the accentricity factor ($omega$) in
the SRK equation of state, a semi empirical EOS. 
The SRK EOS is written as @felder2020elementary:

$ (P + (a alpha)/ (tilde(V)(tilde(V) + b)))(tilde(V) - b)  = R T (kelvin) $

$ alpha = [1 + m (1 -sqrt(T_r))]^2 $

$ m = 0.48508 + 1.55171 omega - 0.1561 omega^2 $

$ T_r = (T(K))/T_c $

$ a = 0.42747 ( R^2 T_c^2) / ( P_c ) $
$ b = 0.08664 (R T_c)/ ( P_c) $


== Peng Robinson EOS

The Peng Robinson EOS (PR-EOS) is an equation of state (EOS)
meant for oil and gas industry use @lopez2017peng 
@adewumi2023phaserelations. It is used 
for molecules such as hydrogen sulfide, oxygen, nitrogen and
several hydrocarbons @sandler2017chemical.

The equation is as follows:

$ (P + (a(T))/ 
(tilde(V)(tilde(V) + b) +
b(tilde(V)-b)))(tilde(V) - b)  = R T (kelvin) $

$ a(T) = 0.45724 ( R^2 T_c^2) / ( P_c ) alpha(T) $
$ alpha(T) = [1 + kappa (1 -sqrt(T_r))]^2 $

$ kappa = 0.37464 + 1.54226 omega - 0.26992 omega^2 $

$ T_r = (T(K))/T_c $

$ b = 0.07780 (R T_c)/ ( P_c) $


== Solving for P, $tilde(V)$ and T

=== For Pressure and Temeprature, quite easy... 

As long as $T_c$, $P_c$, $omega$ and $R$ are known...

For pressure you need $tilde(V)$ and $T$;

For temperature you need $tilde(V)$ and $P$
solve for a quadratic equation in $sqrt(T(K))$


=== What about Volume?

You'll need to solve cubic equations of state.

Using the compressibilty factor $Z$, where

$ Z = (P tilde(V))/(R T) $

For the Peng Robinson EOS,
you can solve for the cubic polynomial @sandler2017chemical:

$ Z^3 - (1-B) Z^2 + (A- 3B^2 -  2B) Z - (A B -B^2 - B^3)=0 $

Where 

$ B = (b P) / (R T) $
$ A = (a P)/(R^2 T^2) $

Solution procedure, solve for Z, then solve for $tilde(V)$

How to solve for Z?

If you want to derive things yourself, try:

$ (x-a_1)(x-a_2)(x-a_3) = 0 $

Multiply throughout, and compare coefficients 
with a typical cubic formula:
$ x^3 + b x^2 + c x + d = 0 $

Note that $a_1$, $a_2$ and $a_3$ are in general, complex numbers.
At least one of them is real if b, c and d are real.


Cubic Equation formula @perry2007perry:

For the equation:

$ x^3 + b x^2 + c x + d = 0 $

We first calculate a term $r$:

$ r = (p/3)^3 + (q/2)^2 $

Where:

$ p = (3c - b^2)/3 $
$ q = (27d - 9 b c + 2b^3)/27 $

If r < 0, we have three real roots:

$ x_k = y_k - b/3 $

Where,

$ y_k = minus.plus 2 sqrt(-p/3) cos [phi/3 + 120k ] $

Where k = 0, 1, 2

and $ phi = cos^(-1) sqrt((q^2 /4)/(-p^3/27)) $

if q  > 0

$ y_k = - 2 sqrt(-p/3) cos [phi/3 + 120k ] $

if q < 0 

$ y_k = + 2 sqrt(-p/3) cos [phi/3 + 120k ] $

Is the 120k a degree or radian??
We'll have to do some testing with GNU octave..


If $r>0$, expect complex roots and one real root,

$ y_1 = A + B $

Where 

$ A = root(3,-q/2 + sqrt(r)) $
$ B = root(3,-q/2 - sqrt(r)) $

$ y_2 = -1/2(A+B) + i sqrt(3)/2 (A-B) $
$ y_3 = -1/2(A+B) - i sqrt(3)/2 (A-B) $


If $r=0$, expect three real roots, at least two are 
equal

=== Test cases for cubic solving

We will need some verification tests to ensure our 
code works correctly
$ x^3 - 7 x + 7 = 0 $

Roots = $mat(-3.0489;
 1.6920;
 1.3569)$



$ x^3 + 8 x^2 - 7 x + 9 = 0 $

Roots = $mat(  -8.9001 ;
   0.4501 + 0.8993i;
   0.4501 - 0.8993i)$

#bibliography("biblio.bib",
style: "harvard-cite-them-right")

