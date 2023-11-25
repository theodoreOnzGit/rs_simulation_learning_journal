= Newton's law of Cooling


$ Q = h A (T_"cool"-T_s)$

#link("https://demonstrations.wolfram.com/NewtonsLawOfCooling/")[
Wolfram's Article on Newton's law of Cooling]

Differential Equation: 

$ m c_p (d T_"solid")/(d t) 
= - h A (T_"solid" - T_"surroundings") $

$ m c_p (T_"solid"^(t+ Delta t) - T_"solid"^(t))/(Delta t) 
= - h A (T_"solid" - T_"surroundings") $

== Implicit Euler Method

$ m c_p (T_"solid"^(t+ Delta t) - T_"solid"^(t))/(Delta t) 
= - h A (T_"solid"^(t+ Delta t) - T_"surroundings") $

== Explicit Euler Method

$ m c_p (T_"solid"^(t+ Delta t) - T_"solid"^(t))/(Delta t) 
= - h A (T_"solid"^(t) - T_"surroundings") $

== Let's do Implicit Method

$ m c_p (T_"solid"^(t+ Delta t) 
- T_"solid"^(t))/(Delta t) +h A T_"solid"^(t+ Delta t)
=    h A T_"surroundings" $

$ m c_p (T_"solid"^(t+ Delta t) )/(Delta t) +h A T_"solid"^(t+ Delta t)
=    h A T_"surroundings" + m c_p (T_"solid"^(t))/(Delta t) $


$ T_"solid"^(t+ Delta t) [ (m c_p )/(Delta t) +h A ]
=    h A T_"surroundings" + m c_p (T_"solid"^(t))/(Delta t) $

$ T_"solid"^(t+ Delta t) 
=   1/[ (m c_p )/(Delta t) +h A ] 
(h A T_"surroundings" + m c_p (T_"solid"^(t))/(Delta t)) $

=== Time stepping...

At t=0, $Delta t = 0.1 s$

$ T_"solid"^(0.1 s) 
=   1/[ (m c_p )/(Delta t) +h A ] 
(h A T_"surroundings" + m c_p (T_"solid"^(0.0 s))/(Delta t)) $

At t=0.1

$ T_"solid"^(0.2 s) 
=   1/[ (m c_p )/(Delta t) +h A ] 
(h A T_"surroundings" + m c_p (T_"solid"^(0.1 s))/(Delta t)) $













