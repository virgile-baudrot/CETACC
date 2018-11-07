# CETACC
Cetacean Accumulation


## Ingestion part

The dynamic of the internal concentration of the contaminant, also known as the toxicokinetics, may be describe by this simple equation:

\[
\frac{dC_{j,in}}{dt} = TrIng_{j,tot} - k_{j, out}  C_{j,in}
\]

\[
\frac{dC_{a,in}}{dt} = TrIng_{a, tot} - k_{a, out} C_{a,in}
\]


where:
- $C_{j, in}$ and $C_{a,in}$ are respectivelly the internal concentration for *juveniles* and *adults*
- $TrIng_{j, tot}$ and $TrIng_{a, tot}$  are the trophical ingestion of contaminant for *juveniles* and *adults* respectivelly (see later for details)
- $k_{j,out}$ and $k_{a, out}$ the excretion rate of the contaminant for *juveniles* and *adults* respectivelly

Note that since $TrIng_{a, tot}$ is constant, we have:

\[
C_{i,in} = \frac{TrIng_{i, tot}}{k_{out}}\left\( 1 - e^{- k_{i,out} t} \right\)
\]


### About $TrIng_{\dot, tot}$

A juvenile is only exposed to the contaminant through the maternal feeding (milk). Note that a new-born has likely been exposed through maternal gestation.

\[
TrIng_{j, tot} = \eta_{j} \times I_{maternal} \times C_{maternal}
\]

And, for $n$ prey species, and adult is exposed to the contaminant through food items:

\[
TrIng_{a, tot} = \eta_{a} \sum_{i=1}^{n} \times I_{i} \times C_{i}
\]

In both equations, we have:
- $\eta_{j}$ and $\eta_{a}$ the assimilation rate in *juveniles* and *adults*,
- $I_i$: the ingestion rate of item $i$ (e.g., $kg.day^{-1}$),
- $C_i$: the concentration in item $i$ (e.g., $mg.kg^{-1}$).






