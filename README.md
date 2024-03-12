# Exact Diagonalization

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)

*Solving fermionic lattice models using exact diagonalization*

The project was originally developed in order to perform exact diagonalization caculations for fermionic models, specifically for Cluster Mott Insulators but can be extended for other models as well. 

The fock space for states is represented in binary
with following layers.

$$
\begin{aligned}

\text{Spin} &\rightarrow [\uparrow][\downarrow] = [0][1] \\

\text{Orbital} &\rightarrow n_{orbital} \; \text{Spin} \\

\text{Site} &\rightarrow n_{site} \; \text{Orbital}
\end{aligned}
$$

A spin 1/2 system occupies $2n_{orbitals} \cdot n_{sites}$ bits.

> Note: The module only supports upto 64 bits at current time.