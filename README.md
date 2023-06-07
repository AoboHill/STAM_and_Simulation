# STAM_and_Simulation

### Author: Aobo Wang from Australian National University

## This repository provides the simulation R code for the STAM algorithm and the source JAVA code for STAM.

```
All of the STAM code is included in the STAM file folder.
```

## STAM
```
1. For running STAM to do the species tree inference, users are required to have BEAST2 version='2.6.7'.
2. Users may also be required to install the SNAPP (Bryant et al., 2012) and SNAPPER (Stolz et al., 2021).
```

## Simulation
```
For the simulation:

1. SimulateSeq.R is to sample the sequences at tips using discrete WF.
2. Simulation_MOM.R is to assess the error of moment approximation approach. 
   The true distribution is sampled using discrete WF.
3. Simulation_test_chebyshev is to asess the performance by grouping allele frequencies into bins
4. Simulation_test_chebyshev_alternative is to asess the performance by grouping allele frequencies into bins using alternative formulas.

Note that to run the simulation code may require load the functions by order.
```

## Scaling

```
About the scaling:

r (number of the generations)
N (population size)
U (transition probability matrix)
Q (transition rate matrix)
I (identity matrix)

sclaing:

1. t = r/N

2. Q = N(U - I)

This scaling both used in STAM and corresponding simulations. Moreover, in the STAM, 
the mu in Q (Kimura two-parameter model) is fixed whereas the kappa is the only parameter required to be estimated

```
