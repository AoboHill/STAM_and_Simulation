# STAM_and_Simulation

## This repository provides the simulation R code for the STAM algorithm.

```
All of the STAM code is included in the STAM file folder.
```

## Simulation
```
For the simulation:

1. simulateSeq.R is to sample the sequences at tips using discrete WF.

2. Simulation_MOM.R is to assess the error of moment approximation approach. 
   The true distribution is sampled using discrete WF.
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
