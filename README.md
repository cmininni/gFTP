# gFTP

MATLAB implementation of the generalised Firing-to-Parameter algorithm (gFTP) introduced in "Constructing neural networks with pre-specified dynamics".  

Function gFTP.m takes as input a structure array named "data" with fields: G (the transition graph in matrix format), mode ('consistency' or 'construction') and N_neu_min.
It returns G_cons, matrices Z_s, Z_t, W_y, W_r, and structure array data_out, with the time spend constructing G_cons and constructing the activation and synaptic weight matrices.

# An example
First, we construct a graph $G$ by defining its associated matrix $\mathbf{G}$
```Matlab
G = [1 1 1
     2 1 2
     3 1 1
     1 2 2
     2 2 3
     3 2 1
     1 3 2
     2 3 3
     3 3 3
     1 4 4
     2 4 1
     3 4 2];
```
We plot the graph with the tgPlot function:
```Matlab
tgPlot(G);
title('$G$','interpreter','latex')
```
