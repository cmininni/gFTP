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
<img src="https://github.com/cmininni/gFTP/assets/167989205/b0a756b3-847c-4dfb-beb8-98ea98d98d2a" width=50% height=50%>


Next, we construct the input structure data_in. We will pass this structure as input argument to gFTP
```Matlab
data_in.G = G;
data_in.mode = 'construction';     %to construct the network; use 'consistency' to only generate $G_{cons}$
data_in.N_neu_min = [];            %no minimum number of neurons imposed
```
Call gFTP
```Matlab
[G_cons, Y, Z_s, Z_t, W_y, W_r, data_out] = gFTP(data_in);
```
Check if the output network actually follows graph G_cons
```Matlab
is_isomorf=check_dynamics(G_cons,Y',Z_s',Z_t')

is_isomorf =

  logical

   1
```
check consistency and construction times in data_out
```Matlab
data_out

data_out = 

  struct with fields:

       consistency_time: 0.2358
      construction_time: 0.5145
    perceptron_training: 1
```

Plot consistent graph G_cons
```Matlab
tgPlot(G_cons)
title('$G_{cons}$','interpreter','latex')
```
![fig2](https://github.com/cmininni/gFTP/assets/167989205/74bf02d0-b29a-4b17-82da-06982733dcb7)
<img src="https://github.com/cmininni/gFTP/assets/167989205/b0a756b3-847c-4dfb-beb8-98ea98d98d2a" width=50% height=50%>

Construct and plot auxiliary graph D with function make_D and plot_D
```Matlab
D = make_D(G);
D_cons = make_D(G_cons);
plot_D(D);
title('D(G)','interpreter','latex')
plot_D(D_cons);
title('$D(G_{cons})$','interpreter','latex')
```
![fig3](https://github.com/cmininni/gFTP/assets/167989205/2823b6c7-ee8d-44d3-800d-ffb9c3e15bb1)
![fig4](https://github.com/cmininni/gFTP/assets/167989205/60b0f25c-8f54-4fb2-9a93-c043437928a0)

