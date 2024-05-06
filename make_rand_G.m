%Generates random graph of N_v nodes and N_s inputs
function g = make_rand_G(N_s, N_v)
y=[];
vi=[];
vd=[];
vd_obligados=randperm(N_v);
for i=1:N_v
    y=[y;(1:N_s)'];
    vi=[vi;repmat(i,N_s,1)];
    aux=mod(i-1+linspace(-2,2,5),N_v)+1;
    vd=[vd;[randsample(aux,N_s-1,true)';vd_obligados(i)]];
end
g = [y, vi, vd];
