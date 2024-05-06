% gFTP (Generalized Firing to Parameter).
% MATLAB implementation of the algorithm introduced in "Constructing neural
% networks with pre-specified dynamics".
% 
% INPUT: structure array "data_in" with fields:
%
%        G:     Matrix of N_transition rows and 3 columns, specifying
%               stimulus, source node and target node of each transition.

%        mode:  String variable.
%               If data.mode == 'construction', gFTP
%               construct consistent graph G_cons together with matrices Z_s, Z_t,
%               W_y, W_r.
%               If data.mode=='consistency', gFTP only generates G_cons

%        N_neu_min:  Minimum number of neurons the output network will have.
%                  If N_neu_min is empty, no minimum number of neurons is imposed.

% OUTPUT:   G_cons: consistent graph computed from input graph G.

%           Y:      Matrix with N_stimuli rows and N_transition columns.
%           The kth column is the one-hot encoded vector associated with the
%           stimulus presented in the kth transition.

%           Z_s: Matrix with N_neu rows and N_transition columns.
%           The kth column is the network activation vector of the source node of the
%           kth transition.

%           Z_t: Matrix with N_neu rows and N_transition columns.
%           The kth column is the network activation vector of the target node of the
%           kth transition.

%           W_y: Matrix with N_neu rows and N_stimuli columns. The ith
%           row holds the synaptic weights from stimuli to the ith neuron in
%           the recurrent network.

%           W_r: Matrix with N_neu rows and N_neu columns. The ith row
%           holds the synaptic weights from all neurons to the ith neuron in
%           the recurrent network.

%           data_out: structure array with fields:
                    %consistency_time: time in seconds to build a
                    %consistent graph construction_time: time(sec) to build
                    %a consistent graph and construct matrices Z_s, Z_t,
                    %W_y, W_r.
                    
function [G_cons, Y, Z_s, Z_t, W_y, W_r, data_out] = gFTP(data_in)
g = remueve_v_excedentes(data_in.G);
[G_cons, consistencia_y, Y, Z_s, Z_t, W, U, isomorf_1, isomorf_2, nCambios, data_out] = f_global(g, data_in.mode, data_in);
Y=Y';
Z_s=Z_s';
Z_t=Z_t';
W_y=W(1:size(Y,1),:)';
W_r=W(size(Y,1)+1:end,:)';
end

function [g, consistencia_y, y, Zi, Zd, w, U, isomorf_1, isomorf_2, nCambios, data_out_2]=f_global(g,modo, data_in)

clear global data_out
global data_out

salida=false;
c_solap=[];
d=[];
iter=0;
nS=max(g(:,1));
while ~salida
    iter=iter+1;
    tic
    [G,T]=g2D(g);
    C=detecta_solap_1(G, nS);
    C=eye(size(C))==1;
    if ~isempty(G)
        while 1
            [c_solap,d]=solap_indirecto(G,T,C);
            [G_inv,T,C2]=aplico_solap_ind(G,T,C,c_solap,d);
            if all(C2(:)==C(:))
                break
            else
                C=C2;
            end
        end
        C=C2;
        G=G_inv;
        [g, hay_ciclos, nPasos]=remueve_ciclos_2(g,G,T,C);
    end
    data_out.consistency_time=toc;
    salida=true;
end

tic
switch modo
    case 'construction'
        [g, consistencia_y, y, Zi, Zd, w, U, isomorf_1, isomorf_2, nCambios]=make_zd(g, data_in);
        data_out.construction_time=toc;
        data_out.perceptron_training=nCambios;
    case 'consistency'
        data_out.construction_time=toc;
        data_out.perceptron_training=true;
        consistencia_y=[];
        y=[];
        Zi=[];
        Zd=[];
        w=[];
        U=[];
        isomorf_1=1;
        isomorf_2=1;
        nCambios=1;
end
data_out_2=data_out;
end

function C=detecta_solap_1(G,nS)
colores=unique(G(:,3));
solap=[];
nColores=(nS^2-nS)/2;
C=eye(nColores)==1;

for c1=1:length(colores)
    ind1=G(:,3)==colores(c1);
    for c2=c1+1:length(colores)
        ind2=G(:,3)==colores(c2);
        aux=ismember(G(ind1,[1,2]),G(ind2,[1,2]),'rows')';
        if any(aux)
            solap=[solap;[colores(c1),colores(c2)]];
            C(colores(c1),colores(c2))=true;
            C(colores(c2),colores(c1))=true;
        end
    end
end

while 1
    C_aux=C;
    for c1=1:size(C,1)
        C(c1,:)=C(c1,:) | any(C(C(c1,:),:),1);
    end
    if ~any(C_aux(:)~= C(:))
        break
    end
end
end

function [c_solap,d]=solap_indirecto(G,T,C)
G_ori=G;
T_ori=T;
rows_G=1:size(G,1);
ind_ini=randi(size(G,1));
v=G(ind_ini,1);
c=G(ind_ini,3);
vi_f=T(ind_ini,1);
n_c=unique(G(:,3));

v_h_f=[];
c_h_f=[];
vi_h=[];

v_h_b=[];
c_h_b=[];
vi_b=[];

k1=1;
ind_last_trans=[];

c_solap=[];
d=[];
terminal=false;
salida=false;
if length(C)>1
    while ~salida
        
        k1=k1+1;
        
        c_tran=find(C(c,:));
        ind1=find(G(:,1)==v & any(G(:,3)==c_tran,2));
        
        if ~terminal
            if ~isempty(ind1) %si puedo avanzar, avanzo
                v_h_f=[v_h_f,v];
                c_h_f=[c_h_f,c];
                vi_h=[vi_h,vi_f];
                
                ind2=randi(length(ind1));
                v=G(ind1(ind2),2);
                c=G(ind1(ind2),3);
                vi_f=T(ind1(ind2),1);
                
                
                ind_last_trans=[ind_last_trans,rows_G(ind1(ind2))];
                
                G(ind1(ind2),:)=[];
                T(ind1(ind2),:)=[];
                rows_G(ind1(ind2))=[];
                
                if sum(v_h_f==v)>0
                    terminal=true;
                end
            else 
                terminal=true;
                v_h_f(end+1)=v;
                c_h_f(end+1)=c;
                v=v_h_f(1);
            end
        else
            
            c_tran=find(C(c,:));
            ind1=find(G(:,2)==v & any(G(:,3)==c_tran,2));
            if ~isempty(ind1)
                
                ind2=randi(length(ind1));
                v=G(ind1(ind2),1);
                c=G(ind1(ind2),3);
                vi=T(ind1(ind2),1);
                
                G(ind1(ind2),:)=[];
                T(ind1(ind2),:)=[];
                rows_G(ind1(ind2))=[];
                
                v_h_b=[v,v_h_b];
                c_h_b=[c,c_h_b];
                vi_b=[vi,vi_b];
                
                
            else
                p=[v_h_b,v_h_f];
                c_valid=find(~C(c,:));
                m_c_valid=any(G(:,3)==c_valid,2);
                for c1=1:length(p)
                    for c2=c1+1:length(p)
                        v1=p(c1);
                        v2=p(c2);
                        [c_solap_aux,d_aux,p_arc]=busca_camino(G,v1,v2,m_c_valid,rows_G);
                        if ~isempty(c_solap_aux)
                            c_solap=[c_solap;[c,c_solap_aux]];
                            c_usados=find(any(C(c_solap(:),:),1));
                            d=[d;d_aux];
                            if isempty(setdiff(c_valid,c_usados))
                                salida=true;
                                break
                            end
                        end
                    end
                    if salida
                        break
                    end
                end
                ind_ini=randi(size(G,1));
                
                v=G(ind_ini,1);
                c=G(ind_ini,3);
                terminal=false;
                
                ind_ini=randi(size(G,1));
                v=G(ind_ini,1);
                c=G(ind_ini,3);
                vi_f=T(ind_ini,1);
                
                v_h_f=[];
                c_h_f=[];
                vi_h=[];
                
                v_h_b=[];
                c_h_b=[];
                vi_b=[];
                
                ind_last_trans=[];
            end
        end
        if isempty(G)
            break
        end
    end
end
end

function [G,T,C]=aplico_solap_ind(G,T,C,c_solap,d)
for c1=1:size(c_solap,1)
    C(C(c_solap(c1,1),:),C(c_solap(c1,2),:))=true;
    C(C(c_solap(c1,2),:),C(c_solap(c1,1),:))=true;
    if d(c1)<0
        ind_c=G(:,3)==c_solap(c1,2);
        G(ind_c,[1,2])=G(ind_c,[2,1]);
        T(ind_c,[2,3])=T(ind_c,[3,2]);
    end
end
end

function [c_solap,d,p]=busca_camino(G,v1,v2,m_c_valid,rows_G)
v=v1;
v_h=[];
c_h=[];
d_h=[];
ind_list=[];
k1=0;
while 1
    k1=k1+1;
    ind1_1=G(:,1)==v & m_c_valid;
    ind1_2=G(:,2)==v & m_c_valid;
    ind1=find(ind1_1 | ind1_2);
    if ~isempty(ind1)
        ind2=randi(length(ind1));
        
        if ind2<=sum(ind1_1)
            d=1;
        else
            d=-1;
        end
        
        if d==1
            v=G(ind1(ind2),2);
            c=G(ind1(ind2),3);
        else
            v=G(ind1(ind2),1);
            c=G(ind1(ind2),3);
        end
        
        v_h=[v_h,v];
        c_h=[c_h,c];
        d_h=[d_h,d];
        ind_list=[ind_list,rows_G(ind1(ind2))];
        if length(v_h)>1 && any(v==v_h(1:end-1))
            v=v_h(end-1);
            G(ind1(ind2),:)=[];
            m_c_valid(ind1(ind2))=[];
            rows_G(ind1(ind2))=[];
            
            v_h(end)=[];
            c_h(end)=[];
            d_h(end)=[];
            ind_list(end)=[];
            
        end
        if v==v2
            ind_def=randi([1,length(v_h)]);
            c_solap=c_h(ind_def);
            d=d_h(ind_def);
            p=[v1,v_h];
            break
        end
        
    elseif ~isempty(v_h)
        v=v_h(end);
        ind_aux=rows_G==ind_list(end);
        G(ind_aux,:)=[];
        m_c_valid(ind_aux)=[];
        rows_G(ind_aux)=[];
        
        v_h(end)=[];
        c_h(end)=[];
        d_h(end)=[];
        ind_list(end)=[];
    else
        c_solap=[];
        d=[];
        p=[];
        break
    end
    
end
end

function [g, hay_ciclos, k1]=remueve_ciclos_2(g,G,T,C)
g_ori=g;
G_ori=G;
T_ori=T;
rows_G=1:size(G,1);
n_s=length(unique(g(:,1)));
ind_ini=randi(size(G,1));

v=G(ind_ini,1);
c=G(ind_ini,3);
vi_f=T(ind_ini,1);
n_c=unique(G(:,3));

v_h_f=[];
c_h_f=[];
vi_h=[];
k1=1;
ind_last_trans=[];
hay_ciclos=false;
while 1
    k1=k1+1;
    if k1==4
        1;
    end
    
    c_tran=find(C(c,:));
    ind1=find(G(:,1)==v & any(G(:,3)==c_tran,2));
    
    if ~isempty(ind1)%si puedo avanzar, avanzo
        v_h_f=[v_h_f,v];
        c_h_f=[c_h_f,c];
        vi_h=[vi_h,vi_f];
        
        ind2=randi(length(ind1));
        v=G(ind1(ind2),2);
        c=G(ind1(ind2),3);
        vi_f=T(ind1(ind2),1);
        
        ind_last_trans=[ind_last_trans,rows_G(ind1(ind2))];
        
        %si hay ciclo expando
        if sum(v_h_f==v)>0
            hay_ciclos=true;
            ciclo_v=[v_h_f(find(v_h_f==v,1,'first'):end),v];
            ciclo_c=[c_h_f(find(v_h_f==v,1,'first'):end),c];
            
            t_v=[];
            t_c=[];
            for i=1:length(ciclo_v)-1
                t_v=[t_v;ciclo_v(i:i+1)];
                t_c=[t_c;ciclo_c(i:i+1)];
            end
            
            ind_v_out=[];
            for i=2:length(ciclo_v)
                indAux1=g(:,3)==ciclo_v(i);
                if all(g(indAux1,2)==ciclo_v(i))
                    ind_v_out=[ind_v_out,i];
                end
            end
            ciclo_v(ind_v_out)=[];
            ciclo_c(ind_v_out)=[];
            
            ind_tran_viable=find(ismember(t_v(:,1),ciclo_v));
            ind_tran_viable=ind_tran_viable(randi(length(ind_tran_viable)));
            
            borde_expand=[t_v(ind_tran_viable,:),t_c(ind_tran_viable,2)];
            ind3=find(ismember(G,borde_expand,'rows'));
            T_expand=T(ind3,2);
            
            g=expand(g,T_expand);
            g=remueve_v_excedentes(g);
            
            %construyo de vuelta G y T
            [G,T]=g2D(g);
            C=eye(size(C))==1;
            
            while 1
                
                if length(C)<max(G(:,3))
                    1;
                end
                [c_solap,d]=solap_indirecto(G,T,C);
                
                [G_inv,T,C2]=aplico_solap_ind(G,T,C,c_solap,d);
                
                if all(C2(:)==C(:))
                    break
                else
                    C=C2;
                end
                
            end
            C=C2;
            G=G_inv;

            G_ori=G;
            T_ori=T;
            rows_G=1:size(G,1);
            
            ind_ini=randi(size(G,1));

            v=G(ind_ini,1);
            c=G(ind_ini,3);
            vi_f=T(ind_ini,1);
            
            v_h_f=[];
            c_h_f=[];
            vi_h=[];
            ind_last_trans=[];
            
            
        end
    elseif ~isempty(v_h_f)
        ind_aux=rows_G==ind_last_trans(end);
        G(ind_aux,:)=[];
        T(ind_aux,:)=[];
        ind_last_trans(end)=[];
        rows_G(ind_aux)=[];
        v=v_h_f(end);
        c=c_h_f(end);
        vi_f=vi_h(end);
        v_h_f(end)=[];
        c_h_f(end)=[];
        vi_h(end)=[];
        
    else
        ind_ini=randi(size(G,1));
        v=G(ind_ini,1);
        c=G(ind_ini,3);
        vi_f=T(ind_ini,1);
        v_h_f=[];
        c_h_f=[];
        vi_h=[];
        ind_last_trans=[];
    end
    
    if isempty(G)
        break
    end
    
end
end

function g=expand(g,T_expand)
V=unique(reshape(g(:,2:3),[],1));
S=unique(g(:,1));
n_s=length(S);
v_new=max(V)+1;
v_to_expand=g(T_expand(1),3);

g(T_expand,3)=v_new;

v_t_expanded=g(g(:,2)==v_to_expand,3);

g_add=[S,ones(n_s,1)*v_new,v_t_expanded];
g=[g;g_add];
end

function [g, consistencia_y, y, Zi, ZdFit, w, U, isomorf_1, isomorf_2, nCambios]=make_zd(g,data_in)
isomorf_1=false;
y=[];
Zi=[];
Zd=[];
ZdFit=[];
nCambios=0;
n_rep_1=1;
n_rep_2=4;
U=[];
rep_1=0;
cont_1_max=1000;
n_generaciones=0;

isomorf_2=false;
y=[];
Zi=[];
ZdFit=[];
nCambios=0;
n_rep_1=1;
n_rep_2=5;
U=[];
v_d=dame_unicos(g(:,3));
nFilas=size(dame_unicos(Zd,'rows'),1);
rep_1=0;

while ~isomorf_2 && rep_1<n_rep_1
    rep_1=rep_1+1;
    rep_2=0;
    
    while ~isomorf_2 && rep_2<=n_rep_2
        n_generaciones=n_generaciones+1;
        rep_2=rep_2+1;
        
        consistencia_y=true;
        isomorf_1=false;
        isomorf_2=false;
        nCambios=0;
        
        lista_pares=nchoosek(1:length(v_d),2);
        lista_pares=lista_pares(randperm(size(lista_pares,1)),:);
        lista_pares_ori=lista_pares;
        while 1
            lista_pares=lista_pares_ori;
            Zd_aux=[];
            c1=1;
            
            while size(lista_pares,1)>0
                v1=v_d(lista_pares(c1,1));
                v2=v_d(lista_pares(c1,2));
                
                
                consistencia_y_aux=false;
                cont_1=0;
                while ~consistencia_y_aux && cont_1<=cont_1_max
                    cont_1=cont_1+1;
                    [consistencia_y_aux, z_aux, v_incons, v_dif_1, v_dif_2] = make_zd_5(g, v1, v2);
                    if (~consistencia_y_aux)
                        display('inconsistencia de delta y que solucione arriba pero requiere varios intentos')
                        return
                        
                    else
                        Zd_aux=[Zd_aux,z_aux];
                        pares_dif_1=ismember(v_d(lista_pares(:,1)), v_dif_1) & ismember(v_d(lista_pares(:,2)), v_dif_2);
                        pares_dif_2=ismember(v_d(lista_pares(:,1)), v_dif_2) & ismember(v_d(lista_pares(:,2)), v_dif_1);
                        lista_pares(pares_dif_1 | pares_dif_2,:)=[];
                        lista_pares_latest=lista_pares;
                        break
                    end
                end
                if cont_1>cont_1_max
                    display('while duro mucho')
                    return
                end
                
                1;
            end
            
            if ~isempty(data_in.N_neu_min)
                nNeuPutative=size(Zd,2)+size(Zd_aux,2);
                if nNeuPutative/length(v_d)>n_generaciones
                    if nNeuPutative>data_in.N_neu_min
                        nNeuDelta=nNeuPutative-data_in.N_neu_min;
                        Zd=[Zd,Zd_aux(:,1:end-nNeuDelta)];
                        n_neu=size(Zd,2);
                        break
                    else
                        Zd=[Zd,Zd_aux];
                        n_neu=size(Zd,2);
                    end
                else
                    Zd=[Zd,Zd_aux];
                    n_neu=size(Zd,2);
                end
            else
                Zd=[Zd,Zd_aux];
                n_neu=size(Zd,2);
                if n_neu/length(v_d)>n_generaciones
                    break
                end
            end
        end
        [y, Zi, ZdFit, w, U, isomorf_1, isomorf_2, nCambios]=use_FTP(g,Zd, 'perceptron');
        if isomorf_2
            break
        end
    end
end
if ~isomorf_2
    ind_cc_neg=find(corr(Zd')==-1);
    disp(['pares zd con corr neg: ',num2str(length(ind_cc_neg))])
end
end

function [consistencia, z_aux, v_incons, v_dif_1, v_dif_2] = make_zd_5(g, v1, v2)
edges=(1:max(g(:,3))+1)-1/2;
counts=histcounts(g(:,3),edges);
[count_sort,v_d_sort]=sort(counts,'descend');
v_d_sort(count_sort==0)=[];
v_undef=setdiff(v_d_sort,[v1,v2],'stable');


z_aux=NaN(size(g,1),1);
z_v1_v2=[0;1];
if rand<0.5
   z_v1_v2=flipud(z_v1_v2);
end
z_aux(g(:,3)==v1)=z_v1_v2(1);
z_aux(g(:,3)==v2)=z_v1_v2(2);
s=dame_unicos(g(:,1));
n_s=length(s);
z_val=repmat([0;1],1,length(v_undef));

for i=1:size(z_val,2)
    if rand<0.5
        z_val(:,i)=flipud(z_val(:,i));
    end
end

z_val_ori=z_val;
v_incons=[0,0];

v_def=[];
L_u=v_undef;
L_a=[];
L_a=[L_a,L_u(1)];
z_backup=z_aux;
nan_g=isnan(g(:,1));

consistencia=true;

while ~isempty(L_a)
    
    c1=find(ismember(v_undef,L_a(end)));
    z_val_aux=z_val(:,c1);
    if isempty(c1)
        break
    end
    if ~all(isnan(z_val_aux))
        ind1=find(~isnan(z_val_aux),1);
        ind2=g(:,3)==v_undef(c1);
        z_aux_ori=z_aux;
        z_aux(ind2)=z_val_aux(ind1);
        [consistencia, t_incons, delta]=check_delta_y(g,z_aux,s);
        switch 1
            case 1  %uso deltas hasta ahora para definir lo que impliquen
                if consistencia
                    
                    z_val_aux(ind1)=NaN;
                    z_val(:,c1)=z_val_aux;
                    
                    g_aux=reshape(g(:,3),n_s,[]);
                    [fil,col]=find(triu(delta)~=0);
                    
                    for c2=1:size(fil,1)
                        z_aux_2=reshape(z_aux,n_s,[]);
                        d=delta(fil(c2),col(c2));
                        z_aux_3=z_aux_2([fil(c2),col(c2)],:);
                        g_aux_3=g_aux([fil(c2),col(c2)],:);
                        for c3=1:size(z_aux_3,2)
                            if isnan(z_aux_3(1,c3)) && z_aux_3(2,c3)==0 && d==-1
                                z_aux_3(1,c3)=0;
                                v_mod=g_aux_3(1,c3);
                                z_new=0;
                                z_aux_2([fil(c2),col(c2)],:)=z_aux_3;
                                z_aux(g(:,3)==v_mod)=z_new;
                            elseif isnan(z_aux_3(2,c3)) && z_aux_3(1,c3)==0 && d==1
                                z_aux_3(2,c3)=0;
                                v_mod=g_aux_3(2,c3);
                                z_new=0;
                                z_aux_2([fil(c2),col(c2)],:)=z_aux_3;
                                z_aux(g(:,3)==v_mod)=z_new;
                            elseif isnan(z_aux_3(1,c3)) && z_aux_3(2,c3)==1 && d==1
                                z_aux_3(1,c3)=1;
                                v_mod=g_aux_3(1,c3);
                                z_new=1;
                                z_aux_2([fil(c2),col(c2)],:)=z_aux_3;
                                z_aux(g(:,3)==v_mod)=z_new;
                            elseif isnan(z_aux_3(2,c3)) && z_aux_3(1,c3)==1 && d==-1
                                z_aux_3(2,c3)=1;
                                v_mod=g_aux_3(2,c3);
                                z_new=1;
                                z_aux_2([fil(c2),col(c2)],:)=z_aux_3;
                                z_aux(g(:,3)==v_mod)=z_new;
                            end
                            
                        end
                    end
                    [consistencia, t_incons, delta]=check_delta_y(g,z_aux,s);
                    if ~consistencia
                        z_aux=z_aux_ori;
                        if v_incons(1)<c1
                            v_incons=[c1, t_incons];
                        end
                    else
                        v_def=dame_unicos(g(~isnan(z_aux),3));
                        L_u = setdiff(v_undef, v_def, 'stable');
                        if isempty(L_u)
                            break
                        else
                            L_a = [L_a, L_u(1)];
                            L_u(1)=[];
                            z_backup=[z_backup,z_aux];
                        end
                        
                    end
                else
                    z_val_aux(ind1)=NaN;
                    z_val(:,c1)=z_val_aux;
                    z_aux=z_aux_ori;
                    if v_incons(1)<c1
                        v_incons=[c1, t_incons];
                    end
                    
                end
        end
        
    else
        L_a(end)=[];
        z_backup(:,end)=[];
        if ~isempty(L_a)
            z_aux=z_backup(:,end);
            v_def=dame_unicos(g(~isnan(z_aux),3));
            L_u = setdiff(v_undef, [v_def;L_a(end)], 'stable');
            ind = find(ismember(v_undef,L_u));
            z_val(:,ind)=z_val_ori(:,ind);
        else
            break
        end
    end
    
end
consistencia = consistencia & ~any(isnan(z_aux) & ~nan_g);
ind_1=z_aux==1;
ind_0=z_aux==0;
v_dif_1=dame_unicos(g(ind_1,3));
v_dif_2=dame_unicos(g(ind_0,3));


end

function [y, Zi, Zd_sol, wPer, U, isomorf_1, isomorf_2, success]=use_FTP(g,Zd,modo_genera_u)
v_d=dame_unicos(g(:,3));
v_i=dame_unicos(g(:,2));
v_all=dame_unicos(reshape(g(:,2:3),[],1));
S=dame_unicos(g(:,1));
n_vi=length(v_i);
n_s=length(S);

y=kron(ones(n_vi,1),eye(n_s));
Zi=NaN(size(Zd));
for i=1:length(v_all)
    ind_aux=find(g(:,3)==v_all(i),1,'first');
    ind_aux_2=g(:,2)==v_all(i);
    if ~isempty(ind_aux)
        Zi(ind_aux_2,:)=repmat(Zd(ind_aux,:),sum(ind_aux_2),1);
    else
        disp('source node is never a target node')
        break
    end
end

G=g2G(g);
ind_g_nan=isnan(g(:,1));
y=y(~ind_g_nan,:);
Zi=Zi(~ind_g_nan,:);
Zd=Zd(~ind_g_nan,:);

g2=cz2g([y,Zi],Zd, n_s);
G2=g2G(g2);
switch modo_genera_u
    case 'FTP'
        [~, filLC, rango, U, F, ZdFit, success]=gFTP([y,Zi],Zd);
    case 'perceptron'
        [wPer,~,EPer,iter_sum]=genera_w(y,Zi,Zd);
        U=[y,Zi]*wPer;
        success=EPer;
    case 'linear_programming'
        [~,U]=linear_prog(y,Zi,Zd);
        success=0;
        ZdFit=Zd;
end
isomorf_1=isisomorphic(G,G2);
if size(U,2)==size(Zd,2)
    Zd_sol=(signo(U)+1)/2;
    isomorf_2=success;
else
    isomorf_2=false;
end

if isomorf_2
    1;
end

end

function consistencia=check_y_cons(g,Zd)
n_neu=size(Zd,2);
consistencia=true;
s=dame_unicos(g(:,1));
for c1=1:n_neu
    [cons_aux, t_incons, delta]=check_delta_y(g,Zd(:,c1),s);
    if ~cons_aux
        consistencia=false;
        break
    end
    
end
end

function [consistencia, t_incons, delta]=check_delta_y(g,zd,s)
consistencia=true;
nS=length(s);
zd=reshape(zd,nS,[]);
vi_f=reshape(g(:,2),nS,[]);
vd=reshape(g(:,3),nS,[]);
delta=zeros(nS);
t_incons=0;
for c1=1:length(s)
    s1=s(c1);
    for c2=c1+1:length(s)
        s2=s(c2);
        delta_aux=zd(s1,:)-zd(s2,:);
        ind1=delta_aux==1;
        ind2=delta_aux==-1;
        if any(ind1) && any(ind2)
            fil=ind1;
        else
            fil=[];
        end
        vi_aux=vi_f(1,fil);
        if ~isempty(vi_aux)
            t_incons=find(g(:,2)==vi_aux & (g(:,1)==s1 | g(:,1)==s2),1);
        end
        op1=any(delta_aux>0) & any(delta_aux<0);
        op1_signo=sign(nansum(delta_aux));
        if ~all(isnan(delta_aux)) && ~op1
            delta(c1,c2)=op1_signo;
            delta(c2,c1)=-op1_signo;
        elseif ~all(isnan(delta_aux)) && op1
            consistencia=false;
            break
        end
    end
    if ~consistencia
        break
    end
end
end

function [w,U,E,iter]=genera_w(y,Zi,Zd)
C=[y,Zi];
Zd=(Zd-1/2)*2;
w=zeros(size(C,2),size(Zd,2));
U=zeros(size(Zd));
e=zeros(size(Zd,2),1);
iter=zeros(size(Zd,2),1);
for c1=1:size(Zd,2)
    [w(:,c1), U(:,c1), e(c1), iter(c1)]= ajusta_perceptron5(C,Zd(:,c1));
    if e(c1)>0
        break
    end
end
E=~any(e>0);
iter=sum(iter);
end

function [v, u, e, iter]=ajusta_perceptron5(C,zd)
%Accelerated perceptron (Wang 2022)
norma_C=norma(C,2);
C=C./norma_C;
A=C.*zd;

n_iter_max=size(A,1)*1000;

v=zeros(size(A,2),1);
q=ones(size(zd,1),1)/size(zd,1);
g=zeros(size(v));
iter=0;
while iter<n_iter_max
    iter=iter+1;
    theta=iter/(2*(iter+1));
    beta=iter/(iter+1);
    v=v-theta*(g-A'*q);
    aux=exp(-A*v);
    q=aux/sum(aux);
    g=beta*(g-A'*q);
    u=C*v;
    z=sign(u);
    e=sum(abs(z-zd));
    if e==0
        break
    end
end
end

function [w,U]=linear_prog(y,Zi,Zd)
A=[y,Zi];
Zd=-2*Zd+1;
options=optimoptions(@linprog, 'display', 'none');
w=zeros(size(A,2),size(Zd,2));
for c1=1:size(Zd,2)
    A_aux=A.*Zd(:,c1);
    A_aux=-A_aux;
    b=ones(size(Zd,1),1).*(-rand(size(Zd,1),1))*1;
    [w_aux, fval, exitflag] = linprog(zeros(1,size(A_aux,2)),A_aux,b,[],[],[],[],options);
    if exitflag~=1
        break
    else
        w(:,c1)=w_aux;
    end
end
U=-A*w;
end

function y=dame_unicos(varargin)
x=varargin{1};
op='flat';
if nargin==2
    op=varargin{2};
end
switch op
    case 'flat'
        y=unique(x(:));
        y(isnan(y))=[];
    case 'rows'
        y=unique(x,'rows');
        y(isnan(sum(y,2)),:)=[];
end
end

function [G,am]=g2G(g)
est_list=dame_unicos(g(:,[2,3]));
gAux = g(:,[2,3]);
gAux2 = gAux;
for i=1:length(est_list)
    gAux2(gAux==est_list(i))=i;
end
g = [g(:,1),gAux2];
nEst=max(reshape(g(:,[2,3]),[],1));
am=zeros(nEst);
idS=g(:,[1,2]);
idT=g(:,3);
for k1=1:nEst
    ind=find(idS(:,2)==k1);
    am(k1,idT(ind))=1;
end
G=digraph(am);
end

function [D,T]=g2D(g)
S=dame_unicos(g(:,1));
n_s=length(S);
v_list=dame_unicos(g(:,2));
D=[];
T=[];
k1=0;
for c1=1:n_s
    s1=S(c1);
    ind1=g(:,1)==s1;
    for c2=c1+1:n_s
        s2=S(c2);
        ind2=g(:,1)==s2;
        k1=k1+1;
        for c3=1:length(v_list)
            ind_v_i=g(:,2)==v_list(c3);
            ind3=find((ind1 | ind2) & ind_v_i);
            if length(ind3)>1
                D=[D;[g(ind3,3)'],k1];
                T=[T;[v_list(c3),ind3(1),ind3(2)]];
            end
        end
    end
end
ind_out=D(:,1)==D(:,2);
D(ind_out,:)=[];
T(ind_out,:)=[];
end

function g=cz2g(C,zTarget,nS)
idS=NaN(size(C,1),2);
idT=NaN(size(C,1),1);
zBase=dame_unicos([C(:,nS+1:end);zTarget],'rows');
y=C(:,1:nS);
yBase=unique(y,'rows');
nEst=size(zBase,1);
yID=[];
for i=1:size(yBase,1)
    ind=find(ismember(y,yBase(i,:),'rows'));
    yID(ind,1)=i*ones(length(ind),1);
end
for k1=1:nEst
    ind1=find(ismember(zTarget,zBase(k1,:),'rows'));
    ind2=find(ismember(C(:,nS+1:end),zBase(k1,:),'rows'));
    idT(ind1,1)=k1;
    idS(ind2,:)=[yID(ind2),k1*ones(length(ind2),1)];
end
g=[idS,idT];
end

function g = remueve_v_excedentes(g)
salida=false;
while ~salida
    g_ori=g;
    g = remueve_v_raiz(g);
    g = remueve_v_hojas(g);
    salida = size(g_ori,1)==size(g,1) && (all(g(:)==g_ori(:) | (isnan(g(:)) | isnan(g(:)))));
end
end

function g=remueve_v_hojas(g)
vi=unique(g(:,2));
vd=unique(g(:,3));
vd_no_vi=setdiff(vd,vi);
vd_no_vi=vd_no_vi(~isnan(vd_no_vi));
for c1=1:length(vd_no_vi)
    ind_out=g(:,3)==vd_no_vi(c1);
    g(ind_out,:)=[];
end
end

function g=remueve_v_raiz(g)
vi=unique(g(:,2));
vd=unique(g(:,3));
vi_no_vd=setdiff(vi,vd);
vi_no_vd=vi_no_vd(~isnan(vi_no_vd));
for c1=1:length(vi_no_vd)
    ind_out=g(:,2)==vi_no_vd(c1);
    g(ind_out,:)=[];
end
end