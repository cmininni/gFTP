%Construct auxiliary graph D from graph G
function [D,T]=make_D(G)
S=dame_unicos(G(:,1));
n_s=length(S);
v_list=dame_unicos(G(:,2));
D=[];
T=[];
k1=0;
for c1=1:n_s
    s1=S(c1);
    ind1=G(:,1)==s1;
    for c2=c1+1:n_s
        s2=S(c2);
        ind2=G(:,1)==s2;
        k1=k1+1;
        for c3=1:length(v_list)
            ind_v_i=G(:,2)==v_list(c3);
            ind3=find((ind1 | ind2) & ind_v_i);
            
            D=[D;[G(ind3,3)'],k1];
            T=[T;[v_list(c3),ind3(1),ind3(2)]];
        end
        
    end
end
ind_out=D(:,1)==D(:,2);
D(ind_out,:)=[];
T(ind_out,:)=[];
end
