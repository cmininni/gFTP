%Check that each node v that appears in G is encoded by one unique vector
%z.
%Outputs "true" if this happens, false otherwise.
function isomorf=check_dynamics(G,Y,Z_s,Z_t)
isomorf=true;
V=unique(reshape(G(:,2:3),[],1));
S=unique(G(:,1));
Z_pool=[];
for c1=1:length(V)
    ind1=G(:,2)==V(c1);
    ind2=G(:,3)==V(c1);
    Z=[Z_s(ind1,:);Z_t(ind2,:)];
    aux=size(unique(Z,'rows'),1);
    if aux>1
        isomorf=false;
        break
    end
    Z_pool=[Z_pool;Z(1,:)];
end

Z_pool=unique(Z_pool,'rows');
if size(Z_pool,1)~= length(V)
    isomorf=false;
end

Y_pool=[];
for c1=1:length(S)
    ind1=find(G(:,1)==S(c1));
    aux=size(unique(Y(ind1,:),'rows'),1);
    if aux>1
        isomorf=false;
        break
    end
    Y_pool=[Y_pool;Y(ind1(1),:)];
end

Y_pool=unique(Y_pool,'rows');
if size(Y_pool,1)~=length(S)
    isomorf=false;
end
