function g=cz2g(C,zTarget,nS)
idT=[];
idS=[];
zBase=unique([C(:,nS+1:end);zTarget],'rows');
y=C(:,1:nS);
yBase=unique(y,'rows');
nEst=size(zBase,1);
yID=[];
for i=1:size(yBase,1)
    ind=find(ismember(y,yBase(i,:),'rows'));
    yID(ind,1)=i*ones(length(ind),1);
end
for k1=1:nEst
    ind1=find(ismember(zTarget,zBase(k1,:),'rows'));%filas donde aparece un zk1 en zTarget
    ind2=find(ismember(C(:,nS+1:end),zBase(k1,:),'rows'));%filas donde aparece un zkq en C
    idT(ind1,1)=k1;
    idS(ind2,:)=[yID(ind2),k1*ones(length(ind2),1)];
end
g=[idS,idT];