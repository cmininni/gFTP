%Rename nodes starting from 1, avoiding gaps
function gAux=gSinGaps(g)
   gAux=g;
   V=unique(reshape(gAux(:,[2,3]),[],1));
   for c1=1:length(V)
       ind=g(:,2)==V(c1);
       gAux(ind,2)=c1;
       ind=g(:,3)==V(c1);
       gAux(ind,3)=c1;
   end
end