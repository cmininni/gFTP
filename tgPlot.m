% Plots labeled multidigraph G. 
% It uses distinguishable_colors.m (Tim Holy (2024). Generate maximally
% perceptually-distinct colors
% (https://www.mathworks.com/matlabcentral/fileexchange/29702-generate-maximally-perceptually-distinct-colors),
% MATLAB Central File Exchange)
function tgPlot(varargin)
layout='force';
if nargin==3
    C=varargin{1};
    zTarget=varargin{2};
    lY=varargin{3};
    g=cz2g(C,zTarget,lY);
elseif nargin==1
    g=varargin{1};
elseif nargin==2
    g=varargin{1};
    layout=varargin{2};
end
S=dame_unicos(g(:,1));
nS=max(S);
g(isnan(sum(g,2)),:)=[];
g=gSinGaps(g);
%[G,am]=g2G(g);
e1=g(:,2);
e2=g(:,3);
l=(1:size(g,1));
G=digraph(e1,e2,l);
colorSet=distinguishable_colors(nS);
t_flecha=10;    %16 ; 20
t_linea=2;      %2  ; 3
t_nodo=5;       %5  : 8
h=plot(G,'Layout',layout,'LineWidth',t_linea,'ArrowSize',t_flecha);
h.MarkerSize=t_nodo;
for i=1:nS
    t=find(g(:,1)==S(i));
    ind=find(ismember(G.Edges.Weight,t));
    highlight(h,'Edges',ind,'EdgeColor',colorSet(i,:))
    annotation('textbox',[0.02,(nS-i+1)*0.1,0.3,0.3],'String',num2str(i),'FitBoxToText','on','Color',colorSet(i,:));
end
end
