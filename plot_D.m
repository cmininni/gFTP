% Generates graph D and plots it.
% It Uses distinguishable_colors.m (Tim Holy (2024). Generate maximally
% perceptually-distinct colors
% (https://www.mathworks.com/matlabcentral/fileexchange/29702-generate-maximally-perceptually-distinct-colors),
% MATLAB Central File Exchange)
function plot_D(varargin)
D=varargin{1};
if nargin==2
    layout=varargin{2};
else
    layout=[];
end
ind_self_loops=D(:,1)==D(:,2);
D(ind_self_loops,:)=[];
S=dame_unicos(D(:,3));
nS=max(S);
D(isnan(sum(D,2)),:)=[];
e1=D(:,1);
e2=D(:,2);
nombres=1:max([e1;e2]);
l=(1:size(D,1));
G=digraph(e1,e2,l);
colorSet=distinguishable_colors(nS);
t_flecha=20;
t_linea=3;
t_nodo=7;
if isempty(layout)
    h=plot(G,'LineWidth',t_linea,'ArrowSize',t_flecha,'NodeLabel',nombres);
else
    h=plot(G,'layout',layout,'LineWidth',2,'ArrowSize',11,'MarkerSize',10,'NodeLabel',nombres);
end
h.MarkerSize=t_nodo;

for i=1:nS
    t=find(D(:,3)==S(i));
    ind=find(ismember(G.Edges.Weight,t));
    highlight(h,'Edges',ind,'EdgeColor',colorSet(i,:))
    annotation('textbox',[0.02,(nS-i+1)*0.1,0.3,0.3],'String',num2str(i),'FitBoxToText','on','Color',colorSet(i,:));
end
1;






