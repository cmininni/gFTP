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