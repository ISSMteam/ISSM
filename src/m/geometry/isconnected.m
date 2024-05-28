function flag=isconnected(elements,A,B)
%ISCONNECTED: are two nodes connected by a triangulation?
%
%   Usage: flag=isconnected(elements,A,B)
%
%

elements=ElementsFromEdge(elements,A,B);
if isempty(elements),
	flag=0;
else
	flag=1;
end
