function [GKNew,fNew] = boudnaryCondition(GK,fixNode,f,sumNode)

nodeNof = unique([fixNode;fixNode+sumNode;fixNode+sumNode*2;fixNode+sumNode*3]);

GKNew = GK;
GKNew(nodeNof,:) = 0;

GKNew = GKNew+sparse(nodeNof,nodeNof,ones(length(nodeNof),1),size(GK,1),size(GK,1));

f(nodeNof) = 0;
fNew = sparse(f);