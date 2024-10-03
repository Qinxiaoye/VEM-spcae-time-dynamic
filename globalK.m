function [GK] = globalK(node,elem,mat,zline,elem3,theta)
% calculate the global stiffness matrix
% input: node,elem,k,ndim
% output: GK

% calculate the geometric quantity for each element
aux = auxgeometry(node,elem);% 计算每个单元的形心，面积，半径，具体什么原理不用管
centroid = aux.centroid;  
diameter = aux.diameter;  
area = aux.area;


sumElem = size(elem,1); % the number of element
sumNode = size(node,1);


elemLen = cellfun('length',elem);
elemLenNew = repmat(elemLen*2,1,length(zline)-1)';
elemLen = elemLenNew(:);
nnz = sum((elemLen*4).^2);
NsumNode = sumNode*length(zline);


ii = zeros(nnz,1); jj = zeros(nnz,1); ss = zeros(nnz,1); 

% calculate element stiffness matrix
Pi_all = cell(sumElem,1);
Pis_all = cell(sumElem,1);
for n = 1:sumElem
    [Pi,Pis] = calculatePi(n,centroid,diameter,area,node,elem,2,1);
    Pi_all{n} = Pi;
    Pis_all{n} = Pis;
end

ia = 0;
for n = 1:sumElem
    Pi = Pi_all{n}; Pis = Pis_all{n};
    for m = 1:(length(zline)-1)
        AK = elemK(n,node,elem,area(n),centroid(n,:),diameter(n),Pi,Pis,mat,[zline(m),zline(m+1)],theta);
        AB = reshape(AK,1,[]); 
        index = [elem{n}+sumNode*(m-1),elem{n}+sumNode*m];
        Nv = length(index);
        elemDof = [index, index+NsumNode, index+NsumNode*2, index+NsumNode*3];
        % --------- assembly index for ellptic projection -----------
        indexDof = elemDof;  Ndof = Nv*4;  % local to global index
        ii(ia+1:ia+Ndof^2) = repmat(indexDof(:), Ndof, 1);
        jj(ia+1:ia+Ndof^2) = reshape(repmat(indexDof, Ndof, 1), [], 1);
        ss(ia+1:ia+Ndof^2) = AB(:);
        ia = ia + Ndof^2;
    end
end

GK = sparse(ii,jj,ss,NsumNode*4,NsumNode*4);