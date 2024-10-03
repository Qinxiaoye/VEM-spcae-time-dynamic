function [Pi,Pis] = calculatePi(elemID,centroid,diameter,Area,node,elem,ndim,k)
% calculate elem stiffness matrix
% for k = 1
% input: elemID,centroid,diameter,node,elem,ndim,k,elem2edge,edge
% output: AK


index = elem{elemID};
Nv = length(index);
hK = diameter(elemID);

x = node(index ,1); y = node(index ,2);

% calculate D
D = myM([x,y],k,hK,centroid(elemID,:),ndim);

% calculate B
Gradm = [0 0; 1./hK*[1, 0]; 1./hK*[0, 1]]; % k = 1
rotid1 = [Nv,1:Nv-1]; rotid2 = [2:Nv,1]; % ending and starting indices
normVec = 0.5*[y(rotid2)-y(rotid1), x(rotid1)-x(rotid2)]'; % a rotation of edge vector
B = Gradm*normVec; % B, Eq.(69)

% constraint
Bs = B;  Bs(1,:) = 1/Nv; % constraint条件，公式24
% consistency relation
G = B*D;  Gs = Bs*D; % Eqs.(55),(70)

% --------- local stiffness matrix ---------
Pis = Gs\Bs;  % Eq.(71)
Pi  = D*Pis;
% I = eye(size(Pi));
% AK  = Pis'*G*Pis + (I-Pi)'*(I-Pi);% Eq.(72)