clear;
zline = 0:0.1:10;
load poly.mat
[node3,elem3] = PolyMesh3Simple(node,elem,zline);

theta = 0.01;

sumNode = size(node,1);
mat = [200000,0.3,10]; % 弹性模量、泊松比, 密度
% global stiffness matrix

[GK] = globalK(node,elem,mat,zline,elem3,theta);

leftNode = find(node3(:,1)<0.01);
bottomNode = find(node3(:,3)<0.01);
rightNode = find(node3(:,1)>4.99);

NsumNode = size(GK,1)/4;

pface = findFace(node,elem,rightNode);
p = -100;
press = [pface,ones(length(pface),1)*p];

auxT = auxstructure(node,elem);

f = getForce(node,sumNode,size(node3,1),elem,auxT,1,press,zline,'y',theta);


[GK,f] = boudnaryCondition(GK,[leftNode;bottomNode],f,NsumNode);

u = GK\f;
u = full(u);

velo = reshape(u(1:NsumNode*2),[],2)/mat(3);
vx = reshape(velo(:,1),sumNode,[]);

u = u(NsumNode*2+1:end);
u = reshape(u,[],2);
node3New = node3;
node3New(:,1:2) = node3New(:,1:2)+u;

figure;
showsolution3D(node3New(:,[1,3,2]),elem3,u(:,2))
view(43,24)

figure;
showsolution3D(node3New(:,[1,3,2]),elem3,velo(:,2))
view(43,24)

% quad 26
% poly 188
% nonconvex 210
dispNode = 188:sumNode:size(node3,1);
figure;
plot(zline,u(dispNode,2));

