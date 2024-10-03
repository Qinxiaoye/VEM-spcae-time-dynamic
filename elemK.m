function AK = elemK(elemID,node,elem,Area,centroid,hK,Pi,Pis,mat,zline,theta)
% calculate elem stiffness matrix
% for k = 1
% input: elemID,centroid,diameter,node,elem,ndim,k,elem2edge,edge
% output: AK


index = elem{elemID};

x = node(index ,1); y = node(index ,2);

EX = mat(1); mu = mat(2); rho = mat(3);
h = zline(2)-zline(1);
D = EX/(1-mu^2)*[1,mu,0;mu,1,0;0,0,(1-mu)/2];


Nv = length(x);
nodeT = [node(index,:);centroid]; % triangulation of K
elemT = [(Nv+1)*ones(Nv,1),(1:Nv)',[2:Nv,1]'];

Int = zeros(Nv);
[lambda,weight] = quadpts(2);
NT = size(elemT,1);
for iel = 1:NT
    vT = nodeT(elemT(iel,:),:);
    area = 0.5*abs(det([[1;1;1],vT]));
    xy = lambda*vT;
    for p = 1:size(xy,1)
        M = polyBisis(xy(p,:),hK,centroid,1,2);
        Pim = M*Pis;
        Int = Int + area*weight(p)*(Pim'*Pim);
    end
end
I = eye(size(Pi,1));
Int = Int+Area*(I-Pi)'*(I-Pi);

nip = 2;
[point,weight] = gaussInt(nip);
NN = zeros(2);
NdN = zeros(2);
dNN = zeros(2);
dNdN = zeros(2);
for n = 1:nip
    f = fun(point(n),0,0,1,2);
    df = dzshape(f'*zline',zline(1),zline(2));
    NN = NN+weight(n)*(f*f')*(zline(2)-zline(1))/2;
    NdN = NdN+weight(n)*(f*df)*(zline(2)-zline(1))/2;
    dNN = dNN+weight(n)*(df'*f')*(zline(2)-zline(1))/2;
    dNdN = dNdN+weight(n)*(df'*df)*(zline(2)-zline(1))/2;
end

Kpp_c = [Int*NN(1,1),Int*NN(1,2);
         Int*NN(2,1),Int*NN(2,2)]+...
         theta*h*[Int*dNN(1,1),Int*dNN(1,2);
         Int*dNN(2,1),Int*dNN(2,2)];
Kpu_c = [Int*NdN(1,1),Int*NdN(1,2);
         Int*NdN(2,1),Int*NdN(2,2)]+...
    theta*h*[Int*dNdN(1,1),Int*dNdN(1,2);
    Int*dNdN(2,1),Int*dNdN(2,2)];
Kup_c = [Int*NdN(1,1),Int*NdN(1,2);
         Int*NdN(2,1),Int*NdN(2,2)]+...
    theta*h*[Int*dNdN(1,1),Int*dNdN(1,2);
    Int*dNdN(2,1),Int*dNdN(2,2)];
Kpp_c = blkdiag(Kpp_c,Kpp_c);
Kpu_c = blkdiag(Kpu_c,Kpu_c);
Kup_c = blkdiag(Kup_c,Kup_c);

dM = myGradmc(centroid,1,hK,centroid,2);
A1 = [1,0;0,0;0,1];
A2 = [0,0;0,1;1,0];
A = [A1,A2];
dM1 = blkdiag(dM',dM');
G0 = Area*dM1'*A'*D*A*dM1;
Pis = blkdiag(Pis,Pis);
Pi = blkdiag(Pi,Pi);
int_uu_c = Pis'*G0*Pis;
int_uu = int_uu_c+0.15*trace(int_uu_c)*(eye(size(Pi,1))-Pi)'*(eye(size(Pi,1))-Pi);
int_uu11 = int_uu(1:Nv,1:Nv);
int_uu12 = int_uu(1:Nv,Nv+1:2*Nv);
int_uu21 = int_uu(Nv+1:2*Nv,1:Nv);
int_uu22 = int_uu(Nv+1:2*Nv,Nv+1:2*Nv);
Kuu_c = [NN(1,1)*int_uu11,NN(2,1)*int_uu11,NN(1,1)*int_uu12,NN(2,1)*int_uu12;
    NN(1,2)*int_uu11,NN(2,2)*int_uu11,NN(1,2)*int_uu12,NN(2,2)*int_uu12;
    NN(1,1)*int_uu21,NN(2,1)*int_uu21,NN(1,1)*int_uu22,NN(2,1)*int_uu22;
    NN(1,2)*int_uu21,NN(2,2)*int_uu21,NN(1,2)*int_uu22,NN(2,2)*int_uu22]+...
    theta*h*[dNN(1,1)*int_uu11,dNN(1,2)*int_uu11,dNN(1,1)*int_uu12,dNN(1,2)*int_uu12;
            dNN(2,1)*int_uu11,dNN(2,2)*int_uu11,dNN(2,1)*int_uu12,dNN(2,2)*int_uu12;
            dNN(1,1)*int_uu21,dNN(1,2)*int_uu21,dNN(1,1)*int_uu22,dNN(1,2)*int_uu22;
            dNN(2,1)*int_uu21,dNN(2,2)*int_uu21,dNN(2,1)*int_uu22,dNN(2,2)*int_uu22];

AK = [-1/rho*Kpp_c,Kpu_c;Kup_c,Kuu_c];


function f = zshape(t,t1,t2)
f1 = (t2-t)/(t2-t1);
f2 = (t-t1)/(t2-t1);
f = [f1,f2];

function df = dzshape(t,t1,t2)

df1 = -1/(t2-t1);
df2 = 1/(t2-t1);

df = [df1,df2];

