function N = dfun(ks,yt,zita,ndim,mnode)
% 四节点四边形的形函数
N = zeros(mnode,1);
if ndim == 1
    N(1) = -1/2;
    N(2) = 1/2;
end