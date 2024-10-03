function M = myM(node,k,hK,c,ndim)
% 计算scaled monomials
% c: center

x = node(:,1);
y = node(:,2);
xK = c(1);
yK = c(2);


% M = zeros(nnode,)

if k == 1
    M = [1.0+0*x, (x-xK)./hK, (y-yK)./hK];
elseif k == 2
    M = [1.0+0*x, (x-xK)./hK, (y-yK)./hK, (x-xK).^2/hK^2, (x-xK).*(y-yK)./hK^2, (y-yK).^2./hK^2];
elseif k == 3
    M = [1.0+0*x, (x-xK)./hK, (y-yK)./hK, (x-xK).^2/hK^2, (x-xK).*(y-yK)./hK^2, (y-yK).^2./hK^2, ...
        (x-xK).^3/hK^3,(x-xK).^2.*(y-yK)/hK^3,(x-xK).*(y-yK).^2/hK^3,(y-yK).^3/hK^3];
end