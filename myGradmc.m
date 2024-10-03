function dM = myGradmc(node,k,hK,c,ndim)

x = node(:,1);
y = node(:,2);
xK = c(1);
yK = c(2);

if ndim == 3
    z = node(:,3);
    zK = c(3);
end

if k == 1
    dM = [0,0; 1./hK*[1, 0]; 1./hK*[0, 1]];
elseif k == 2
    dM = [[0+0*x, 0+0*x];[1+0*x, 0+0*x]./hK;[0+0*x, 1+0*x]./hK;[2*(x-xK), 0+0*x]./hK^2;[(y-yK), (x-xK)]./hK^2;[0+0*x, 2*(y-yK)]./hK^2];
elseif k == 3
    dM = [[0+0*x, 0+0*x];[1+0*x, 0+0*x]./hK;[0+0*x, 1+0*x]./hK;[2*(x-xK), 0+0*x]./hK^2;[(y-yK), (x-xK)]./hK^2;[0+0*x, 2*(y-yK)]./hK^2;...
        [3*(x-xK).^2/hK^3, 0+0*x];[2*(x-xK).*(y-yK)/hK^3, (x-xK).^2/hK^3];[(y-yK).^2/hK^3, 2*(x-xK).*(y-yK)/hK^3];[0+0*x, 3*(y-yK).^2/hK^3]];
end

