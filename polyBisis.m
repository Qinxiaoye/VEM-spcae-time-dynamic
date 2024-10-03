function M = polyBisis(node,hK,c,k,ndim)

if ndim == 2
    x = node(:,1);
    y = node(:,2);
    xK = c(1);
    yK = c(2);
    if k == 1
        M = [1.0+0*x, (x-xK)./hK, (y-yK)./hK];
    elseif k == 2
        M = [1.0+0*x, (x-xK)./hK, (y-yK)./hK, (x-xK).^2/hK^2, (x-xK).*(y-yK)./hK^2, (y-yK).^2./hK^2];
    elseif k == 3
        M = [1.0+0*x, (x-xK)./hK, (y-yK)./hK, (x-xK).^2/hK^2, (x-xK).*(y-yK)./hK^2, (y-yK).^2./hK^2, ...
            (x-xK).^3/hK^3,(x-xK).^2.*(y-yK)/hK^3,(x-xK).*(y-yK).^2/hK^3,(y-yK).^3/hK^3];
    elseif k == 4
        M = [1.0+0*x, (x - xK)./hK, (y - yK)./hK, (x - xK).^2./hK^2, ((x - xK)*(y - yK))./hK^2, (y - yK).^2./hK^2, (x - xK).^3./hK^3, ((x - xK).^2*(y - yK))./hK^3, ((x - xK)*(y - yK).^2)./hK^3, (y - yK).^3./hK^3, (x - xK).^4./hK^4, ((x - xK).^3*(y - yK))./hK^4, ((x - xK).^2*(y - yK).^2)./hK^4, ((x - xK)*(y - yK).^3)./hK^4, (y - yK).^4./hK^4];
    end
else
    x = node(:,1);
    y = node(:,2);
    z = node(:,3);
    xK = c(1);
    yK = c(2);
    zK = c(3);
    
    xi = (x-xK)./hK;
    eta = (y-yK)./hK;
    zeta = (z-zK)./hK;

    if k == 1
        M = [1.0+0*x, xi,eta,zeta];
    elseif k == 2
        M = [1.0+0*x,xi,eta,zeta,xi.*eta,eta.*zeta,zeta.*xi,xi.^2,eta.^2,zeta.^2];
    elseif k == 3
        M = [1.0+0*x,xi,eta,zeta,xi.*eta,eta.*zeta,zeta.*xi,xi.^2,eta.^2,zeta.^2,...
            xi.^2.*eta,xi.^2.*zeta,eta.^2.*xi,eta.^2.*zeta,zeta.^2.*xi,zeta.^2.*eta,...
            xi.*eta.*zeta,xi.^3,eta.^3,zeta.^3];
    end

end
