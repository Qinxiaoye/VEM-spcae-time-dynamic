function ff = getForce(nodeNew,sumNode,sumNodeNew,elem,auxT,k,press,zline,direction,theta)

sumT = length(zline);


ff = zeros(sumNodeNew*2,1);
sumP = size(press,1);

nip = k+2;
[x,w] = gaussInt(nip);

for n = 1:sumP
    p = zeros(k+1,1);
    elemID = press(n,1);
    faceID = press(n,2);
    value = press(n,3);
    index = elem{elemID};
    indexEdge = auxT.elem2edge{elemID};
    Nv = length(index);
    v1 = 1:Nv; v2 = [2:Nv,1]; % loop index for vertices or edges
    
    if k == 1
        elem1 = [v1(:), v2(:)];
        faceNode = elem1(faceID,:);
        faceNodeID = index(faceNode);
        edgeNodeCoor = nodeNew(faceNodeID,:);
    elseif k == 2
        elem1 = [v1(:), v2(:), v1(:)+Nv];
        faceNode = elem1(faceID,:);
        midNode = sumNode+indexEdge(faceID);
        faceNodeID = [index([faceNode(1:2)]),midNode];
        edgeNodeCoor = nodeNew(faceNodeID,:);
    end
    
    L = nodeNew(faceNodeID(1:2),:);
    L = (L(1,:)-L(2,:));
    Normal = [L(2),-L(1)]/norm(L);
    for m = 1:nip
        N = fun(x(m),0,0,1,k+1);
        dN = dfun(x(m),0,0,1,k+1);
        if k == 2
            N = [N(1);N(3);N(2)];
            dN = [dN(1);dN(3);dN(2)];
        end
        J = sqrt((dN'*edgeNodeCoor(:,1))^2+(dN'*edgeNodeCoor(:,2))^2);
        p = p+w(m)*N*value*J;
    end
    
    % integral along time
    for t = 1:sumT-1
        h = zline(t+1)-zline(t);
        T = [0;0];
        intDf = [0;0];
        for m = 1:nip
            N = fun(x(m),0,0,1,2);
            J = h/2;
            T = T+w(m)*N*J;
            df = dzshape(N'*[zline(t),zline(t+1)]',zline(t),zline(t+1));
            intDf = intDf+w(m)*df'*J;
        end
        P = [p*T(1);p*T(2)]+theta*h*[p*intDf(1);p*intDf(2)];
    
        faceNodeIDNew = [faceNodeID+sumNode*(t-1),faceNodeID+sumNode*t];

        if strcmp(direction,'normal')
            ff(faceNodeIDNew) = ff(faceNodeIDNew)+P*Normal(1);
            ff(faceNodeIDNew+sumNodeNew) = ff(faceNodeIDNew+sumNodeNew)+P*Normal(2);
        elseif strcmp(direction,'x')
            ff(faceNodeIDNew) = ff(faceNodeIDNew)+P;
        elseif strcmp(direction,'y')
            ff(faceNodeIDNew+sumNodeNew) = ff(faceNodeIDNew+sumNodeNew)+P;
        end
    end
end
ff = sparse([zeros(sumNodeNew*2,1);ff]);

function df = dzshape(t,t1,t2)

df1 = -1/(t2-t1);
df2 = 1/(t2-t1);

df = [df1,df2];