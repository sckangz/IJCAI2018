function [result]=selfweightmkl(H,s,alpha,beta,gamma)
% s is the true class label.
warning off
[m,n,nn]=size(H);
sumK=zeros(n);
for i=1:nn
    sumK=sumK+H(:,:,i);
end
K=sumK/nn;
Z=eye(n);
c=length(unique(s));
itt = 0;
%options = optimset( 'Algorithm','interior-point-convex','Display','off');
for i=1:100
    itt = i;
    Zold=Z;

    D = diag(sum(Z));
    L = D-Z;
    
    [F, temp, ev]=eig1(L, c, 0);
    sumH=zeros(n);
    sum1=0;
    for j=1:nn
        temp(j)=norm(H(:,:,j)-K,'fro');
        omega(j)=1/(2*temp(j));
        sumH=sumH+omega(j)*H(:,:,j);
        sum1=sum1+omega(j);
    end

    for ij=1:n
        for ji=1:n
            all(ji)=norm(F(ij,:)-F(ji,:));
        end

        Z(:,ij)=(K+gamma*eye(n))\(K(ij,:)'-alpha/4*all');
        % we use the free package to solve quadratic equation: http://sigpromu.org/quadprog/index.html
        %    [Z(:,ij),err,lm] = qpas(H,ff,[],[],ones(1,n),1,zeros(n,1),ones(n,1));
        % Z(:,ij)=quadprog(H,(beta/2*all'-(alpha-2)*K(:,ij))',[],[],ones(1,n),1,zeros(n,1),ones(n,1),Z(:,ij),options);
    end
    Z(find(Z<0))=0;
    Z= (Z+Z')/2;
    %      K=(Z'+beta*sumH)/(Z+beta*sum1*eye(n));
    K=(2*Z'-Z*Z'-eye(n)+2*beta*sumH)/(2*beta*sum1);
    K(find(K<0))=0;
    K=(K'+K)/2;
    if i>5 &((norm(Z-Zold)/norm(Zold))<1e-5)
        break
    end
end
actual_ids= kmeans(F, c, 'emptyaction', 'singleton', 'replicates', 100, 'display', 'off');
[result] = [ClusteringMeasure( actual_ids,s) itt];
