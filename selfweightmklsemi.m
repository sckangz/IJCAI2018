function [result]=selfweightmklsemi(H,s,labelind,alpha,beta,gamma)
    % s is the true class label.
    warning off
    [m,n,nn]=size(H);
    sumK=zeros(n);
    for i=1:nn
        sumK=sumK+H(:,:,i);
    end
    K=sumK/nn;
    Z=eye(n);
    %Z=rand(n);
    c=length(unique(s));
    ulabel= setdiff(1:n,labelind);%index of unlabeled data

    ll=length(labelind);
    Fl=zeros(ll,c);
    for i=1:ll
        Fl(i,s(labelind(i)))=1;
    end
        Fu=zeros(n-ll,c);
    for i=1:30
        Zold=Z;

        sumH=zeros(n);
        sum1=0;
        for j=1:nn
            temp(j)=norm(H(:,:,j)-K,'fro');
            omega(j)=1/(2*temp(j));
            sumH=sumH+omega(j)*H(:,:,j);
            sum1=sum1+omega(j);
        end
        K=(2*Z'-Z*Z'-eye(n)+2*beta*sumH)/(2*beta*sum1);
        K(find(K<0))=0;
        K=(K'+K)/2;
        F=[Fl;Fu];
        parfor ij=1:n
        %         for ji=1:n
        %             all(ji)=norm(F(ij,:)-F(ji,:));
        %         end
            [all]=veccomp2(ij,n,F);
            Z(:,ij)=(K+gamma*eye(n))\(K(ij,:)'-alpha/4*all');
        end
        Z(find(Z<0))=0;
        Z= (Z+Z')/2;

        D = diag(sum(Z));
        L = D-Z;
        uu=zeros(n-ll,n-ll);
        for ii=1:(n-ll)
            for jj=1:(n-ll)
                uu(ii,jj)=L(ulabel(ii),ulabel(jj));
            end
        end
        ul=zeros(n-ll,ll);
        for ii=1:(n-ll)
            for jj=1:ll
                ul(ii,jj)=L(ulabel(ii),labelind(jj));
            end
        end
        Fu=-uu\(ul*Fl);
        
        if i>5 &((norm(Z-Zold)/norm(Zold))<1e-5)
                   break
        end

    end
    [ur,uc]=size(ulabel);
    [max_value,max_ind] = max(Fu,[],2);
    cnt = 0;
    for i = 1:uc
        if max_ind(i) == s(ulabel(i))
            cnt = cnt+1;
        end
    end
    result = cnt/uc;
end

function [all]=veccomp2(ij,n,F)
    warning off
    for ji=1:n
        all(ji)=norm(F(ij,:)-F(ji,:));
    end
end