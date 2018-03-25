load YALE_165n_1024d_15c_zscore_uni.mat
load YALE_165n_1024d_15c_zscore_uni_allkernel.mat

alpha=1e-5;
beta=25;
mu=.1;

r=0.1;%rate of labeled data

[m,n,rr]=size(K);
c=length(unique(y)); % number of class
numperc=floor(n/c); % number of data per class
labelperc=floor(r*numperc); % number of labeled data per class
labelindperc=sort(randperm(numperc,labelperc)); % index of labeled data selected
labelind=[]; % labelind: index of known label
for i=1:c
    labelind=[labelind labelindperc+(i-1)*numperc];
end

[result]=selfweightmklsemi(K,y,labelind,alpha,beta,mu)
