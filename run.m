warning off

load YALE_165n_1024d_15c_zscore_uni.mat
load YALE_165n_1024d_15c_zscore_uni_allkernel.mat

alpha=1e-5;
beta=15;
gamma=1;

result=selfweightmkl(K,y,alpha,beta,gamma)

