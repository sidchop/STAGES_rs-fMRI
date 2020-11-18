pcs = importdata("a2_illness_cca_pcs.csv") ; 
clin = importdata("a2_illness_cca_clin1.csv") ;

%optional PCA
%[~,PCscores] = pca(zscore(clin.data));

% [wcoeff,score,latent,tsquared,explained] = pca(zscore(clin.data)) ;
%pareto(explained)

 clinPC = PCscores(:,1:8);   
 clinPC = clin.data(:,1:5); 

[pfwer,r,A,B,U,V] = permcca(zscore(pcs.data),clin.data,1000)

%[pfwer,r,A,B,U,V] = permcca(pcs.data, clin.data,5000)

corr(V(:,1), PCscores)

scatter(U(:,1), V(:,1))

corr(U(:,1),PCscores)