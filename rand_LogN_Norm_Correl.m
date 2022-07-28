function X = rand_LogN_Norm_Correl(meanLogN, stdLogN, meanNorm, stdNorm, correl,n)

muconv = @(m,v) log(m/sqrt(1+v/m^2)); %convert mean and std of logN in its 1st parameter
sigmaconv = @(m,v) sqrt(log(1+v/m^2)); %convert mean and std of logN in its 2sd parameter

Z = mvnrnd([0 0], [1 correl; correl 1], n);
U = normcdf(Z,0,1);
X = [logninv(U(:,1),muconv(meanLogN,stdLogN),sigmaconv(meanLogN,stdLogN)) norminv(U(:,2),meanNorm,stdNorm)];
end