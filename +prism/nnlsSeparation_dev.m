function sources = nnlsSeparation_dev(mixing, data)

addpath(genpath('C:\Users\Federico\Documents\MATLAB\nnlslab'));

opt_gen1 = initopt_general('perfmeas', {1,2}, 'maxcpu', 5, 'stopcrit', 2, 'tol', 1E-12);
pg_1 = projgradient(mixing, data', opt_gen1, opt_projgradient());
sources = pg_1.xopt;

end