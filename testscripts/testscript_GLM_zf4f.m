function [numcalls, peakpos, peakval, neglogli, dc] = testscript_GLM_zf4f(dataname, variantname, indexmapper, startsecs, endsecs, nlfun, resimuldur)
% [numcalls, peakpos, peakval, neglogli, dc] = testscript_GLM_zf4f(dataname, variantname, indexmapper, startsecs, endsecs)
%
% load some zf4f data and analyse "as if" it were cell spiking data. writes out a plot.

if nargin < 6
	nlfun = @softplus
end
if nargin < 7
	resimuldur = 0
end

printf('testscript_GLM_zf4f(%s, %s, %s, %i, %i, %s)\n', dataname, variantname, mat2str(indexmapper), startsecs, endsecs, func2str(nlfun));

runlabel = sprintf('%s%s%s', dataname, variantname, func2str(nlfun)(1:3));
csvpath = sprintf('~/git/stored_docs/python/zftranscribe/output/annotreconciledproofed/zcompiled_%s.csv', dataname);

disp(sprintf('Fitting with nonlin %s on %s', func2str(nlfun), csvpath));

k = length(indexmapper);
regln = -1; % NOTE default regularisation strength here
[numcalls, peakpos, peakval, neglogli, dc] = dofit_fromcsv_GLM_zf4f(csvpath, runlabel, k, indexmapper, startsecs, endsecs, regln, 'outplot', 'outcsv', resimuldur, nlfun);

