function [numcalls, peakpos, peakval, neglogli] = testscript_GLM_zf4f(dataname, variantname, indexmapper, startsecs, endsecs, nonlin)
% [peakpos, peakval] = testscript_GLM_zf4f(dataname, variantname, indexmapper, startsecs, endsecs)
%
% load some zf4f data and analyse "as if" it were cell spiking data. writes out a plot.

if nargin < 6
	nonlin = 'exp'
end

printf('testscript_GLM_zf4f(%s, %s, %s, %i, %i)\n', dataname, variantname, mat2str(indexmapper), startsecs, endsecs);

runlabel = sprintf('%s%s', dataname, variantname);
if nonlin == 'exp'
	nlfun = @expfun;
elseif nonlin == 'sof'
	nlfun = @softplus;
	runlabel = sprintf('%s%s', runlabel, nonlin);
else
	error('Unknown nonlin choice "%s"', nonlin);
end
csvpath = sprintf('~/git/stored_docs/python/zftranscribe/output/annotreconciledproofed/zcompiled_%s.csv', dataname);

disp(sprintf('Fitting with nonlin %s on %s', nonlin, csvpath));

k = 4;
regln = -1; % NOTE default regularisation strength here
[numcalls, peakpos, peakval, neglogli] = dofit_fromcsv_GLM_zf4f(csvpath, runlabel, k, indexmapper, startsecs, endsecs, regln, 'outplot', 'outcsv', 0, nlfun);
% TODO maybe plotpath and csvoutpath as args?
% TODO maybe allow this to do the resynth? maybe as an arg?

