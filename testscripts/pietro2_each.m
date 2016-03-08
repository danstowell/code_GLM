%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
use_daysgrps = { ...
'd3_S', 'd5_S', 'd7_S', ...
'd3_T', 'd5_T', 'd7_T', ...
};
use_vers = {'ind'};

for whichver=use_vers
for whichdaytrio=use_daysgrps
	csvpath = sprintf('~/datasets/pietro_mpio_duos/translated/pietro2_%s_%s.csv', whichdaytrio{1}, whichver{1});
	disp(sprintf('Fitting %s', csvpath));
	if whichver{1}=='ind'
		k = 2;
	end
	for startsecs=0:3600:28800-1
		endsecs   = startsecs + 3600;
		runlabel = sprintf('%s_h%i_%s', whichdaytrio{1}, startsecs/3600, whichver{1});
		disp(sprintf('   run %s', runlabel));
		regln = -1; % NOTE default regularisation strength here
		resimuldur = 0;
		[numcalls, peakpos, peakval, neglogli, dc] = dofit_fromcsv_GLM_zf4f(csvpath, runlabel, k, [1:k], startsecs, endsecs, regln, 'pietro2/outplot', 'pietro2/outcsv', resimuldur);
	end
end
end

