%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
use_daysgrps = { ...
'15.4', ...
};
use_vers = {'ind'};

for whichver=use_vers
for whichdaytrio=use_daysgrps
	csvpath = sprintf('~/datasets/pietro_mpio_tetra/translated/pietro4_%s_%s.csv', whichdaytrio{1}, whichver{1});
	disp(sprintf('Fitting %s', csvpath));
	if whichver{1}=='ind'
		k = 4;
	end
	for startsecs=0:3600:32400-1
		endsecs   = startsecs + 3600;
		runlabel = sprintf('%s_h%i_%s', whichdaytrio{1}, startsecs/3600, whichver{1});
		disp(sprintf('   run %s', runlabel));
		regln = -1; % NOTE default regularisation strength here
		resimuldur = 0;
		[numcalls, peakpos, peakval, neglogli, dc] = dofit_fromcsv_GLM_zf4f(csvpath, runlabel, k, [1:k], startsecs, endsecs, regln, 'pietro4/outplot', 'pietro4/outcsv', resimuldur);
	end
end
end

