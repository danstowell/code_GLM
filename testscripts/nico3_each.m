%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
use_days = [1, 2, 4, 5];
use_trios = {'A', 'B', 'C', 'D', 'E'};
use_vers = {'ind', 'typ'};

for whichver=use_vers
for whichtrio=use_trios
for whichday=use_days

	runlabel = sprintf('d%i_t%s_%s', whichday, whichtrio{1}, whichver{1});
	csvpath = sprintf('~/datasets/nico_mpio_trios/translated/nico3_d%i_t%s_%s.csv', whichday, whichtrio{1}, whichver{1});

	disp(sprintf('Fitting %s', csvpath));

	if whichver{1}=='ind'
		k = 3;
	else
		k = 18;
	end
	startsecs = 0;
	endsecs   = 99999;
	regln = -1; % NOTE default regularisation strength here
	resimuldur = 0;
	[numcalls, peakpos, peakval, neglogli, dc] = dofit_fromcsv_GLM_zf4f(csvpath, runlabel, k, [1:k], startsecs, endsecs, regln, 'nico3/outplot', 'nico3/outcsv', resimuldur);

end
end
end

