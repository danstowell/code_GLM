%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%use_days = [1, 2, 4, 5];
%use_trios = {'A', 'B', 'C', 'D', 'E'};
use_daystrios = { ...
'd1_tA', 'd2_tA', 'd4_tA', 'd5_tA',  ...
'd1_tB', 'd2_tB', 'd5_tB', 'd6_tB',  ...
'd1_tC', 'd2_tC', 'd4_tC', 'd5_tC',  ...
'd1_tD', 'd2_tD', 'd4_tD', 'd7_tD',  ...
'd1_tE', 'd2_tE', 'd4_tE', 'd5_tE',  ...
};
use_vers = {'ind', 'typ'};

for whichver=use_vers
%for whichtrio=use_trios
%for whichday=use_days
for whichdaytrio=use_daystrios

	runlabel = sprintf('%s_%s', whichdaytrio{1}, whichver{1});
	csvpath = sprintf('~/datasets/nico_mpio_trios/translated/nico3_%s_%s.csv', whichdaytrio{1}, whichver{1});

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

%end
end
end

