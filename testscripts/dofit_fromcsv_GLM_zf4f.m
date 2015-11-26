function [numcalls, peakpos, peakval, neglogli, dc] = dofit_fromcsv_GLM_zf4f(csvpath, runlabel, k, indexmapper, startsecs, endsecs, regln, plotpath, csvoutpath, resimuldur, nlfun)
% [numcalls, peakpos, peakval, neglogli, dc] = dofit_fromcsv_GLM_zf4f(csvpath, runlabel, k, indexmapper, startsecs, endsecs, regln, plotpath, csvoutpath, resimuldur, nlfun)
%
% load some zf4f-format data and analyse "as if" it were cell spiking data. returns analysed data.
% also does a plot and writes it to a file in the folder named by 'plotpath'. if 'plotpath' is '' or 0 it DOESN'T plot. to plot in cwd use '.'
% 'csvoutpath' parameter is analogous, and is about writing CSV data out to file
% 'regln' is regularisation strength. use 0 for no regln (ML rather than MAP), or maybe a val like 0.1, or -1 for default setting.
% 'resimuldur' is 0 if you don't want to re-simulate from the fitted model; otherwise the num seconds worth of data to synthesise

global RefreshRate;
RefreshRate = 1;  % the "RefreshRate" is the samplerate of the stimulus (in Hz). I don't currently use stimulus so I set it to 1. Below, "DTsim" sets the time-resultion used in the model.

plotcols = {'r', 'b', 'g', 'm', 'y'};

numcalls = zeros(k,1);

printf('dofit_fromcsv_GLM_zf4f(%s, %s, %i, %s, %i, %i, %g, %s, %s, %i)\n', csvpath, runlabel, k, mat2str(indexmapper), startsecs, endsecs, regln, plotpath, csvoutpath, resimuldur);

if regln == -1
	regln = 0.01;
end

if nargin < 11
	nlfun = @softplus;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load the CSV data
events = csvread(csvpath);
%size(events)
events = events((events(:,1) >= startsecs) & (events(:,1) <= endsecs), [1,3]);
%size(events)

tsp = cell(1,k);
for whichn=1:k
	matchid = indexmapper(whichn)-1;
	tsp{whichn} = (events(events(:,2)==matchid, 1) - startsecs) * RefreshRate;   % subtracting startsecs, and converting to units of RefreshRate
	numcalls(whichn) = size(tsp{whichn}, 1);
	printf('Bird %i has %i events\n', whichn, numcalls(whichn));
end;
fflush(stdout); % NB octave-only

endsecs_actual = max(events(:,1)) + 0.1;   % this tells us how far into the distance we actually need to look.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%  1.  Set parameters

DTsim = .005; % Bin size for simulating model & computing likelihood (in units of stimulus frames)
nkt = 40;    % Number of time bins in filter;
ttk = [-nkt+1:0]';
ggsimsolo = cell(k,1);
for whichn = 1:k
	ggsimsolo{whichn} = makeSimStruct_GLM(nkt,DTsim);  % Create GLM struct with default params
end;
ggsim = makeSimStruct_GLMcpl(ggsimsolo{1:k});
ggsim.nlfun = nlfun;


%% 3. Set up the "stimulus" appropriately (here it's zeros)
slen = round((endsecs_actual-startsecs) * RefreshRate);  % Stimulus length (frames)
if slen==0
	error('Stim is zero length - no data in selected time window?');
end
swid = 1;
Stim = zeros(slen,swid);




% -------------- Compute STAs------------
nsp = length(tsp{1});
stas = cell(k,1);
for whichn = 1:k
	sta0 = simpleSTC(Stim,tsp{whichn},nkt);
	stas{whichn} = reshape(sta0,nkt,1);
	%stas{whichn} = reshape(sta0,nkt,[]);
end;


%% 4. Do ML fitting of GLM params
gg = cell(k, 1);
neglogli = 0;
opts = {'display', 'iter', 'maxiter', 100};
for whichn = 1:k
	%fprintf('Fitting bird #%i\n', whichn);
	gg0 = makeFittingStruct_GLM(stas{whichn},DTsim,ggsim,whichn);  % Initialize params for fitting struct w/ sta
	gg0.ih = gg0.ih*0;  % Initialize to zero
	gg0.dc = gg0.dc*0;  % Initialize to zero
	gg0.nlfun = nlfun;

	gg0.tsp = tsp{whichn};   % cell spike times (vector)
	gg0.tsp2 = tsp(setdiff(1:k, whichn));  % spike trains from "coupled" cells (cell array of vectors)
	gg0.tspi = 1; % 1st spike to use for computing likelihood (eg, can ignore 1st n spikes)
	[gg{whichn}, neglogli_each] = MLfit_GLM(gg0, Stim, regln, opts); % do ML (requires optimization toolbox)
	printf('MLfit #%i gets neglogli %g\n', whichn, neglogli_each);
	neglogli += neglogli_each;
end


%% --- Calc summary stats - used for csv and for returning ----------------------------
% the init of 1.23456789 is to flag any errors that may creep in with failure to fill values in. a blue plaster.
peakpos = ones(k) * 1.23456789;
peakval = ones(k) * 1.23456789;
dc = ones(k,1) * 1.23456789;
plotx = gg{1}.iht;
kernels_discret = ones(k,k,size(plotx)) * 1.23456789;
for whichn = 1:k
	for fromn = 1:k
		if whichn==fromn
			ihdata = gg{whichn}.ih;
		elseif whichn<fromn
			ihdata = gg{whichn}.ih2(:, fromn-1);
		else
			ihdata = gg{whichn}.ih2(:, fromn);
		end
		ploty = gg{whichn}.ihbas*ihdata;
		if fromn==whichn
			[ignore_this, peakposraw] = max(-ploty);  % negative so that for self-self it's the inhibitions we're looking at
		else
			[ignore_this, peakposraw] = max(ploty);
		endif
		peakpos(fromn,whichn) = plotx(peakposraw) / RefreshRate;
		peakval(fromn,whichn) = ploty(peakposraw); % re-grab the peakval, because now we don't want the abs
		kernels_discret(fromn,whichn,:) = ploty;
	end
	dc(whichn) = gg{whichn}.dc;
end

%% --- Plot results ----------------------------
if plotpath
	h = figure();
	clf;

	set (h,'papertype', '<custom>');
	set (h,'paperunits','inches');
	set (h,'papersize',[6 5]);
	set (h,'paperposition', [0,0,[6 5]]);
	set (0,'defaulttextfontsize', 10);
	set (0,'defaultaxesfontsize', 10);

	legendargs = cell(k+2,1);
	for whichn = 1:k
		legendargs{whichn} = sprintf('to %i', whichn);
	end
	legendargs{k+1} = 'location';
	legendargs{k+2} = 'northeast';

	ttk = -nkt+1:0;

	numrows = ceil(sqrt(k));
	numcols = ceil(k/numrows);
	ymin = -0.5;
	ymax =  1;
	for fromn = 1:k
		subplot(numrows, numcols, fromn); % ----------------------------------
		hold on;
		for whichn = 1:k
			plotcol = plotcols{mod(whichn-1, numel(plotcols))+1};
			if whichn==fromn
				plotcol = 'k--';
			else
				% we only allow self-other kernels to stretch the y-axis
				ymin = min(ymin, min(kernels_discret(fromn,whichn,:)));
				ymax = max(ymax, max(kernels_discret(fromn,whichn,:)));
			end
			plot(plotx, kernels_discret(fromn,whichn,:), plotcol);
		end;
		title(sprintf('Bird %i: kernels (%s) %s', fromn, func2str(nlfun), runlabel));
		legend(legendargs{1:k+2});
		axis tight;
		ylim([ymin, ymax]);
	end;
	xlabel('time (frames)')

	disp(sprintf('Saving %s/zf4f_glm_kernels_%s.png', plotpath, runlabel));
	saveas(h, sprintf('%s/zf4f_glm_kernels_%s.png', plotpath, runlabel));

	sleep(1);
	close();
else
	disp '  (not plotting)';
end



if csvoutpath
	outfnamestem = sprintf('%s/zf4f_glm_stats_%s', csvoutpath, runlabel);
	csvfp_0d = fopen(sprintf('%s_0d.csv', outfnamestem), 'w');
	csvfp_1d = fopen(sprintf('%s_1d.csv', outfnamestem), 'w');
	csvfp_2d = fopen(sprintf('%s_2d.csv', outfnamestem), 'w');
	csvfp_tx = fopen(sprintf('%s_kernels_timeaxis.csv', outfnamestem), 'w');
	csvfp_kd = fopen(sprintf('%s_kernels_discret.csv', outfnamestem), 'w');
	csvfp_ba = fopen(sprintf('%s_basis.csv', outfnamestem), 'w');
	csvfp_bb = fopen(sprintf('%s_basis_orth.csv', outfnamestem), 'w');
	% headers
	fprintf(csvfp_0d, 'runname,neglogli\n');
	fprintf(csvfp_1d, 'runname,individ,numcalls,dc\n');
	fprintf(csvfp_2d, 'runname,frm,too,peakval,peakpos\n');
	% data
	fprintf(csvfp_0d, '%s,%g\n', runlabel, neglogli);
	for whichn = 1:k
		fprintf(csvfp_1d, '%s,%i,%i,%g\n', runlabel, whichn, numcalls(whichn), dc(whichn));
		for fromn = 1:k
			fprintf(csvfp_2d, '%s,%i,%i,%g,%g\n', runlabel, fromn, whichn, peakval(fromn,whichn), peakpos(fromn,whichn));

			fprintf(csvfp_kd, '%s,%i,%i', runlabel, fromn, whichn);
			for xpos=1:size(plotx)
				fprintf(csvfp_kd, ',%g', kernels_discret(fromn,whichn,xpos));
			end
			fprintf(csvfp_kd, '\n');
		end
	end
	fprintf(csvfp_tx, '%g', plotx(1));
	for xpos=2:size(plotx)
		fprintf(csvfp_tx, ',%g', plotx(xpos));
	end
	fprintf(csvfp_tx, '\n');
	% Here we're getting an explicit copy of the non-orthogonalised basis, using the "ihbasprs" parameters that were used to initialise the fit
	[iht,ihbas,ihbasis] = makeBasis_PostSpike(gg{1}.ihbasprs,DTsim);
	for whichbas=1:size(ihbasis, 2)
		fprintf(        csvfp_ba, '%g',  ihbasis(1   , whichbas));
		fprintf(        csvfp_bb, '%g',  ihbas(  1   , whichbas));
		for xpos=2:size(plotx)
			fprintf(csvfp_ba, ',%g', ihbasis(xpos, whichbas));
			fprintf(csvfp_bb, ',%g', ihbas(xpos, whichbas));
		end
		fprintf(csvfp_ba, '\n');
		fprintf(csvfp_bb, '\n');
	end
	fflush(csvfp_0d);
	fclose(csvfp_0d);
	fflush(csvfp_1d);
	fclose(csvfp_1d);
	fflush(csvfp_2d);
	fclose(csvfp_2d);
	fflush(csvfp_kd);
	fclose(csvfp_kd);
	fflush(csvfp_tx);
	fclose(csvfp_tx);
	fflush(csvfp_ba);
	fclose(csvfp_ba);
	fflush(csvfp_bb);
	fclose(csvfp_bb);

	if resimuldur > 0
		disp("*** NOTE: resimulation from fitted model. not validated yet.");

		swid = 1;  % Stimulus width  (pixels).  Must match # pixels in stim filter
		Stim = zeros(resimuldur * RefreshRate,swid);

		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% full-group resimulation

		% copying fitted parameters into simulation struct
		%disp('dc pre');
		%disp(ggsim.dc);
		for whichn = 1:k
			ggsim.dc(whichn) = gg{whichn}.dc;
		end
		%disp('dc fitted');
		%disp(ggsim.dc);
		for whichn = 1:k
			for fromn = 1:k
				%disp(size(ggsim.ih(:,fromn,whichn)));
				%disp(size(gg{whichn}.ih));
				%disp(size(gg{whichn}.ih2));
				
				if whichn==fromn
					kernelprs = gg{whichn}.ih;
				elseif whichn < fromn
					kernelprs = gg{whichn}.ih2(:,fromn-1);
				else
					kernelprs = gg{whichn}.ih2(:,fromn);
				end
				kernel = gg{whichn}.ihbas * kernelprs;
				ggsim.ih(:,fromn,whichn) = kernel;
			end
		end

		[tsp, vmem, Ispk] = simGLMcpl(ggsim, Stim);  % Simulate GLM response
		% since this is a multi thing, it returns "tsp" as a cell array of 3 items, each one being a list of timestamps

		% convert the tsp to a sortable list
		for whichtsp=1:size(tsp, 2)
			atsp = tsp{whichtsp};
			atsp = horzcat(atsp, repmat(whichtsp, size(atsp), 1));
			if whichtsp==1
				tspdata = atsp;
			else
				tspdata = vertcat(tspdata, atsp);
			end;
		end;

		tspdata = sortrows(tspdata);

		% write out a CSV file of the generated timestamps
		csvfp = fopen(sprintf('%s_resimulated.csv', outfnamestem), 'w');
		fprintf(csvfp, 'time,dursecs,individ\n');
		for row=tspdata'
			fprintf(csvfp, '%g,%g,%i\n', row(1), 0.1, row(2)-1);
		end
		fclose(csvfp);

		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% single-solo-individual resimulation

		% create a brand new gg object
		ggsim_just_solo = makeSimStruct_GLM(nkt,DTsim);
		% set gg parameters: nlfun, dc, ih
		whichn = 1
		ggsim_just_solo.nlfun = nlfun;
		ggsim_just_solo.dc = gg{whichn}.dc;
		ggsim_just_solo.ih = gg{whichn}.ihbas * gg{whichn}.ih;
		% call the main simul thing
		[tsp, vmem, Ispk] = simGLM(ggsim_just_solo, Stim);
		% write out a CSV file of the generated timestamps
		csvfp = fopen(sprintf('%s_resimulated_asolo.csv', outfnamestem), 'w');
		fprintf(csvfp, 'time,dursecs,individ\n');
		for atimepos=tsp'
			fprintf(csvfp, '%g,%g,%i\n', atimepos, 0.1, whichn-1);
		end
		fclose(csvfp);

	end

	sleep(2.5);
else
	disp '  (not writing csv)';
end

