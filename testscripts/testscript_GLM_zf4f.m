function [numcalls, peakpos, peakval] = testscript_GLM_zf4f(dataname, variantname, indexmapper, startsecs, endsecs)
% [peakpos, peakval] = testscript_GLM_zf4f(dataname, variantname, indexmapper, startsecs, endsecs)
%
% load some zf4f data and analyse "as if" it were cell spiking data. writes out a plot.

global RefreshRate;
RefreshRate = 1;  % the "RefreshRate" is the samplerate of the stimulus (in Hz). I don't currently use stimulus so I set it to 1. Below, "DTsim" sets the time-resultion used in the model.

plotcols = {'r', 'b', 'g', 'm'};

numcalls = zeros(4,1);
peakpos = zeros(4);
peakval = zeros(4);

printf('testscript_GLM_zf4f(%s, %s, %s, %i, %i)\n', dataname, variantname, mat2str(indexmapper), startsecs, endsecs);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load the CSV data
events = csvread(sprintf('~/git/stored_docs/python/zftranscribe/output/annotreconciledproofed/zcompiled_%s.csv', dataname));
%size(events)
events = events((events(:,1) >= startsecs) & (events(:,1) <= endsecs), [1,3]);
%size(events)

tsp = cell(1,4);
for whichn=1:4
	matchid = indexmapper(whichn)-1;
	tsp{whichn} = (events(events(:,2)==matchid, 1) - startsecs) * RefreshRate;   % subtracting startsecs, and converting to units of RefreshRate
	numcalls(whichn) = size(tsp{whichn}, 1);
	printf('Bird %i has %i events\n', whichn, numcalls(whichn));
end;
fflush(stdout); % NB octave-only

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%  1.  Set parameters

DTsim = .005; % Bin size for simulating model & computing likelihood (in units of stimulus frames)
nkt = 40;    % Number of time bins in filter;
ttk = [-nkt+1:0]';
ggsimsolo = cell(4,1);
for whichn = 1:4
	ggsimsolo{whichn} = makeSimStruct_GLM(nkt,DTsim);  % Create GLM struct with default params
end;
ggsim = makeSimStruct_GLMcpl(ggsimsolo{1},ggsimsolo{2},ggsimsolo{3},ggsimsolo{4});


%% 3. Set up the "stimulus" appropriately (here it's zeros)
slen = round((endsecs-startsecs) * RefreshRate);  % Stimulus length (frames)
swid = 1;
Stim = zeros(slen,swid);





% -------------- Compute STAs------------
nsp = length(tsp{1});
stas = cell(4,1);
for whichn = 1:4
	sta0 = simpleSTC(Stim,tsp{whichn},nkt);
	stas{whichn} = reshape(sta0,nkt,[]); 
end;


%% 4. Do ML fitting of GLM params
gg = cell(4, 1);
opts = {'display', 'iter', 'maxiter', 100};
for whichn = 1:4
	%fprintf('Fitting bird #%i\n', whichn);
	gg0 = makeFittingStruct_GLM(stas{whichn},DTsim,ggsim,whichn);  % Initialize params for fitting struct w/ sta
% TODO: the call above imposes an "absref" (absolute refactory period) of 10 * DTsim. I'd quite like to be able to shrink that away.
	gg0.ih = gg0.ih*0;  % Initialize to zero
	gg0.dc = gg0.dc*0;  % Initialize to zero

	gg0.tsp = tsp{whichn};   % cell spike times (vector)
	gg0.tsp2 = tsp(setdiff(1:4, whichn));  % spike trains from "coupled" cells (cell array of vectors)
	gg0.tspi = 1; % 1st spike to use for computing likelihood (eg, can ignore 1st n spikes)
	[gg{whichn}, negloglivalx] = MLfit_GLM(gg0, Stim, opts); % do ML (requires optimization toolbox)
	printf('MLfit #%i gets neglogli %g\n', whichn, negloglivalx);
end


%% --- Plot results ----------------------------
h = figure(3);
clf;

set (h,'papertype', '<custom>');
set (h,'paperunits','inches');
set (h,'papersize',[6 5]);
set (h,'paperposition', [0,0,[6 5]]);
set (0,'defaulttextfontsize', 10);
set (0,'defaultaxesfontsize', 10);

ttk = -nkt+1:0;

for whichn = 1:4
	subplot(2,2, whichn); % ----------------------------------
	hold on;
	for fromn = 1:4
		plotcol = plotcols{fromn};
		if whichn==fromn
			ihdata = gg{whichn}.ih;
			plotcol = 'k--';
		elseif whichn<fromn
			ihdata = gg{whichn}.ih2(:, fromn-1);
		else
			ihdata = gg{whichn}.ih2(:, fromn);
		end
		plotx = gg{whichn}.iht;
		ploty = exp(gg{whichn}.ihbas*ihdata);
		plot(plotx, ploty, plotcol);
		ylim([0, 5]);

		% For collating, we'll calc the peak pos&sizes
		[peakvalraw, peakposraw] = max(ploty);
		peakpos(fromn,whichn) = plotx(peakposraw) / RefreshRate;
		peakval(fromn,whichn) = peakvalraw;
	end;
	title(sprintf('Bird %i: exp(kernels) %s%s', whichn, dataname, variantname));
	legend('from 1', 'from 2', 'from 3', 'from 4', 'location', 'northeast');
	axis tight;
end;
xlabel('time (frames)')

%saveas(h, 'zf4f_glm_kernels_session2a1.pdf')
%print(h, sprintf('zf4f_glm_kernels_%s%s.eps', dataname, variantname), '-depsc')
saveas(h, sprintf('zf4f_glm_kernels_%s%s.png', dataname, variantname))

