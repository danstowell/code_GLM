% testscript_GLM_zf4f.m
%
% load some zf4f data and analyse "as if" it were cell spiking data

global RefreshRate;
RefreshRate = 2;

plotcols = {'r', 'b', 'g', 'm'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load the CSV data
events = csvread('~/git/stored_docs/python/zftranscribe/output/annotreconciledproofed/zcompiled_session2a.csv');
starttime = 300;
endtime = 1200;
%size(events)
events = events((events(:,1) >= starttime) & (events(:,1) <= endtime), [1,3]);
%size(events)

tsp = cell(1,4);
for whichn=1:4
	tsp{whichn} = (events(events(:,2)==whichn-1, 1) - starttime) * RefreshRate;   % subtracting starttime, and converting to units of RefreshRate
	printf('Bird %i has %i events\n', whichn, size(tsp{whichn}, 1));
end;
fflush(stdout); % NB octave-only

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%  1.  Set parameters

DTsim = .01; % Bin size for simulating model & computing likelihood (in units of stimulus frames)
nkt = 40;    % Number of time bins in filter;
ttk = [-nkt+1:0]';
ggsimsolo = cell(4,1);
for whichn = 1:4
	ggsimsolo{whichn} = makeSimStruct_GLM(nkt,DTsim);  % Create GLM struct with default params
end;
ggsim = makeSimStruct_GLMcpl(ggsimsolo{1},ggsimsolo{2},ggsimsolo{3},ggsimsolo{4});


%% 3. Set up the "stimulus" appropriately (here it's zeros)
slen = round((endtime-starttime) * RefreshRate);  % Stimulus length (frames)
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
	fprintf('Fitting neuron %i\n', whichn);
	gg0 = makeFittingStruct_GLM(stas{whichn},DTsim,ggsim,whichn);  % Initialize params for fitting struct w/ sta
	gg0.ih = gg0.ih*0;  % Initialize to zero
	gg0.dc = gg0.dc*0;  % Initialize to zero

	gg0.tsp = tsp{whichn};   % cell spike times (vector)
	gg0.tsp2 = tsp(setdiff(1:4, whichn));  % spike trains from "coupled" cells (cell array of vectors)
	gg0.tspi = 1; % 1st spike to use for computing likelihood (eg, can ignore 1st n spikes)
	[gg{whichn}, negloglivalx] = MLfit_GLM(gg0, Stim, opts); % do ML (requires optimization toolbox)
end


%% --- Plot results ----------------------------
figure(3);
clf;
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
		plot(gg{whichn}.iht, exp(gg{whichn}.ihbas*ihdata), plotcol);
		ylim([0, 5]);
	end;
	title(sprintf('Bird %i: exponentiated post-call kernels', whichn));
	legend('from 1', 'from 2', 'from 3', 'from 4', 'location', 'northeast');
	axis tight;
end;
xlabel('time (frames)')

