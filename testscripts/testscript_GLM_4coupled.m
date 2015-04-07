% testscript_GLM_4coupled.m
%
% Test code for simulating and fitting a coupled GLM (4 neurons).
%
% Notes:
%   Fitting code uses same functions as for single-cell responses.
%   Simulation code requires new structures / functions
%     (due to the need to pass activity between neurons)
%
% Code Blocks:  
%   1. Set up model params and plot
%   2. Show samples from simulated model
%   3. Generate training data
%   4. Fit simulated dataset w/ maximum-likelihood (ML)
%   5. Plot results

global RefreshRate;  % Stimulus refresh rate (Stim frames per second)
RefreshRate = 100;


%%  1.  Set parameters and display for GLM  ============ %

DTsim = .01; % Bin size for simulating model & computing likelihood (in units of stimulus frames)
nkt = 40;    % Number of time bins in filter;
ttk = [-nkt+1:0]';
ggsimsolo = cell(4,1);
for whichn = 1:4
	ggsimsolo{whichn} = makeSimStruct_GLM(nkt,DTsim);  % Create GLM struct with default params
end;
k = ggsimsolo{1}.k;  % Stimulus filter
ggsimsolo{2}.k = shift(k,-3)*1.2;
ggsimsolo{3}.k = shift(k,-5)*1.3;
ggsimsolo{4}.k = shift(k,-7)*1.4;
ggsim = makeSimStruct_GLMcpl(ggsimsolo{1},ggsimsolo{2},ggsimsolo{3},ggsimsolo{4});

% Make some coupling kernels
[iht,ihbas,ihbasis] = makeBasis_PostSpike(ggsim.ihbasprs,DTsim);
hhcpl = ihbasis*[.6;.47;.25;0;0]*2;
hhcpl(:,2) = ihbasis*[-1;-1;0;0;.25]*2;
% NOTE: here we are setting our simulation up as TWO coupled pairs, with the pairs independent of each other
ggsim.ih(:,2,1) = hhcpl(:,2); % 2nd cell coupling to first
ggsim.ih(:,1,2) = hhcpl(:,1); % 1st cell coupling to second
ggsim.ih(:,4,3) = hhcpl(:,2); % 4th to 3rd
ggsim.ih(:,3,4) = hhcpl(:,1); % 3rd to 4th

% === Make Fig: model params =======================
figure(1);
clf;
plotcols = {'r', 'b', 'g', 'c'};
for whichn = 1:4
	subplot(4,2,whichn*2-1);
	plot(ttk, ggsim.k(:,:,whichn));
	title(sprintf('neuron %i stim kernel', whichn));
end;
xlabel('time (frames)');

for whichn = 1:4
	subplot(4,2,whichn*2); % --------
	hold on;
	for fromn = 1:4
		plot(ggsim.iht, exp(ggsim.ih(:,fromn,whichn)), plotcols{fromn});
	end
	title(sprintf('spike kernels into Cell %i (exponentiated)', whichn));
	legend('from 1', 'from 2', 'from 3', 'from 4', 'location', 'northeast');
	ylabel('multiplicative gain');
end;
xlabel('time after spike (frames)');


%% 2. Make GWN stimulus & simulate the glm model response. ========= %
% 
slen = 50; % Stimulus length (frames) & width (# pixels)
swid = size(ggsim.k,2);
Stim = 2*randn(slen,swid);  % Gaussian white noise stimulus
[tsp, vmem, Ispk] = simGLM(ggsim, Stim); % Simulate GLM response

% ==== Make Figure ========
figure(2);
clf;
tt = [DTsim:DTsim:slen]';
subplot(5,2,1); %------------------------
plot(1:slen, Stim, 'k', 'linewidth', 2); 
title('GWN stimulus');
axis tight;
for whichn = 4:-1:1
	printf("Plotting currents for cell %i. minvmem is %g, maxvmem is %g, size(tsp) is %s\n", whichn, min(vmem(:,whichn)), max(vmem(:,whichn)), mat2str(size(tsp{whichn})));
	subplot(5,2,whichn*2+1); %------------------------
	plot(tt, vmem(:,whichn), tsp{whichn}, max(vmem(:,whichn))*ones(size(tsp{whichn})), 'ro');
	title(sprintf('cell %i: net voltage + spikes', whichn));
	axis tight;
	subplot(5,2,whichn*2+2)
	plot(tt, vmem(:,whichn)-Ispk(:,whichn), 'k', tt, Ispk(:,whichn), 'r');
	title('stim-induced & spike-induced currents'); 
	axis tight;
	xlabel('time (frames)');
end;



%% 3. Generate some training data
slen = 2500;  % Stimulus length (frames);  More samples gives better fit
Stim = round(rand(slen,swid))*4-2;  %  Run model on long, binary stimulus
[tsp,vmem,ispk] = simGLM(ggsim,Stim);  % run model

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
	subplot(4,2,whichn*2-1);  % Filters cell 1 % ---------------
	plot(ttk, ggsimsolo{whichn}.k, 'k', ttk, gg{whichn}.k, 'r');
	title(sprintf('Cell %i stim filt (True=blck, ML=red)', whichn));
end;
xlabel('time (frames)')



for whichn = 1:4
	subplot(4,2, whichn*2); % ----------------------------------
	hold on;
	for fromn = 1:4
		if whichn==fromn
			ihdata = gg{whichn}.ih;
		elseif whichn<fromn
			ihdata = gg{whichn}.ih2(:, fromn-1);
		else
			ihdata = gg{whichn}.ih2(:, fromn);
		end
		plot(ggsim.iht, exp(ggsim.ih(:,fromn,whichn)), sprintf('%s--', plotcols{fromn}), gg{whichn}.iht, exp(gg{whichn}.ihbas*ihdata), plotcols{fromn});
	end;
	title(sprintf('Cell %i: exponentiated post-spk kernels', whichn));
	axis tight;
end;
xlabel('time (frames)')

