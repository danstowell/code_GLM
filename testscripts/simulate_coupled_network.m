
% simulate from a coupled model, pure A->B->C structure

global RefreshRate;  % Stimulus refresh rate (Stim frames per second)
RefreshRate = 1;

DTsim = .005; % Bin size for simulating model & computing likelihood (in units of stimulus frames)
nkt = 20;    % Number of time bins in stimulus filter
ttk = [-nkt+1:0]';  % time relative to spike of stim filter taps


ggsim1 = makeSimStruct_GLM(nkt,DTsim); % Create GLM structure with default params
ggsim2 = makeSimStruct_GLM(nkt,DTsim); % Create GLM structure with default params
ggsim3 = makeSimStruct_GLM(nkt,DTsim); % Create GLM structure with default params
ggsimcpl = makeSimStruct_GLMcpl(ggsim1, ggsim2, ggsim3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NOTE: at this point you modify ggsim.ih to have the kernels you want. it's [timepos x individ x individ], and each kernel is originally made by a method like:
%    [iht,ihbas,ihbasis] = makeBasis_PostSpike(ihbasprs,dt);
%    ih = ihbasis*[-10 -5 0 2 -2]';  % h current
% A simple way to inspect the dependence structure quickly is:
%    mean(abs(ggsim.ih))
% And to plot one of the kernels:
%    plot(ggsim.iht, ggsim.ih[:,1,1])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Let's construct an A->B->C influence structure
ggsim.ih *= 0;
[iht,ihbas,ihbasis] = makeBasis_PostSpike(ggsim.ihbasprs, ggsim.dt);
toneighbour = ihbasis*[0 0 0 1 0.5]';
%plot(ggsim.iht, toneighbour);
ggsim.ih(:,1,2) = toneighbour;
ggsim.ih(:,2,3) = toneighbour;


% OK, now we can sample from the model.
slen = 60; % Stimulus length (frames) 
swid = 1;  % Stimulus width  (pixels).  Must match # pixels in stim filter
Stim = zeros(slen,swid);
[tsp, vmem, Ispk] = simGLMcpl(ggsimcpl, Stim);  % Simulate GLM response
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


% write out a CSV file
csvfp = fopen('data_simulate_coupled_network.csv', 'w');
fprintf(csvfp, 'time,individ,dursecs\n');
for row=tspdata'
	fprintf(csvfp, '%g,%i,%g\n', row(1), row(2), 0.1);
end
fclose(csvfp);

