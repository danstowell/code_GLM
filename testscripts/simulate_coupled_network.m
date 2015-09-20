
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

% Let's construct an A->B->C influence structure, with manually-designed message-passing and self-inhibiting kernels
[iht,ihbas,ihbasis] = makeBasis_PostSpike(ggsimcpl.ihbasprs, ggsimcpl.dt);
%toself      = ihbasis*[0 0 -1 0 0 0 0 0]' * 10;
toself      = ihbasis*[0 0 -1 -0.8 -0.6 -0.4 -0.2 0]' * 10;
toneighbour = ihbasis*[0 0 0 -0.2 1.5 -0.2 0 0]' * 10;
%plot(ggsimcpl.iht, toneighbour);

%disp(ggsimcpl);

if false
	disp('ihbas:');
	disp(size(ihbas));
	disp('toneighbour:');
	disp(size(toneighbour));
	disp('ggsimcpl.ih (originally filled in with eg gg.ihbas*gg.ih(:,1)):');
	disp(size(ggsimcpl.ih));
	disp('ggsimcpl.iht (toneighbour should match this length):');
	disp(size(ggsimcpl.iht));
end;

for which=1:3
	ggsimcpl.ih(:,which,which) = toself;
end
for which=1:2
	ggsimcpl.ih(:,which,which+1) = toneighbour;
end

% the DC. note that we want the "originator" to spike unstimulated now and again, but not the others.
ggsimcpl.dc = ggsimcpl.dc * 0 - [0, 10, 10];


% Having designed the kernels, we had better plot them to file so we understand them!
h = figure();
% orient(h, 'landscape');  % ugh, fails on my octave
clf();
for frm=1:3
	for too=1:3
		subplot(4, 3, (frm-1)*3+too);
		plot(ggsimcpl.iht, exp(ggsimcpl.iht*0), 'k--', ggsimcpl.iht, exp(ggsimcpl.ih(:, frm, too)));
		if frm==1
			title('exponentiated post-spike kernels');
		end
		if frm==3
			xlabel('time after spike (frames)');
		end
		if too==1
			ylabel('gain');
		end
		axis tight;
	end;
end;
subplot(4, 3, 11);
plot(ggsimcpl.iht, ihbasis);
axis tight;
xlabel('time after spike (frames)');
title('basis for h');
saveas(h, 'plot_simulate_coupled_network_kernels.pdf');


% for each coupling kernel, find the peak (pos and val), and write it out to file
[peakvals, peakposs] = max(ggsimcpl.ih,    [], 1);
peakvals = squeeze(peakvals);
peakposs = squeeze(peakposs) * DTsim;
csvfp = fopen('data_simulate_coupled_network_params.csv', 'w');
fprintf(csvfp, 'frm,too,peakval,peakpos\n');
for frm=1:3
	for too=1:3
		fprintf(csvfp, '%i,%i,%g,%g\n', frm, too, peakvals(frm, too), peakposs(frm, too));
	end
end
fclose(csvfp);


% OK, now we can sample from the model.
slen = 6000; %60; %TODO 600; % Stimulus length (frames) 
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


% write out a CSV file of the generated timestamps
csvfp = fopen('data_simulate_coupled_network.csv', 'w');
fprintf(csvfp, 'time,individ,dursecs\n');
for row=tspdata'
	fprintf(csvfp, '%g,%i,%g\n', row(1), row(2), 0.1);
end
fclose(csvfp);

