function [numcalls, peakpos, peakval] = dofit_fromcsv_GLM_zf4f(csvpath, runlabel, indexmapper, startsecs, endsecs, plotpath, csvpath)
% [peakpos, peakval] = dofit_fromcsv_GLM_zf4f(csvpath, runlabel, indexmapper, startsecs, endsecs, plotpath, csvpath)
%
% load some zf4f-format data and analyse "as if" it were cell spiking data. returns analysed data.
% also does a plot and writes it to a file in the folder named by 'plotpath'. if 'plotpath' is '' or 0 it DOESN'T plot. to plot in cwd use '.'
% 'csvpath' parameter is analogous, and is about writing CSV data out to file

% TODO add 'csvpath' parameter, treated analogously to plotpath

% TODO adapt the code so that it correctly adapts to varying k

global RefreshRate;
RefreshRate = 1;  % the "RefreshRate" is the samplerate of the stimulus (in Hz). I don't currently use stimulus so I set it to 1. Below, "DTsim" sets the time-resultion used in the model.

plotcols = {'r', 'b', 'g', 'm'};

numcalls = zeros(4,1);
peakpos = zeros(4);
peakval = zeros(4);

printf('dofit_fromcsv_GLM_zf4f(%s, %s, %s, %i, %i, %s)\n', csvpath, runlabel, mat2str(indexmapper), startsecs, endsecs, plotpath);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load the CSV data
events = csvread(csvpath);
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
	% NOTE: the call above imposes an "absref" (absolute refactory period) of 10 * DTsim. Happily this is extremely reasonable in the zf4f case.
	gg0.ih = gg0.ih*0;  % Initialize to zero
	gg0.dc = gg0.dc*0;  % Initialize to zero

	gg0.tsp = tsp{whichn};   % cell spike times (vector)
	gg0.tsp2 = tsp(setdiff(1:4, whichn));  % spike trains from "coupled" cells (cell array of vectors)
	gg0.tspi = 1; % 1st spike to use for computing likelihood (eg, can ignore 1st n spikes)
	[gg{whichn}, negloglivalx] = MLfit_GLM(gg0, Stim, opts); % do ML (requires optimization toolbox)
	printf('MLfit #%i gets neglogli %g\n', whichn, negloglivalx);
end


%% --- Plot results ----------------------------
if plotpath
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
		title(sprintf('Bird %i: exp(kernels) %s', whichn, runlabel));
		legend('from 1', 'from 2', 'from 3', 'from 4', 'location', 'northeast');
		axis tight;
	end;
	xlabel('time (frames)')

	disp(sprintf('Saving %s/zf4f_glm_kernels_%s.png', plotpath, runlabel));
	saveas(h, sprintf('%s/zf4f_glm_kernels_%s.png', plotpath, runlabel));

	sleep(2);
else
	disp '  (not plotting)';
end



if csvpath
	% TODO - refactor the "plotdata" stuff out of zf4f_glm_each and put it in a separate statsgetter func file
	%  2d stats: peak posses, peak mags
	%  1d stats: num calls
	%  0d stats: likelihood under fitted model
	disp('           TODO');
	sleep(2);
else
	disp '  (not writing csv)';
end

