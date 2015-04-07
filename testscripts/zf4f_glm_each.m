
% script to run the GLM fitting on a sequence of zf4f 15-min data chunks

dofit = 0;   % if 0 it doesn't rerun the fit, but uses whatever we got in the last run

plotcols = {'r', 'b', 'g', 'm'};

setses = { ...
struct('dataname', 'session2a', 'variantname', '1', 'indexmapper', [4,3,2,1], 'startsecs',  300, 'endsecs', 1200), ...
struct('dataname', 'session2a', 'variantname', '2', 'indexmapper', [4,3,2,1], 'startsecs', 1200, 'endsecs', 2100), ...
struct('dataname', 'session2b', 'variantname', '1', 'indexmapper', [4,3,2,1], 'startsecs', 2100, 'endsecs', 3000), ...
struct('dataname', 'session2b', 'variantname', '2', 'indexmapper', [4,3,2,1], 'startsecs', 3000, 'endsecs', 3900), ...
struct('dataname', 'session3a', 'variantname', '1', 'indexmapper', [1,2,3,4], 'startsecs',  350, 'endsecs', 1250), ...
struct('dataname', 'session3a', 'variantname', '2', 'indexmapper', [1,2,3,4], 'startsecs', 1250, 'endsecs', 2150), ...
struct('dataname', 'session3b', 'variantname', '1', 'indexmapper', [1,2,3,4], 'startsecs', 2150, 'endsecs', 3050), ...
struct('dataname', 'session3b', 'variantname', '2', 'indexmapper', [1,2,3,4], 'startsecs', 3050, 'endsecs', 3950), ...
};

if dofit
	resultspos = struct();
	resultsval = struct();
	for whichset=1:size(setses,2)
		d = setses{whichset};
		runname = sprintf('%s%s', d.dataname, d.variantname);
		[resultspos.(runname), resultsval.(runname)] = testscript_GLM_zf4f(d.dataname, d.variantname, d.indexmapper, d.startsecs, d.endsecs);
	end
end


seqlist = {'session2a1', 'session2a2', 'session2b1', 'session2b2'};

h = figure(4);
clf();
for whichn=1:4
	% plot peak level and peak time
	plotdata_p = zeros(4);
	plotdata_t = zeros(4);
	for fromn=1:4
		for whichsess=1:4
			plotdata_pos(whichsess, fromn) = resultspos.(seqlist{whichsess})(fromn, whichn);
			plotdata_val(whichsess, fromn) = resultsval.(seqlist{whichsess})(fromn, whichn);
		end;
	end;
	subplot(4, 2, whichn*2-1);
	hold on;
	for fromn=1:4
		if whichn==fromn
			plotcol = 'kx--';
		else
			plotcol = sprintf('%sx-', plotcols{fromn});
		end
		plot(plotdata_val(:, fromn), plotcol);
	end
	hold off;
	ylabel(sprintf('Bird %i', whichn));
	%ylim([0, 200]);
	xlim([0.5, 4.5]);
	set(gca,'XTick',[]);
	legend('from 1', 'from 2', 'from 3', 'from 4', 'location', 'northeast');
	if whichn==1
		title('Value at kernel peak');
	end;
	subplot(4, 2, whichn*2);
	hold on;
	for fromn=1:4
		if whichn==fromn
			plotcol = 'kx--';
		else
			plotcol = sprintf('%sx-', plotcols{fromn});
		end
		plot(plotdata_pos(:, fromn), plotcol);
	end
	hold off;
	%ylim([0, 10]);
	xlim([0.5, 4.5]);
	set(gca,'XTick',[]);
	legend('from 1', 'from 2', 'from 3', 'from 4', 'location', 'northeast');
	if whichn==1
		title('Time of kernel peak (sec)');
	end;
end;
saveas(h, sprintf('zf4f_glm_evolution_%s.png', 'session2_15'))

