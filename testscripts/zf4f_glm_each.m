
% script to run the GLM fitting on a sequence of zf4f 15-min data chunks

dofit = 1;   % if 0 it doesn't rerun the fit, but uses whatever we got in the last run

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
struct('dataname', 'session2full', 'variantname', '', 'indexmapper', [4,3,2,1], 'startsecs',  300, 'endsecs', 3900), ...
struct('dataname', 'session3full', 'variantname', '', 'indexmapper', [1,2,3,4], 'startsecs',  350, 'endsecs', 3950), ...
};

if dofit
	numcalls   = struct();
	resultspos = struct();
	resultsval = struct();
	negloglis  = struct();
	for whichset=1:size(setses,2)
		d = setses{whichset};
		runname = sprintf('%s%s', d.dataname, d.variantname);
		[numcalls.(runname), resultspos.(runname), resultsval.(runname), negloglis.(runname)] = testscript_GLM_zf4f(d.dataname, d.variantname, d.indexmapper, d.startsecs, d.endsecs);
	end
end


seqlists = { ...
{'session2_15', {'session2a1', 'session2a2', 'session2b1', 'session2b2'}}, ...
{'session3_15', {'session3a1', 'session3a2', 'session3b1', 'session3b2'}}, ...
{'session2full60', {'session2full'}}, ...
{'session3full60', {'session3full'}}, ...
};

for whichseq=1:size(seqlists, 2)
	seqlist = seqlists{whichseq};
	oursetses = seqlist{2};
	numsesses = length(oursetses);

	outfnamestem = sprintf('zf4f_glm_evolution_%s', seqlist{1});

	csvfname = sprintf('outcsv/%s_1d.csv', outfnamestem);
	disp(csvfname);
	csvfp_1d = fopen(csvfname, 'w');
	fprintf(csvfp_1d, 'runname,individ,numcalls\n');

	csvfname = sprintf('outcsv/%s_2d.csv', outfnamestem);
	disp(csvfname);
	csvfp_2d = fopen(csvfname, 'w');
	fprintf(csvfp_2d, 'runnname,frm,too,peakval,peaklag\n');

	h = figure(4);
	clf();
	for whichn=1:4
		% plot peak level and peak time
		plotdata_num = zeros(numsesses,1);
		plotdata_pos = zeros(numsesses,4);
		plotdata_val = zeros(numsesses,4);
		for fromn=1:4
			for whichsess=1:numsesses
				plotdata_num(whichsess)        =     numcalls.(oursetses{whichsess})(whichn);
				plotdata_pos(whichsess, fromn) = 1 / max(1e-1, resultspos.(oursetses{whichsess})(fromn, whichn));
				plotdata_val(whichsess, fromn) =     resultsval.(oursetses{whichsess})(fromn, whichn);
			end;
		end;
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		subplot(4, 3, whichn*3-2);
		plotcol = sprintf('%sx-', plotcols{whichn});
		plot(plotdata_num, plotcol);
		xlim([0.5, 4.5]);
		ylim([0, max(plotdata_num)+10]);
		set(gca,'XTick',[]);
		ylabel(sprintf('Bird %i', whichn));
		if whichn==1
			title('Num calls emitted in each 15min period');
		end;
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		subplot(4, 3, whichn*3-1);
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
			title('Relative intensity at kernel peak');
		end;
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		subplot(4, 3, whichn*3);
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
		ylabel(sprintf('Bird %i', whichn));
		%ylim([0, 10]);
		xlim([0.5, 4.5]);
		set(gca,'XTick',[]);
		legend('from 1', 'from 2', 'from 3', 'from 4', 'location', 'northeast');
		if whichn==1
			title('1/latency of kernel peak (1/sec)');
		end;

		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% write csv data
		for whichsess=1:numsesses
			fprintf(csvfp_1d, '%s,%i,%i\n', oursetses{whichsess}, whichn, plotdata_num(whichsess));
			for fromn=1:4
				fprintf(csvfp_2d, '%s,%i,%i,%g,%g\n', oursetses{whichsess},fromn, whichn, plotdata_val(whichsess, fromn), 1/plotdata_pos(whichsess, fromn));
			end;
		end;

	end;
	saveas(h, sprintf('outplot/%s.png', outfnamestem));
	fflush(csvfp_1d);
	fflush(csvfp_2d);
	fclose(csvfp_1d);
	fclose(csvfp_2d);
end;

