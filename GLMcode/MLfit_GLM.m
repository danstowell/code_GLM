function [gg, fval,H] = MLfit_GLM(gg,Stim,regln,optimArgs);
%  [ggnew,fval,H] = MLfit_GLM(gg,Stim,regln,optimArgs);
% 
%  Computes the ML estimate for GLM params, using grad and hessians.
%  Assumes basis for temporal dimensions of stim filter
%
%  Inputs: 
%     gg = param struct
%     Stim = stimulus
%     regln = regularisation strength
%     optimArgs = cell array of optimization params (optional)
%
%  Outputs:
%     ggnew = new param struct (with ML params);
%     fval = negative log-likelihood at ML estimate

MAXSIZE  = 1e6 ; % TODO OCTAVE GIVES A MEM ERROR IF 1e7;    1e7;  % Maximum amount to be held in memory at once;

if nargin < 3
	regln = 0
end

% Set optimization parameters 
if nargin > 3
    opts = optimset('Gradobj','on','Hessian','on', optimArgs{:});
else
    opts = optimset('Gradobj','on','Hessian','on','display','iter');
end

% Set initial params
prs0 = extractFitPrs_GLM(gg,Stim,MAXSIZE);

lossfunc = @(x)Loss_GLM_logli(x,regln);  % a wrapper function to pass the regln parameter through on our behalf

% minimize negative log likelihood
[prs,fval] = fminunc(lossfunc,prs0,opts);
if nargout > 2 % Compute Hessian if desired
    [fv,gradval,H] = lossfunc(prs);
end

% Put returned vals back into param structure ------
gg = reinsertFitPrs_GLM(gg,prs);

%----------------------------------------------------
% % ------ Check analytic gradients, Hessians -------
% HessCheck(@lossfunc,prs0,opts);
% HessCheck_Elts(@lossfunc, [1 12],prs0,opts);
% tic; [lival,J,H]=lossfunc(prs0); toc;

