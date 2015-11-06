function [f,df,ddf] = softplus(x)
%  [f,df,ddf] = softplus(x)
%
% The softplus nonlinearity, a smooth approximation of recitification.
% NOTE: Internally we have added a rescaling factor to make the elbow sharper. This highlights that the nonlinearity is scale-dependent.

k = 10;  % rescaling factor -- higher means the elbow is sharper, bringing it closer to a "true" rectification nonlinearity

f   = log(1 + exp(k * x))/k;
df  = 1  ./ (1 + exp(-k * x));
ddf = k * df .* (1 - df);

