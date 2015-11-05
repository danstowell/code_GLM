function [f,df,ddf] = softplus(x)
%  [f,df,ddf] = softplus(x)

f   = log(1 + exp(x));
df  = 1  ./ (1 + exp(-x));
ddf = df .* (1 - df);

