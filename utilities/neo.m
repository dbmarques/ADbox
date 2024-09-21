% NEO() - Non-linear Energy Operator
%   Usage
%       [Xneo] = neo(X)
% 
%   Inputs 
%       X = data (samples,channels)
%       d = sample difference to calculate [default: 1]
% 
%   Output 
%       Xneo = NEO data
% 
% Author: Danilo Benette Marques, 2024

function [Xneo] = neo(X,d)

if nargin<2 | isempty(d)
    d = 1;
end

for i = 1:size(X,2)
    x = X(:,i);

    x = shiftdim(x);
    % xneo = x(2:end-1).^2 - x(1:end-2) .* x(3:end);
    xneo = x((d+1):end-(d)).^2 - x(1:end-(2*d)) .* x((2*d+1):end);
    xneo = [x(1:d); xneo; x(end-d+1:end)];

    Xneo(:,i) = xneo;
end
   
end