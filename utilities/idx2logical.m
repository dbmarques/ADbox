% IDX2LOGICAL() - Converts indices to logical
%   Usage
%       [idx01] = idx2logical(idx,N)
% 
%   Inputs 
%       idx = indices (vector)
%       N   = Length of logical vector [default: last index]
%   Output 
%       idx01 = logical vector
% 
% Author: Danilo Benette Marques, 2024

function [idx01] = idx2logical(idx,N)

if nargin<2 | isempty(N)
    N = length(idx);
end

idx01 = zeros(N,1);
idx01(idx) = 1;

idx01 = logical(idx01);

end