function y = permfacto(x,indPerm)
% PERMFACTO : permutes cpd factors
%
% *** Input Arguments ***
%
%   - x : cpd factors (cell of size 1xM)
% x contains the M cpd factors. x is of cell type and of size 1xM where M
% is the order of the cpd and thus the order of the tensor y (see output
% argument). Each x{m} is a IxR matrix of double containing the factors for
% the m-th variable (stored in columns),
%
%   - indPerm : permutation vector (double of size Rx1)
% indPerm represent the permutation vector. If indPerm(r) = s, y{m}(:,r)
% will be equal to x{m}(:,s).
%
% *** Output Argument ***
% 
%   - y : permuted cpd factors (cell of size 1xM)
% y contains the factor matrices of x but permuted with the permutation 
% indPerm (see Input Arguments).
%
% *** Summary ***
% This function permute cpd factors:
% y{m}(:,r) = x{m}(:,indPerm(r)).
%

M = size(x,2);
R = size(x{1},2);

y = cell(1,M);

for m = 1:M
    y{m} = zeros(size(x{m}));
    for r=1:R
        y{m}(:,r) = x{m}(:,indPerm(r));
    end
end

end