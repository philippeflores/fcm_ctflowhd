function y = permfacto(x,perm)
% PERMFACTO : permutes cpd factors
%
% Author: 
% name : Philippe Flores
% e-mail : flores.philipe@gmail.com
% github : github.com/philippeflores/fcm_ctflowhd
%
% *** Input Arguments ***
%
%   - x : cpd factors (cell of size 1xM)
% x contains the M cpd factors. x is of cell type and of size 1xM where M
% is the order of the CPD and thus the order of the tensor y (see output
% argument). Each x{m} is a IxR matrix of double containing the factors for
% the m-th variable (stored in columns),
%
%   - perm : permutation vector (double of size Rx1)
% perm represent the permutation vector. If perm(r) = s, y{m}(:,r) will be 
% equal to x{m}(:,s).
%
% *** Output Argument ***
% 
%   - y : permuted CPD factors (cell of size 1xM)
% y contains the factor matrices of x but permuted with the permutation 
% perm (see Input Arguments).
%
% *** Summary ***
% This function permute CPD factors:
% y{m}(:,r) = x{m}(:,perm(r)).

M = size(x,2);
R = size(x{1},2);

y = cell(1,M);

for m = 1:M
    y{m} = zeros(size(x{m}));
    for r=1:R
        y{m}(:,r) = x{m}(:,perm(r));
    end
end

end