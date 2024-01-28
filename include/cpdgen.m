function y = cpdgen(x,lambda)
% CPDGEN : computes the tensor from factors and weights
% 
% *** Input Arguments ***
%
%   - x : cpd factors (cell of size 1xM)
% x contains the M cpd factor matrices. x is of cell type and of size 1xM 
% where M is the order of the cpd and thus the order of the tensor y (see 
% output argument). Each x{m} is a IxR matrix of double containing the 
% factors for the m-th variable (stored in columns),
%
%   -not mandatory- lambda : cpd weights (double of size Rx1)
% lambda contains the weights for each rank one component. This function 
% support the non presence of lambda and if so considers that weights are 
% all equal to 1.
% 
% *** Output Argument ***
%                                             ____ M ____  
%   - y : recontructed tensor (double of size I x ... x I)
% y represents the order-M reconstructed tensor from the factors contained
% with weight lambda (see input arguments).
% 
% *** Summary *** 
% This function returns the recontrusted tensor from cpd factors:
% y_{i_1...i_M} = \sum_{r=1}^R lambda(r) \prod_{m=1}^M x{m}(i_m,r).
% 
% *** Requirements *** 
% This function requires the installation of the N-way toolbox for MATLAB.
% Reference: Anderson and Bro, "The N-way Toolbox for MATLAB", 2000
%

A = x{1};
R = size(A,2);
if exist('lambda','var')==0
    lambda = ones(R,1);
end

for j = 1:R
    A(:,j) = lambda(j)*A(:,j);
end

x{1} = A;

y = nmodel(x);

end