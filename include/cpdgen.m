function Y = cpdgen(y,lambda)
% CPDGEN : compute a tensor from its factor matrices and its loading
% vector.
% 
% *** Input Arguments ***
%
%   - y : CPD factors (cell of size 1xM)
% *y* contains the *M* CPD factor matrices. *y* is of cell type and of size 
% 1xM where *M* is the order of the CPD and thus the order of the tensor 
% *Y* (see output argument). Each y{m} is a IxR matrix of double containing 
% the factors for the m-th variable (stored in columns),
%
%   -not mandatory- lambda : CPD weights (double of size Rx1)
% *lambda* contains the weights for each rank one component. This function 
% supports the non presence of *lambda* and if so considers that weights 
% are all equal to 1/R.
% 
% *** Output Argument ***
%                                             ____ M ____  
%   - Y : recontructed tensor (double of size I x ... x I)
% *Y* represents the order-*M* reconstructed tensor from the factors 
% contained in *y* with weight *lambda* (see input arguments).
% 
% *** Summary *** 
% This function returns the recontrusted tensor from CPD factors:
% Y_{i_1...i_M} = \sum_{r=1}^R lambda(r) \prod_{m=1}^M y{m}(i_m,r).
% 
% Author: 
% name : Philippe Flores
% e-mail : flores.philipe@gmail.com
% github : github.com/philippeflores/fcm_ctflowhd
%
% *** Requirements *** 
% This function requires the installation of the N-way toolbox for MATLAB.
% Reference: Anderson and Bro, "The N-way Toolbox for MATLAB", 2000

A1 = y{1};
R = size(A1,2);
if exist('lambda','var')==0
    lambda = ones(R,1)/R;
end

for j = 1:R
    A1(:,j) = lambda(j)*A1(:,j);
end

y{1} = A1;

Y = nmodel(y);

end