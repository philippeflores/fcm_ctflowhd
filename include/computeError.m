function [error3,error1] = computeError(y,lambda,marg,calT,indVar)
% COMPUTEERROR : compute errors of a CPD model over 1D and 3D marginals.
% 
% *** Input Arguments ***
%
%   - y : CPD factors (cell of size 1xM)
% *y* contains the *M* CPD factor matrices. *y* is a cell variable of size 
% 1xM where *M* is the number of variables in *indVar*. Each y{m} is a 
% IxR matrix of containing the factors for the m-th variable (stored in 
% columns).
%
%   - lambda : CPD loading vector (double of size Rx1)
% *lambda* contains the weights for each rank-one component of the CPD.
%
%   - marg : 3D marginals (cell of size 1xnchoosek(M,3))
% *marg* contains *T* 3D marginals as order-3 arrays of size IxIxI.
%
%   - calT : coupling (cell of size 1xT)
% *calT* contains *T* triplets of variables which represents the coupling
% of variables in *indVar*.
%
%   - indVar : vector of considered variables (double of size 1xM)
% *indVar* contains the *M* variables selected for the following analysis.
%
% *** Output Argument ***
%
%   - error3 : error on 3D marginals (double)
% *error3* is the error computed on 3D marginals. 
% error3 = \sum_{{j,k,l}\in calT} 
% ||H^{(jkl)}-cpdgen(lambda,y{j},y{k},y{l})||_F^2.
%
%   - error1 : error on 1D marginals (double)
% *error1* is the error on the *M* estimated marginals. 
% error1 = \sum_{m=1}^M ||H^{(m)}-cpdgen(lambda,y{m})||_F^2
%
% Author: 
% name : Philippe Flores
% e-mail : flores.philipe@gmail.com
% github : github.com/philippeflores/fcm_ctflowhd
%
% *** Requirements *** 
% This function requires the installation of the N-way toolbox for MATLAB.
% Reference: Anderson and Bro, "The N-way Toolbox for MATLAB", 2000

M = size(y,2);

error3 = 0;
for n = 1:size(calT,2)
	posVar = findPosVar(indVar,calT{n});
	error3 = error3 + frob(marg{n}-cpdgen({y{posVar(1)},y{posVar(2)},y{posVar(3)}},lambda))^2;
end

error1 = -inf(M,1);
nTrip = 1;
while sum(error1)<0
	posVar = findPosVar(indVar,calT{nTrip});
	for i = 1:3
		if error1(posVar(i))<0
			hist1d = squeeze(sum(sum(permute(marg{nTrip},[1:i-1 i+1:3 i]))));
			marg1d = sum((y{posVar(i)}*diag(lambda))')';
			error1(posVar(i)) = frob(hist1d-marg1d)^2;
		end
	end
	nTrip = nTrip+1;
end
error1 = sum(error1);

end