function [error3,error1] = computeError(y,lambda,marg,calT,indVar)
% COMPUTEERROR : compute errors on 1D and 3D marginals of a cpd
% 
% *** Input Arguments ***
%
%   - y : cpd factors (cell of size 1xM)
% y contains the M cpd factor matrices. y is of cell type and of size 1xM 
% where M is the order of the cpd. Each y{m} is a IxR matrix of double 
% containing the factors for the m-th variable (stored in columns).
%
%   - lambda : cpd weights (double of size Rx1)
% lambda contains the weights for each rank-one components of the cpd.
%
%   - marg : 3D marginals (cell of size 1xT)
% marg contains all 3D marginals. Even if the method was applied to a
% subset of triplet (which corresponds to T triplets). The error on 3D 
% marginals is computed on every possible 3D marginals.
%
%   - calT : coupling (cell of size 1xT)
% calT contains T triplets of variables which represents the coupling
% of variables in indVar.
%
%   - indVar : vector of considered variables (double of size 1xM)
% indVar contains the M variables selected for the following analysis.
%
% *** Output Argument ***
%
%   - error3 : error on 3D marginals (double)
% error3 is the error on marginals computed on all possible 3D marginals. 
% A = \sum_{{j,k,l}\in calTfull} 
% ||H^{(jkl)}-cpdgen(lambda,y{j},y{k},y{l})||_F^2.
%
%   - error1 : error on 1D marginals (double)
% error1 is the error on the M estimated marginals. B = \sum_{m=1}^M 
% ||H^{(m)}-cpdgen(lambda,y{m})||_F^2
%
% Author: 
% name : Philippe Flores
% e-mail : flores.philipe@gmail.com
% github : github.com/philippeflores/fcm_ctflowhd

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