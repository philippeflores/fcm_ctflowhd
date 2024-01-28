function [A,B] = computeError(y,lambda,marg,calT,indVar)
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
%   - histX : 3D marginals (cell of size 1xNtriplets)
% histX contains all 3D marginals. Even if the method was applied to a
% subset of triplet. The error on 3D marginals is computed on every
% possible 3D marginals.
%
%   - triplets : coupling strategy (cell of size 1xNtriplets)
% triplets contains Ntriplets triplets of variables which represents the coupling
% strategy of variables [1:M].
%
%   - indVar : vector of considered variables (double of size Mx1)
% indVar contains the M variables selected for the following analysis.
%
% *** Output Argument ***
%
%   - A : error on 3D marginals (double)
% A is the error on marginals computed on all possible 3D marginals. A =
% \sum_{{j,k,l}\in Tall} ||H^{(jkl)}-cpdgen(lambda,y{j},y{k},y{l})||_F^2.
%
%   - B : error on 1D marginals (double)
% B is the error on the M estimated marginals. B = \sum_{m=1}^M ||H^{(m)}
% -cpdgen(lambda,y{m})||_F^2
%

M = size(y,2);

A = 0;
for n = 1:size(calT,2)
	posVar = findPosVar(indVar,calT{n});
	A = A + frob(marg{n}-cpdgen({y{posVar(1)},y{posVar(2)},y{posVar(3)}},lambda))^2;
end

B = -inf(M,1);
nTrip = 1;
while sum(B)<0
	posVar = findPosVar(indVar,calT{nTrip});
	for i = 1:3
		if B(posVar(i))<0
			hist1d = squeeze(sum(sum(permute(marg{nTrip},[1:i-1 i+1:3 i]))));
			marg1d = sum((y{posVar(i)}*diag(lambda))')';
			B(posVar(i)) = frob(hist1d-marg1d)^2;
		end
	end
	nTrip = nTrip+1;
end
B = sum(B);

end