function val = innerprod(X,Y)
% INNERPROD : computes the inner product between two tensors
%
% *** Input Arguments ***
%                                               ____ M ____  
%   - (X,Y) : 2 order-M tensors (double of size I x ... x I)
%
% *** Output Argument ***
% 
%   - val : innerproduct (double)
% val represents the innerproduct between order-M tensors X and Y.
% 
% *** Summary *** 
% This function returns the innerproduct between two tensors :
% val = <X,Y> or val = (X|Y).
%
% Author: 
% name : Philippe Flores
% e-mail : flores.philipe@gmail.com
% github : github.com/philippeflores/fcm_ctflowhd


if length(size(X))~=length(size(Y)) || sum(size(X)==size(Y))~=length(size(X))
    error("*** Error in the function innerprod ***\nThe two elements of the inner product must have the same size.\n");
end

val = X.*Y;

while length(val) > 1
    val = sum(val);
end

end