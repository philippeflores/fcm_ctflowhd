function y = frob(X)
% FROB : computes the frobenius norm of a tensor
% 
% *** Input Arguments ***
%                                        ____ M ____
%   - X : order-M tensor (double of size I x ... x I)
%
% *** Output Argument ***
%
%   - y : frobenius norm (double)
% y represents the frobenius norm of the order-M tensor.
% 
% *** Summary *** 
% This function returns the frobenius norm of a tensor:
% y = ||X||_F.
%

y = sqrt(innerprod(X,X));

end