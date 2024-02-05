function val = frob(X)
% FROB : computes the frobenius norm of a tensor
%
% *** Input Arguments ***
%                                        ____ M ____
%   - X : order-M tensor (double of size I x ... x I)
%
% *** Output Argument ***
%
%   - val : frobenius norm (double)
% *val* represents the frobenius norm of the order-*M* tensor *X*.
% 
% *** Summary *** 
% This function returns the frobenius norm of a tensor: val = ||X||_F.
%
% Author: 
% name : Philippe Flores
% e-mail : flores.philipe@gmail.com
% github : github.com/philippeflores/fcm_ctflowhd

val = sqrt(innerprod(X,X));

end