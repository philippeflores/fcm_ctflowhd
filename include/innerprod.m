function y = innerprod(X,Y)
% INNERPROD : computes the inner product between two tensors
% 
% *** Input Arguments ***
%                                               ____ M ____  
%   - (X,Y) : 2 order-M tensors (double of size I x ... x I)
%
% *** Output Argument ***
% 
%   - y : innerproduct (double)
% y represents the innerproduct between order-M tensors X and Y.
% 
% *** Summary *** 
% This function returns the innerproduct between two tensors :
% y = <X,Y> y = (X|Y).
% 

if length(size(X))~=length(size(Y)) || sum(size(X)==size(Y))~=length(size(X))
    error("*** Erreur dans la fonction myInnerProduct ***\nLes deux éléments du produit ne sont pas de même taille. Fermeture du programme.\n");
end

y = X.*Y;

while length(y) > 1
    y = sum(y);
end

end