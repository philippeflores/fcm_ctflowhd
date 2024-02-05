function indTrip = findPosVar(indVar,jkl)
% FINDPOSVAR : find the variable position for one triplet
% 
% *** Input Arguments ***
%
%   - indVar : vector of considered variables (double of size 1xM)
% *indVar* contains the *M* variables considered in the analysis. For a fcs 
% file (flow cytometry), *indVar* contains the variables for further 
% analysis.
%
%   - jkl : triplet (double of size 1x3)
% *jkl* is a vector that contains the 3 indices [j k l] where {j,k,l} are 3
% variables taken from *indVar*.
%
% *** Output Argument ***
%
%   - *indTrip* : position triplet (double of size 1x3)
% *indTrip* contains the position in *indVar* of the jkl vector.
%
% *** Example ***
%
%   >>> indVar = [6 8 10 12];
%   >>> jkl = [6 8 12];
%   >>> findPosVar(indVar,jkl)
%           will return indTrip = [1 2 4].
%
% Author: 
% name : Philippe Flores
% e-mail : flores.philipe@gmail.com
% github : github.com/philippeflores/fcm_ctflowhd




temp = indVar==jkl';
indTrip = [find(temp(1,:)) find(temp(2,:)) find(temp(3,:))];

end