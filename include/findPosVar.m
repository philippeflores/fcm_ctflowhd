function indTrip = findPosVar(indVar,T)
% FINDPOSVAR : find the variable position for one triplet
% 
% *** Input Arguments ***
%
%   - indVar : vector of considered variables (double of size Mx1)
% indVar contains the M variables considered in the analysis. For a fcs 
% file (flow cytometry), indVar selects the variables for further analysis.
%
%   - T : triplet (double of size 1x3)
% T is a vector that contains the 3 indices [j k l] where (j,k,l) are 3
% variables taken from indVar.
%
% *** Output Argument ***
%
%   - v : position triplet (double of size 1x3)
% v contains the position in indVar of the T vector.
%
% *** Example ***
%
%   - indVar = [6 8 10 12];
%   - T = [6 8 12];
%   - FINDPOSVAR will return v = [1 2 4].
% 

temp = indVar==T';
indTrip = [find(temp(1,:)) find(temp(2,:)) find(temp(3,:))];

end