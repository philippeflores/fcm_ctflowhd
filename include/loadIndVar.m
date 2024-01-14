function [indVarOut,M,strLabel] = loadIndVar(X,indVarIn,fcsHdr)
% LOADINDVAR : load variables in a flow cytometry fcs file
% 
% *** Input Arguments ***
%
%   - X : flow cytometry data (double of size NxM)
% X is the observation matrix contained in a fcs file. 
%
%   - indVarIn : vector of considered variables (double of size Mx1)
% indVarIn contains the M variables required for the following analysis. 
% All those variables which are available for further analysis will be 
% stored in indVarOut (see Output Arguments).
%
%   - fcsHdr : fcs file properties (1x1 struct)
% fcsHdr contains all properties stored inside the header of a fcsfile.
% For example, compensation matrix marker labels, number of events, etc.
%
% *** Output Argument ***
% 
%   - indVarOut : vector of considered variables (double of size Mx1)
% indVarOut contains the M variables actually selected for the following 
% analysis. If indVarIn is empty, all available variables are selected.
%
%   - M : number of considered variables (double)
% M represents the number of variables selected at the end of this
% function. M is the length of indVarOut and also the length of indVarIn if
% all variables in indVarIn are available for further analysis.
%
%   - strLabel : labels for flow cytometry variables (cell of size 1xM)
% strLabel store the labels for each variable of indVarOut. If a
% name is present in the fcsHdr, this name will be the name stored in
% strLabel.
% 
% *** Instructions *** 
% If you want to select all possibly analyzable variables, assign
%  >>> indVarIn = [];
% If so, all variables will be selected for further analysis.
% 
% If you want to analyze only a subset of variables inside a fcs file,
% assign to indVarIn the indices of those selected variables. To find those
% indices, see fcsHdr.par (see Input Arguments).
%

if isempty(indVarIn) || any(indVarIn>size(X,2)-1)
	strFluoTemp = {fcsHdr.par.name};
	indVarOut = zeros(1,size(X,2));
	for m = 1:size(X,2)
		nameTemp = strFluoTemp{m};
        if strcmpi(nameTemp(end-1:end),'-A')
			indVarOut(m) = 1;
            if strcmpi(nameTemp(1:3),'FSC')
                indVarOut(m) = 0;
            end
            if strcmpi(nameTemp(1:3),'SSC')
                indVarOut(m) = 0;
            end
        end
	end
	temp = 1:size(X,2);
	indVarOut = temp(logical(indVarOut));
else
	indVarOut = indVarIn;
end

M = length(indVarOut);

strLabel = {fcsHdr.par.marker};
strLabel = strLabel(indVarOut);

for m = 1:size(strLabel,2)
    temp = strLabel{m};
    if strcmp(temp(end-1:end),'-A')
		strLabel{m} = temp(1:end-2);
    end
end

end