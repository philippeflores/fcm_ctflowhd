function strLabelOut = changeLabel(strLabelIn)
% This function permits to change the labels of a flow cytometry 
% experiment.
%
% Author: 
% name : Philippe Flores
% e-mail : flores.philipe@gmail.com
% github : github.com/philippeflores/fcm_ctflowhd

m = 1;
M = size(strLabelIn,2);
strLabelOut = cell(size(strLabelIn));
while m<=M
    strInput = sprintf("Variable #%d\nOld Label %s\nNew Label : ",m,strLabelIn{m});
    strLabelOut{m} = input(strInput,'s');
    m = m+1;
end

end