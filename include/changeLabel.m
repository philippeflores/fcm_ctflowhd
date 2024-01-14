function strLabelOut = changeLabel(strLabelIn,varargin)

strChange = 'man';

if (rem(length(varargin),2)==1)
    error('Optional parameters should always go by pairs');
else
    for i=1:2:(length(varargin)-1)
        switch upper(varargin{i})
            case 'STRCHANGE'
                strChange = varargin{i+1};
            case 'STRLABELINPUT'
                strLabelInput = varargin{i+1};
            otherwise
                error(['Unrecognized option: ''' varargin{i} '''']);
        end
    end
end

m = 1;
M = size(strLabelIn,2);
strLabelOut = cell(size(strLabelIn));
if strcmpi(strChange,'MAN')
	while m<=M
		strInput = sprintf("Variable #%d\nOld Label %s\nNew Label : ",m,strLabelIn{m});
		strLabelOut{m} = input(strInput,'s');
		m = m+1;
	end
elseif strcmpi(strChange,'LOAD')
	if exist('strLabelInput','var') && size(strLabelInput,2)==size(strLabelIn,2)
		strLabelOut = strLabelInput;
	else
		fprintf("Wrong input parameter 'strLabelInput'. Labelling changes are cancelled.\n")
		strLabelOut = strLabelIn;
	end
else
	fprintf("Wrong input parameter 'strChange'. 'MAN' is applied by default.\n")
	while m<=M
		strInput = sprintf("Variable #%d\nOld Label %s\nNew Label : ",m,strLabelIn{m});
		strLabelOut{m} = input(strInput,'s');
		m = m+1;
	end
end
	

end