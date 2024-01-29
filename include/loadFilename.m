function [X,fcsHdr,strFile] = loadFilename(filename)

listAvailable = {dir("data/").name};
listAvailable = listAvailable(3:end);
for i = 1:size(listAvailable,2), listAvailable{i} = listAvailable{i}(1:end-4); end

if isempty(listAvailable)
    error("There are no valid '.fcs' to be selected.\n")
end

flag = 0;

filename = char(filename);

if any(ismember([listAvailable],filename))
	strFile = filename;
	flag = 1;
elseif flag==0 && length(filename)>7 && strcmpi(filename(1:7),'./data/')
	strFile = filename(8:end);
	if length(strFile)>4 && strcmpi(filename(end-3:end),'.fcs')
		strFile = strFile(1:end-4);
	end
	if any(ismember([listAvailable],strFile)), flag = 1; end
elseif flag==0 && length(filename)>5 && strcmpi(filename(1:5),'data/')
	strFile = filename(6:end);
	if length(strFile)>4 && strcmpi(filename(end-3:end),'.fcs')
		strFile = strFile(1:end-4);
	end
	if any(ismember([listAvailable],strFile)), flag = 1; end
elseif flag==0 && length(filename)>4 && strcmpi(filename(end-3:end),'.fcs')
	strFile = filename(1:end-4);
	if any(ismember([listAvailable],strFile)), flag = 1; end
end

while flag == 0
	strInput = strcat("File(s) available: \n", sprintf("%s\t",listAvailable{:}), "\nChoose a file to analyze: ");
	strFile = input(strInput,'s');
	if any(ismember([listAvailable],strFile))
		flag = -1;
	else
		fprintf("Wrong entry. Try again.\n");
	end
end

filename = strcat('data/',strFile,'.fcs');

[X,fcsHdr] = fca_readfcs(filename);

end