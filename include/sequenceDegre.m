function d = sequenceDegre(T)
	
matT = [T{:}];
d = zeros(max(max(matT)),1);
for m = 1:size(d,1)
	d(m) = size(find(matT(:)==m),1);
end
d = d';
end