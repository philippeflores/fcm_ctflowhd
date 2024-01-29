function d = sequenceDegre(calT)
% This function outputs the sequence of degrees for a coupling calT.
%
% Author: 
% name : Philippe Flores
% e-mail : flores.philipe@gmail.com
% github : github.com/philippeflores/fcm_ctflowhd

matT = [calT{:}];
d = zeros(max(max(matT)),1);
for m = 1:size(d,1)
	d(m) = size(find(matT(:)==m),1);
end
d = d';
end