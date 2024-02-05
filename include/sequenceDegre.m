function d = sequenceDegre(calT)
% This function outputs the sequence of degrees for a coupling *calT*.
%
% Author: 
% name : Philippe Flores
% e-mail : flores.philipe@gmail.com
% github : github.com/philippeflores/fcm_ctflowhd

calTvec = [calT{:}];
d = zeros(max(max(calTvec)),1);
for m = 1:size(d,1)
	d(m) = size(find(calTvec(:)==m),1);
end
d = d';
end