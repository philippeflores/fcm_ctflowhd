function calTfull = fullCoupled(M)
% This function outputs the coupling that contains all possible triplets of
% variables in [1:M]. 
% 
% Author: 
% name : Philippe Flores
% e-mail : flores.philipe@gmail.com
% github : github.com/philippeflores/fcm_ctflowhd

calTfull = cell(1,nchoosek(M,3));
countT = 1;
for j = 1:(M-2)
	for k = (j+1):(M-1)
		for l = (k+1):M
			calTfull{countT} = [j k l];
			countT = countT+1;
		end
	end
end

end