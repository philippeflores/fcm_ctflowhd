function V = TboolMat(calT)
% This function outputs the incidence matrix *V* of the hypergraph *calT*.
%
% Author: 
% name : Philippe Flores
% e-mail : flores.philipe@gmail.com
% github : github.com/philippeflores/fcm_ctflowhd

calTmat = reshape([calT{:}],3,[])';
M = max(calTmat(:));

V = zeros(size(calT,2),M);
for t = 1:size(calT,2)
    V(t,calTmat(t,:)) = 1;
end

end