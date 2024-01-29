function V = TboolMat(calT)
% This function outputs the incidence matrix V of the hypergraph calT.
%
% Author: 
% name : Philippe Flores
% e-mail : flores.philipe@gmail.com
% github : github.com/philippeflores/fcm_ctflowhd

matCalT = reshape([calT{:}],3,[])';
M = max(matCalT(:));

V = zeros(size(calT,2),M);
for t = 1:size(calT,2)
    V(t,matCalT(t,:)) = 1;
end

end