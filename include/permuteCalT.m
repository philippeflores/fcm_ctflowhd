function calTout = permuteCalT(calT)
% This function permits to perform a permutation of variables of a
% coupling. Let us define a triplet calT :
%   calT = {{j,k,l}\subset [1:M]} such that Card(calt) = T.
% This coupling admits a sequence of degrees d. Let us define the
% permutation perm such that perm(d) is in ascending order.
% Finally, calTout is defined as follows:
%   calTout = {{perm(j),perm(k),perm(l) | {j,k,l} \in calT}
% Therefore, calTout has a sequence of degrees equal to perm(d) which has
% increasing values.
%
% Author: 
% name : Philippe Flores
% e-mail : flores.philipe@gmail.com
% github : github.com/philippeflores/fcm_ctflowhd

calTout = calT; for t = 1:size(calT,2), calTout{t} = pi*calT{t}; end

d = sequenceDegre(calT);
[~,indSort] = sort(d,'ascend');

for i = 1:size(indSort,2)
    m = indSort(i);
    for t = 1:size(calT,2)
        calTout{t}(calTout{t}==m*pi) = i;
    end
end

for t = 1:size(calT,2)
    calTout{t} = sort(calTout{t},'ascend');
end
matCalT = reshape([calTout{:}],3,[])';
[~,temp] = sortrows(matCalT,'ascend');
calTout = calTout(temp);

end