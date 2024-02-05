function [flag,varLink] = testStratCoupling(calT,M)
% This test function permits to check if the hypergraph calT is connex.
% This means that calT contains all variables and all pairs of variables 
% are connected. More information on hypergraph connectivity can be found
% in the thesis manuscript of the author.
%
% Author:
% name : Philippe Flores
% e-mail : flores.philipe@gmail.com
% github : github.com/philippeflores/fcm_ctflowhd

T = size(calT,2);
matTriplets = reshape([calT{:}],3,[])';
varLink = {};

if size(unique(matTriplets(:)),1)==M
    flag = 1;

    m = 1;
    varLink = cell(1,M);
    sizeLink = 0;
    flag2 = 0;
    while flag2==0
        indTemp = [];
        for nTrip = 1:T
            if any(matTriplets(nTrip,:)==m)
                indTemp = [indTemp nTrip];
            end
        end
        varLink{sizeLink+1} = unique(reshape(matTriplets(indTemp,:),1,[]));
        sizeLink = sizeLink+1;
        newSizeLink = 0;
        while newSizeLink~=sizeLink
            sizeLink = newSizeLink;
            for j = sizeLink:-1:2
                if size(varLink{1},2)+size(varLink{j},2)~=size(unique([varLink{1} varLink{j}]))
                    varLink{1} = unique([varLink{1} varLink{j}]);
                    varLink{j} = [];
                    newSizeLink = newSizeLink-1;
                end
            end
            newSizeLink = 0;
            for j = 1:M
                if ~isempty(varLink{j})
                    newSizeLink = newSizeLink+1;
                end
            end
        end
        if m>M-1
            flag2 = 1;
        else
            m = m+1;
        end
    end

    if ~isempty(varLink{2})
        flag = 2;
    end

elseif size(unique(matTriplets(:)),1)<M
    flag = -M+size(unique(matTriplets(:)),1);
elseif size(unique(matTriplets(:)),1)>M
    flag = M+size(unique(matTriplets(:)),1);
end

end