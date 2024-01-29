function calTout = permuteCalT(calT)

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