function [calT,d] = balancedCoupling(M,T)

calT = cell(1,T);

matT = zeros(T,M);
nMLeq = floor(T/M);
[ML,nPerm] = genLyndon3(M,nMLeq);
nML = size(ML,1);
indDispo = 1:nML;


if mod(M,3)==0
    if nMLeq==(nML-1)
        indMLeq = 1:nMLeq;
    else
        indMLeq = randperm(nML-2,nMLeq)+1;
    end
else
    indMLeq = randperm(nML-1,nMLeq)+1;
end
indDispo(indMLeq) = [];

countT = 1;
for i = 1:size(indMLeq,2)
    for m = 0:(nPerm(indMLeq(i))-1)
        matT(countT,:) = circshift(ML(indMLeq(i),:),m);
        countT = countT+1;
    end
end

Treste = T-countT+1;
if mod(M,3)==0 && Treste>=(M/3)
    i = indDispo(end);
    indDispo(end) = [];
    for m = 0:(nPerm(i)-1)
        matT(countT,:) = circshift(ML(i,:),m);
        countT = countT+1;
    end
    Treste = T-countT+1;
end

mot = ML(1,:);
if mod(M,3)==0
    for i = 1:Treste
        alpha = floor(i/3);
        matT(countT,:) = circshift(mot,(i-1)*3+alpha);
        countT = countT+1;
    end
else
    for i = 1:Treste
        matT(countT,:) = circshift(mot,(i-1)*3);
        countT = countT+1;
    end
end

matT = sortrows(matT,'descend');

for t = 1:T
    calT{t} = (1:M).*(matT(t,:));
    calT{t}(calT{t}==0) = [];
end

calT = permuteCalT(calT);
d = sequenceDegre(calT);

end