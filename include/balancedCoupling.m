function [calT,d] = balancedCoupling(M,T)
% This function is generating a balanced coupling. It corresponds to the
% algorithm 4 of the thesis manuscript of the author.
%
% Author: 
% name : Philippe Flores
% e-mail : flores.philipe@gmail.com
% github : github.com/philippeflores/fcm_ctflowhd

calT = cell(1,T);

V = zeros(T,M);
nLWbal = floor(T/M);
[LW,nPerm] = genLyndon3(M,nLWbal);
nLW = size(LW,1);
indAvail = 1:nLW;


if mod(M,3)==0
    if nLWbal==(nLW-1)
        indMLeq = 1:nLWbal;
    else
        indMLeq = randperm(nLW-2,nLWbal)+1;
    end
else
    indMLeq = randperm(nLW-1,nLWbal)+1;
end
indAvail(indMLeq) = [];

countT = 1;
for i = 1:size(indMLeq,2)
    for m = 0:(nPerm(indMLeq(i))-1)
        V(countT,:) = circshift(LW(indMLeq(i),:),m);
        countT = countT+1;
    end
end

Trem = T-countT+1;
if mod(M,3)==0 && Trem>=(M/3)
    i = indAvail(end);
    indAvail(end) = [];
    for m = 0:(nPerm(i)-1)
        V(countT,:) = circshift(LW(i,:),m);
        countT = countT+1;
    end
    Trem = T-countT+1;
end

word = LW(1,:);
if mod(M,3)==0
    for i = 1:Trem
        alpha = floor(i/3);
        V(countT,:) = circshift(word,(i-1)*3+alpha);
        countT = countT+1;
    end
else
    for i = 1:Trem
        V(countT,:) = circshift(word,(i-1)*3);
        countT = countT+1;
    end
end

V = sortrows(V,'descend');

for t = 1:T
    calT{t} = (1:M).*(V(t,:));
    calT{t}(calT{t}==0) = [];
end

calT = permuteCalT(calT);
d = sequenceDegre(calT);

end