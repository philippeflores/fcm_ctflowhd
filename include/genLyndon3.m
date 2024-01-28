function [y,p] = genLyndon3(M,varargin)

if size(varargin,2)==0
    nT = nchoosek(M,3);
    T = TboolMat(fullCoupled(M));

    indGarde = zeros(nT,1);

    for i = 1:nT
        y = (1:nT)'.*(indGarde==1); y(y==0) = [];
        flag = 1;
        for j = 1:size(y,1)
            for m = 1:(M-1)
                temp = circshift(T(y(j),:),m);
                if flag>0 && size(unique([temp;T(i,:)],'rows'),1)==1
                    flag = 0;
                end
            end
        end
        if flag==1
            indGarde(i) = 1;
        end
    end

    y = T(indGarde==1,:);

elseif size(varargin,2)==1
    nML = varargin{1};
    y = zeros(nML,M);
    for n = 1:nML
        flag = 0;
        while flag==0
            temp = sort(randi(M,[1 3]),'ascend');
            if size(unique(temp),2)==3
                a = zeros(1,M); a(temp) = 1;
                temp = zeros(M);
                for m = 0:(M-1)
                    temp(m+1,:) = circshift(a,m);
                end
                temp = sortrows(temp,'descend');
                temp = temp(1,:);
                if sum(find(temp==1))~=6
                    y(n,:) = temp;
                    if size(unique(y(1:n,:),"rows"),1)==n
                        flag = 1;
                    end
                end
                if mod(M,3)==0
                    temp2 = find(temp==1); temp3 = ((0:2).*(M/3))+1;
                    if sum(temp2==temp3)==3
                        flag = 0;
                        y(n,:) = zeros(1,M);
                    end
                end
            end
        end
    end
    temp = zeros(1,M); temp(1:3) = 1;
    y = [temp; y];
    if mod(M,3)==0
        temp = zeros(1,M); temp(((0:2).*(M/3))+1) = 1;
        y = [y; temp];
    end
else

end

p = zeros(1,size(y,1));

for indP = 1:size(y,1)
    temp = zeros(M);
    temp(1,:) = y(indP,:);
    for m = 1:(M-1)
        temp(m+1,:) = circshift(temp(1,:),m);
    end
    p(indP) = size(unique(temp,'rows'),1);
end


end