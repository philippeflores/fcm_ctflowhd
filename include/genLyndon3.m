function [matLW,nPerm] = genLyndon3(M,varargin)

if size(varargin,2)==0
    Tfull = nchoosek(M,3);
    Vfull = TboolMat(fullCoupled(M));

    indKeep = zeros(Tfull,1);

    for i = 1:Tfull
        matLW = (1:Tfull)'.*(indKeep==1); matLW(matLW==0) = [];
        flag = 1;
        for j = 1:size(matLW,1)
            for m = 1:(M-1)
                temp = circshift(Vfull(matLW(j),:),m);
                if flag>0 && size(unique([temp;Vfull(i,:)],'rows'),1)==1
                    flag = 0;
                end
            end
        end
        if flag==1
            indKeep(i) = 1;
        end
    end

    matLW = Vfull(indKeep==1,:);

elseif size(varargin,2)==1
    nML = varargin{1};
    matLW = zeros(nML,M);
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
                    matLW(n,:) = temp;
                    if size(unique(matLW(1:n,:),"rows"),1)==n
                        flag = 1;
                    end
                end
                if mod(M,3)==0
                    temp2 = find(temp==1); temp3 = ((0:2).*(M/3))+1;
                    if sum(temp2==temp3)==3
                        flag = 0;
                        matLW(n,:) = zeros(1,M);
                    end
                end
            end
        end
    end
    temp = zeros(1,M); temp(1:3) = 1;
    matLW = [temp; matLW];
    if mod(M,3)==0
        temp = zeros(1,M); temp(((0:2).*(M/3))+1) = 1;
        matLW = [matLW; temp];
    end
else

end

nPerm = zeros(1,size(matLW,1));

for indLW = 1:size(matLW,1)
    temp = zeros(M);
    temp(1,:) = matLW(indLW,:);
    for m = 1:(M-1)
        temp(m+1,:) = circshift(temp(1,:),m);
    end
    nPerm(indLW) = size(unique(temp,'rows'),1);
end


end