function V = TboolMat(calT)

matCalT = reshape([calT{:}],3,[])';
M = max(matCalT(:));

V = zeros(size(calT,2),M);
for t = 1:size(calT,2)
    V(t,matCalT(t,:)) = 1;
end

end