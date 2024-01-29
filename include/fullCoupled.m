function calTfull = fullCoupled(M)

calTfull = cell(1,nchoosek(M,3));
countT = 1;
for j = 1:(M-2)
	for k = (j+1):(M-1)
		for l = (k+1):M
			calTfull{countT} = [j k l];
			countT = countT+1;
		end
	end
end

end