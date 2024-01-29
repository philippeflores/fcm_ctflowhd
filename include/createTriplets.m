function [calTout,Tout] = createTriplets(M,indVar,strStrategy,varargin)
% This function returns a coupling calT from a set of variables in indVar.
% To choose triplets, a coupling strategy must be entered as an input of
% this function. Please refer to the thesis manuscript of the author for
% more information on coupling strategies.
%
% Author: 
% name : Philippe Flores
% e-mail : flores.philipe@gmail.com
% github : github.com/philippeflores/fcm_ctflowhd

if (rem(length(varargin),2)==1)
    error('Optional parameters should always go by pairs');
else
    for i=1:2:(length(varargin)-1)
        switch upper(varargin{i})
            case 'CALTIN'
                calTin = varargin{i+1};
            case 'TIN'
                Tin = varargin{i+1};
            otherwise
                error(['Unrecognized option: ''' varargin{i} '''']);
        end
    end
end

Tmax = nchoosek(M,3);
calTmax = cell(1,Tmax);
count = 1;
for i = 1:M
	for j = i+1:M
		for k = j+1:M
			calTmax{count} = [indVar(i) indVar(j) indVar(k)];
			count = count+1;
		end
	end
end

if strcmpi(strStrategy,'full')

    if exist('calTin','var')==1
        fprintf("\nWarning: the entry 'calTin' is not going to be used in the following.\nPlease change the strStrategy argument if necessary.\n\n")
    end
    if exist('Tin','var')==1
        fprintf("\nWarning: the entry 'Tin' is not going to be used in the following.\nPlease change the strStrategy argument if necessary.\n\n")
    end

	Tout = Tmax;
	calTout = calTmax;
elseif strcmpi(strStrategy,'+1')

    if exist('calTin','var')==1
        fprintf("\nWarning: the entry 'calTin' is not going to be used in the following.\nPlease change the strStrategy argument if necessary.\n\n")
    end
    if exist('Tin','var')==1
        fprintf("\nWarning: the entry 'Tin' is not going to be used in the following.\nPlease change the strStrategy argument if necessary.\n\n")
    end
    
    Tout = M; calTout = cell(1,Tout);

    for nTrip = 1:Tout-2
        if nTrip==1
			calTout{1} = [indVar(1) indVar(2) indVar(3)];
			calTout{2} = [indVar(1) indVar(2) indVar(M)];
			calTout{3} = [indVar(1) indVar(M-1) indVar(M)];
		else
			calTout{nTrip+2} = [indVar(nTrip) indVar(nTrip+1) indVar(nTrip+2)];
        end
    end
elseif strcmpi(strStrategy,'+2')

    if exist('calTin','var')==1
        fprintf("\nWarning: the entry 'calTin' is not going to be used in the following.\nPlease change the strStrategy argument if necessary.\n\n")
    end
    if exist('Tin','var')==1
        fprintf("\nWarning: the entry 'Tin' is not going to be used in the following.\nPlease change the strStrategy argument if necessary.\n\n")
    end
    
    Tout = floor(M/2); calTout = cell(1,Tout);
    calTout{1} = [1 2 3];
    for t = 2:Tout
        calTout{t} = calTout{t-1}+2;
    end
    for t = 1:Tout
        for m = 1:3
            if calTout{t}(m)>M, calTout{t}(m) = calTout{t}(m)-M; end
        end
        calTout{t} = indVar(sort(calTout{t}));
    end

    V = reshape([calTout{:}],3,[])';
    [~,indSort] = sortrows(V);
    calTout = calTout(indSort);

elseif strcmpi(strStrategy,'rng')

    if exist('calTin','var')==1
        fprintf("\nWarning: the entry 'calTin' is not going to be used in the following.\nPlease change the strStrategy argument if necessary.\n\n")
    end
    if exist('Tin','var')==0
        error("To use the strStrategy option 'rng', you need to enter the number of triplets parameter 'Tin' as an input of this function.\n")
    end
    
    if Tin>Tmax
        error("The number of triplets 'Tin' entered (%d) is superior to the maximal number of triplets possible (%d). Try again with another 'Tin'.", Tin, Tmax)
    end
    if Tin<floor(M/2)
        error("The number of triplets 'Tin' entered (%d) is inferior to the least amount of triplets possible (%d). Try again with another 'Tin'.", Tin, floor(M/2))
    end

    Tout = Tin;
	calTout = calTmax(randperm(Tmax,Tout));
    while testStratCoupling(calTout,M)~=1
        calTout = calTmax(randperm(Tmax,Tout));
    end
    
    V = reshape([calTout{:}],3,[])';
    [~,indSort] = sortrows(V);
    calTout = calTout(indSort);

elseif strcmpi(strStrategy,'bal')

    if exist('calTin','var')==1
        fprintf("\nWarning: the entry 'calTin' is not going to be used in the following.\nPlease change the strStrategy argument if necessary.\n\n")
    end
    if exist('Tin','var')==0
        error("To use the strStrategy option 'rng', you need to enter the number of triplets parameter 'Tin' as an input of this function.\n")
    end
    
    if Tin>Tmax
        error("The number of triplets 'Tin' entered (%d) is superior to the maximal number of triplets possible (%d). Try again with another 'Tin'.", Tin, Tmax)
    end
    if Tin<floor(M/2)
        error("The number of triplets 'Tin' entered (%d) is inferior to the least amount of triplets possible (%d). Try again with another 'Tin'.", Tin, floor(M/2))
    end

    Tout = Tin;
    calTout = balancedCoupling(M,Tout);
    for t = 1:Tout
        calTout{t} = indVar(calTout{t});
    end


elseif strcmpi(strStrategy,'man')

    if exist('Tin','var')==1
        fprintf("\nWarning: the entry 'Tin' is not going to be used in the following.\nPlease change the strStrategy argument if necessary.\n\n")
    end
    if exist('calTin','var')==0
        error("To use the strStrategy option 'man', you need to enter the subset of triplets parameter 'calTin' as an input of this function.\n")
    end
    
    if testStratCoupling(calTin,M)==1
        calTout = calTin;
    else
        error("The coupling 'calTin' is not valid and cannot be used in the following. Please enter a new valid 'calTin'.\n")
    end

else
    
    flag = 0;
    strStrategy = char(strStrategy);

    if strcmpi(strStrategy(end-2:end),'bal')
        tau = str2num(strStrategy(1:end-3));
        Tout = round(tau*Tmax);
        if Tout<=Tmax && Tout>=floor(M/2)
            calTout = balancedCoupling(M,Tout);
            for t = 1:Tout
                calTout{t} = indVar(calTout{t});
            end
            flag = 1;
        else
            error("The strStrategy '%s' is not an available triplet strategy.\n",strStrategy)
        end
    end
    
    if flag==0
        tau = str2num(strStrategy);
        Tout = round(tau*Tmax);
        if Tout<=Tmax && Tout>=floor(M/2)
            calTout = calTmax(randperm(Tmax,Tout));
            while testStratCoupling(calTout,M)~=1
                calTout = calTmax(randperm(Tmax,Tout));
            end
            flag = 1;
        else
            error("The strStrategy '%s' is not an available triplet strategy.\n",strStrategy)
        end
    end

    if flag==0
        error("The strStrategy '%s' is not an available triplet strategy.\n",strStrategy)
    end
end

end