function t = supportDist(X,indVar,I,varargin)
% This function permits to define the support of the distribution of
% variables in *indVar*. Because histograms are used afterwards, each
% interval is defined as a set of *I* points equally spread which create an
% interval support [a,b].
%
% Author: 
% name : Philippe Flores
% e-mail : flores.philipe@gmail.com
% github : github.com/philippeflores/fcm_ctflowhd

strMethod = 'none';

minFluo = -5;
maxFluo = 10;
minScat = 0;
maxScat = 1e+6;

if (rem(length(varargin),2)==1)
    error('Optional parameters should always go by pairs');
else
    for i=1:2:(length(varargin)-1)
        switch upper(varargin{i})
            case 'MINFLUO'
                minFluo = varargin{i+1};
            case 'MAXFLUO'
                maxFluo = varargin{i+1};
            case 'MINSCAT'
                minScat = varargin{i+1};
            case 'MAXSCAT'
                maxScat = varargin{i+1};
            case 'STRMETHOD'
                strMethod = varargin{i+1};
            case 'STRLABEL'
                strLabel = varargin{i+1};
            case 'EDGES'
                edges = varargin{i+1};
            otherwise
                error(['Unrecognized option: ''' varargin{i} '''']);
        end
    end
end

if strcmpi(strMethod,'none')

    if exist('strLabel','var')==1
        fprintf("\nWarning: the entry strLabel is not going to be used in the following.\nPlease change the strMethod argument if necessary.\n\n")
    end
    if exist('edges','var')==1
        fprintf("\nWarning: the entry edges is not going to be used in the following.\nPlease change the strMethod argument if necessary.\n\n")
    end

    M = size(indVar,2);
    t = cell(1,size(X,2));
    for m = 1:M
        t{indVar(m)} = linspace(-0.5,4.5,I);
    end

elseif strcmpi(strMethod,'fluo')

    if exist('strLabel','var')==1
        fprintf("\nWarning: the entry strLabel is not going to be used in the following.\nPlease change the strMethod argument if necessary.\n\n")
    end
    if exist('edges','var')==1
        fprintf("\nWarning: the entry edges is not going to be used in the following.\nPlease change the strMethod argument if necessary.\n\n")
    end

    M = size(indVar,2);
    t = cell(1,size(X,2));
    for m = 1:M
        [a,b] = findHistogramEdges(X(:,indVar(m)),'a0',minFluo,'b0',maxFluo);
    	t{indVar(m)} = linspace(a,b,I);
    end

elseif strcmpi(strMethod,'label')

    if exist('strLabel','var')==0
        error("To use the strMethod option 'label', you need to enter the strLabel parameter as an input of this function.\n")
    end
    if exist('edges','var')==1
        fprintf("Warning: the entry edges is not going to be used in the following.\n Please change the strMethod argument if necessary.\n")
    end

    M = size(indVar,2);
    t = cell(1,size(X,2));
    for m = 1:M
        if strcmpi(strLabel{m}(1:2),'FS') || strcmpi(strLabel{m}(1:2),'SS')
	        [a,b] = findHistogramEdges(X(:,indVar(m)),'a0',minScat,'b0',maxScat);
        else
	        [a,b] = findHistogramEdges(X(:,indVar(m)),'a0',minFluo,'b0',maxFluo);
        end
        t{indVar(m)} = linspace(a,b,I);
    end

elseif strcmpi(strMethod,'man')

    if exist('strLabel','var')==1
        fprintf("\nWarning: the entry strLabel is not going to be used in the following.\nPlease change the strMethod argument if necessary.\n\n")
    end

    if exist('edges','var')==0
        error("To use the strMethod option 'man', you need to enter the edges parameter as an input of this function.\n")
    end

    M = size(indVar,2);
    t = cell(1,size(X,2));
    for m = 1:M
        a = edges{m}(1);
        b = edges{m}(2);
        t{indVar(m)} = linspace(a,b,I);
    end
    
else
    error("The option %s for 'strMethod' is wrong. Please try again.\n",strMethod)
end

end

