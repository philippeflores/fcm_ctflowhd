function scatterColor(X,Y,varargin)

limitX = [-0.5 4.5];
limitY = [-0.5 4.5];
vecScreen = [0 0 500 600];
numBins = 200;
colMap = jet(256); colMap = [colMap(80:2:200,:); colMap(201:end-20,:)];
colMap = [0.8, 0.8, 0.85;colMap];
threshHisto = 0.05;
strMarker = 'square';
sizeMarker = 10;
stepCloud = 1;
axFontSize = 15;
newFig = 'on';

if (rem(length(varargin),2)==1)
    error('Optional parameters should always go by pairs');
else
    for i=1:2:(length(varargin)-1)
        switch upper(varargin{i})
            case 'LIMITX'
                limitX = varargin{i+1};
            case 'LIMITY'
                limitY = varargin{i+1};
            case 'VECSCREEN'
                vecScreen = varargin{i+1};
            case 'NUMBINS'
                numBins = varargin{i+1};
            case 'COLMAP'
                colMap = varargin{i+1};
            case 'THRESHHISTO'
                threshHisto = varargin{i+1};
            case 'STRMARKER'
                strMarker = varargin{i+1};
            case 'SIZEMARKER'
                sizeMarker = varargin{i+1};
            case 'STEPCLOUD'
                stepCloud = varargin{i+1};
            case 'AXFONTSIZE'
                axFontSize = varargin{i+1};
			case 'NEWFIG'
                newFig = varargin{i+1};
            otherwise
                error(['Unrecognized option: ''' varargin{i} '''']);
        end
    end
end

if sum(size(X)==size(Y))~=size(size(X),2)
	error("X and Y have not the same size\n")
end

h = histcounts2(X,Y,'NumBins',numBins,'XBinLimits',limitX,'YBinLimits',limitY);
maxH = max(max(h));
if threshHisto>1
	h(h<threshHisto) = 0;
else
	h(h<(threshHisto*maxH))= 0;
end

X = X(1:stepCloud:end);
Y = Y(1:stepCloud:end);
nPoints = size(X,1);

matColor = zeros(nPoints,3);
for i = 1:nPoints
	temp = ceil(numBins*[(X(i)-limitX(1))/(limitX(2)-limitX(1)) (Y(i)-limitY(1))/(limitY(2)-limitY(1))]);
    temp(temp<1) = 1; temp(temp>size(h,1)) = size(h,1);
	temp = h(temp(1),temp(2));
	if temp>0
		matColor(i,:) = colMap(ceil(temp*size(colMap,1)/maxH),:);
	else
		matColor(i,:) = colMap(1,:);
	end
end

[matColor,indPerm] = sortrows(matColor,'ascend');
X = X(indPerm);
Y = Y(indPerm);

if strcmpi(newFig,'on')==1
	figure
	set(gcf,'OuterPosition',vecScreen)
	scatter(X,Y,sizeMarker,matColor,"filled",'Marker',strMarker)
	colormap(colMap)
	xlim(limitX), ylim(limitY)
	ax = gca; ax.FontSize = axFontSize;
else
	scatter(X,Y,sizeMarker,matColor,"filled",'Marker',strMarker)
	colormap(colMap)
end

end