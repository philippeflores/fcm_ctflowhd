function plotTsne(Y,indVar,X,strLabel,varargin)

stepCloud = 50;
M = size(indVar,2);
xlimMan = [-0.5 4.5];
colMax = 4;
strScreen = 'halfTop';
flagScreen = 0;

if (rem(length(varargin),2)==1)
    error('Optional parameters should always go by pairs');
else
    for i=1:2:(length(varargin)-1)
        switch upper(varargin{i})
            case 'STEPCLOUD'
                stepCloud = varargin{i+1};
            case 'STRSCREEN'
                flagScreen = 1;
                strScreen = varargin{i+1};
            case 'XLIMMAN'
                xlimMan = varargin{i+1};
            case 'COLMAX'
                colMax = varargin{i+1};
            otherwise
                error(['Unrecognized option: ''' varargin{i} '''']);
        end
    end
end

if strcmpi(xlimMan,'auto')
    xlimMan = [inf -inf];
end

col = min([colMax M]);
row = ceil(M/col);

figure
cmapJet = [50,136,189; 102,194,165; 171,221,164; 230,245,152; 255,255,191; 254,224,139; 253,174,97; 244,109,67; 213,62,79]/255;
if flagScreen==0, bigScreen('demiTop'), else, bigScreen(strScreen); end
ax = cell(1,M);
for m = 1:M
    ax{m} = subplot(row,col,m);
    scatter(Y(1:stepCloud:end,2),Y(1:stepCloud:end,1),100,X(1:stepCloud:end,indVar(m)),'.')
    if xlimMan(1)~=inf
        clim(xlimMan)
    end
    colormap(cmapJet)
    colorbar
    axe = gca; axe.FontSize = 15;
    title(strLabel{m},'FontSize',30)
end
linkaxes([ax{:}],'xy')

end

