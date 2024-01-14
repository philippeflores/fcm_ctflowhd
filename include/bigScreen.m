function bigScreen(varargin)

if nargin==0
    vScreen = [0 0 1 1];
end
if nargin==1
    if strcmpi(varargin{1},'full')
        vScreen = [0 0 1 1];
    elseif strcmpi(varargin{1},'demiL')
        vScreen = [0 0 0.5 1];
    elseif strcmpi(varargin{1},'demiR')
        vScreen = [0.5 0 0.5 1];
    elseif strcmpi(varargin{1},'demiTop')
        vScreen = [0 0.5 1 0.5];
    elseif strcmpi(varargin{1},'demiBot')
        vScreen = [0 0 1 0.5];
    else
        vScreen = varargin{1};
    end
end

set(gcf,"Units","normalized",'OuterPosition',vScreen)

end