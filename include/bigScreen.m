function bigScreen(varargin)
% This function permits to change a figure position on a screen.
%
% Author: 
% name : Philippe Flores
% e-mail : flores.philipe@gmail.com
% github : github.com/philippeflores/fcm_ctflowhd

if nargin==0
    vScreen = [0 0 1 1];
end
if nargin==1
    if strcmpi(varargin{1},'full')
        vScreen = [0 0 1 1];
    elseif strcmpi(varargin{1},'halfL')
        vScreen = [0 0 0.5 1];
    elseif strcmpi(varargin{1},'halfR')
        vScreen = [0.5 0 0.5 1];
    elseif strcmpi(varargin{1},'halfTop')
        vScreen = [0 0.5 1 0.5];
    elseif strcmpi(varargin{1},'halfBot')
        vScreen = [0 0 1 0.5];
    else
        vScreen = varargin{1};
    end
end

set(gcf,"Units","normalized",'OuterPosition',vScreen)

end