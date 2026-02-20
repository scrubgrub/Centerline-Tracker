%% helper function to create a dialog box for data viewing + input values
% Jade Lariviere | last modified Oct. 20, 2025

function [answer,figSet] = dialogWireStart(InputFigure,viewAngle,varargin)
% formatting figure axes + etc. -------------------------------------------
p = InputFigure.CurrentAxes;
set(p,'Position',[0.1,0.1,0.5,0.8]);
xlabel('x'); ylabel('y'); zlabel('z'); grid minor; axis equal; axis tight;
if isempty(viewAngle), view(3); else, view(viewAngle); end
enableDefaultInteractivity(p);

% tweak figure with editable boxes; show pop-up, ask for inputs -----------
figure(InputFigure);
switch length(varargin)
    case 0 % Manual UI ----------------------------------------------------
    Pts = {};
    fieldP = uicontrol('Style','Edit','Units','Normalized',...
        'Position',[.65 .65 .3 .05],'Tag','myedit','String','100 10 175');
        annotation('textbox',[.65 .71 .3 .05],...
        'String','Enter [X Y Z] coordinates of start point.',...
        'LineStyle','none','Margin',0,'FontSize',10,...
        'VerticalAlignment','bottom');
    fieldN = uicontrol('Style','Edit','Units','Normalized',...
        'Position',[.65 .45 .3 .05],'Tag','myedit','String','0.5 -1 1');
        annotation('textbox',[.65 .51 .3 .05],...
        'String','Enter [X Y Z] vector components of start direction.',...
        'LineStyle','none','Margin',0,'FontSize',10,...
        'VerticalAlignment','bottom');

    case 1 % Automatic UI -------------------------------------------------
    Pts = varargin{1};
    pts_list = strcat(string(1:size(Pts,1))',':  [',num2str(Pts),']');
    fieldP = uicontrol('Style','popupmenu','Units','Normalized',...
        'String',pts_list,'Position',[.65 .65 .3 .05]);
        annotation('textbox',[.65 .71 .3 .05],...
        'String','Select the auto-detected [X Y Z] end point.',...
        'LineStyle','none','Margin',0,'FontSize',10,...
        'VerticalAlignment','bottom');
    fieldN = uicontrol('Style','Edit','Units','Normalized',...
        'Position',[.65 .45 .3 .05],'Tag','myedit','String','0.5 -1 1');
        annotation('textbox',[.65 .51 .3 .05],...
        'String','Enter [X Y Z] vector components of start direction.',...
        'LineStyle','none','Margin',0,'FontSize',10,...
        'VerticalAlignment','bottom');
    otherwise, error('dialogWireStart() has too many input arguments.');
end % ---------------------------------------------------------------------

OKbutton = uicontrol('Style','pushbutton','Units','Normalized',...
    'String','Start Tracking','Position',[.7 .3 .2 .075],...
    'HorizontalAlignment','center');
    OKbutton.Callback = {@OK_callback,[fieldP; fieldN],Pts}; % execute code

uiwait(InputFigure); % prompt user input ----------------------------------

    function OK_callback(~,~,fieldData,points) % extract field inputs
        switch length(varargin)
            case 0
            ans_P = str2num(fieldData(1).String,"Evaluation","restricted");
            case 1
            ans_P = points(fieldData(1).Value,:);
        end
        ans_N = str2num(fieldData(2).String,"Evaluation","restricted");
        answer = [ans_P; ans_N];
        
        fprintf('Wire start and direction entered!\n');
        [figSet{1}(1),figSet{1}(2)] = view; % save fig settings for others
        figSet{2} = get(gcf,"WindowState");
        close(InputFigure);
    end
end
