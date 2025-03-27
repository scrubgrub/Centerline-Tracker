%% helper function to create a dialog box for data viewing + input values
% Jade Lariviere | last modified Mar. 8, 2025

function answer = dialogWireStart(InputFigure)
% formatting figure axes + etc.
p = InputFigure.CurrentAxes;
set(p,'Position',[0.1,0.1,0.5,0.8]);
enableDefaultInteractivity(p);
xlabel('x'); ylabel('y'); zlabel('z'); view(3); grid minor;
enableDefaultInteractivity(p);

% tweak figure with editable boxes; show pop-up & ask for inputs ----------
figure(InputFigure);
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

OKbutton = uicontrol('Style','pushbutton','Units','Normalized',...
    'String','Start Tracking','Position',[.7 .3 .2 .075],...
    'HorizontalAlignment','center');
    OKbutton.Callback = {@OK_callback,[fieldP; fieldN]}; % execute code

uiwait(InputFigure); % prompt user input ----------------------------------

    function OK_callback(~,~,fieldData) % extract field inputs
        ans_P = str2num(fieldData(1).String,"Evaluation","restricted");
        ans_N = str2num(fieldData(2).String,"Evaluation","restricted");
        answer = [ans_P; ans_N];
        
        fprintf('Wire start and direction entered!\n');
        close(InputFigure);
    end
end
