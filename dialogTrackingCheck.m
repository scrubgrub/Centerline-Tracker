%% helper function to create a dialog box for centroid finalization
% Jade Lariviere | last modified Mar. 8, 2025

function answerYN = dialogTrackingCheck(InputFigure)
% formatting figure axes + etc.
p = InputFigure.CurrentAxes;
set(p,'Position',[0.1,0.1,0.5,0.8]);
xlabel('x'); ylabel('y'); zlabel('z'); view(3); grid minor;
legend('','Wire Centroids','Start','End','Location','best');
enableDefaultInteractivity(p);

% tweak figure with editable boxes; show pop-up & ask for inputs ----------
figure(InputFigure);
fieldYN = uicontrol('Style','checkbox','Units','Normalized',...
     'Position',[.65 .45 .3 .05],'Tag','myedit',...
     'String','Flip centroid order');
    annotation('textbox',[.65 .51 .3 .05],...
    'String','Check your final centroids, and proceed if OK.',...
    'LineStyle','none','Margin',0,'FontSize',10,...
    'VerticalAlignment','bottom');

OKbutton = uicontrol('Style','pushbutton','Units','Normalized',...
    'String','Continue to Export',...
    'Position',[.7 .3 .2 .075],'HorizontalAlignment','center');
    OKbutton.Callback = {@OK_callback,fieldYN}; % callback to execute code

uiwait(InputFigure); % prompt user input ----------------------------------

    function OK_callback(~,~,fieldData) % extract field input
        answerYN = logical(fieldData.Value); % return logical
        if answerYN % if user asked to flip centroids...
            fprintf('Centroid order flipped! Continuing to save...\n');
        else
            fprintf('Continuing to save...\n');
        end
        close(InputFigure);
    end
end