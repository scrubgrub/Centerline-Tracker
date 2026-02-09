%% helper function to create a dialog box for centroid finalization
% Jade Lariviere | last modified Dec. 7, 2025

function [answer] = dialogTrackCheck(InputFigure,figSet,TrackStats)
% setting up dynamic labels -----------------------------------------------
num_line = numel(findobj(gca,'Type','Text'));
[labels{1:num_line}] = deal('');
    labels = [{''},{'First Track'},{'Start'},labels(:)',{'End'}];

track_stats = [(1:numel(TrackStats))-1;TrackStats]';
s_track = compose("Track: %d points",track_stats(1,2));

% formatting figure axes + etc. -------------------------------------------
p = InputFigure.CurrentAxes;
set(p,'Position',[0.1,0.1,0.5,0.8]);
xlabel('x'); ylabel('y'); zlabel('z'); grid minor; axis equal; axis tight;
if isempty(figSet), view(3); 
    else, view(figSet{1}), set(gcf,"WindowState",figSet{2});
end
legend(labels,'Location','best');
enableDefaultInteractivity(p);

% tweak figure with editable boxes; show pop-up & ask for inputs ----------
figure(InputFigure);
fieldYN = uicontrol('Style','checkbox','Units','Normalized',...
     'Position',[.65 .45 .3 .05],'Tag','myedit',...
     'String','Flip main centroid sequence');
    annotation('textbox',[.65 .51 .3 .05],...
    'String','Check your final centroids, and proceed if OK.',...
    'LineStyle','none','Margin',0,'FontSize',10,...
    'VerticalAlignment','bottom');
    
if length(TrackStats) > 1 % stats includes branch track point counts
    s_branch = compose("Branch %d: %d points",track_stats(2:end,:));
    annotation('textbox',[.65 .65 .3 .3],'String', [s_track;s_branch],...
        'LineStyle','none','Margin',0,'FontSize',10,...
        'VerticalAlignment','bottom','FontWeight','bold');
else % just one track to display number of points
    annotation('textbox',[.65 .7 .3 .3],'String', s_track,...
        'LineStyle','none','Margin',0,'FontSize',10,...
        'VerticalAlignment','bottom','FontWeight','bold');
end

OKbutton = uicontrol('Style','pushbutton','Units','Normalized',...
    'String','Continue to Export',...
    'Position',[.7 .3 .2 .075],'HorizontalAlignment','center');
    OKbutton.Callback = {@OK_callback,fieldYN}; % callback to execute code

uiwait(InputFigure); % prompt user input ----------------------------------

    function OK_callback(~,~,fieldData) % extract field input
        answer = logical(fieldData.Value); % return logical
        if answer % if user asked to flip centroids...
            fprintf('Centroid order flipped! Continuing to save...\n');
        else
            fprintf('Continuing to save...\n');
        end
        close(InputFigure);
    end
end
