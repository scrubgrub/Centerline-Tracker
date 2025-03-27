%% helper function to produce an isosurface of input volume for plotting
% Jade Lariviere | last modified Mar. 22, 2025

function [patchObj] = makePatch(Volume, Color)
% this function returns a preformatted patch object using the isosurface
% faces and vertices about dimensions and data. the Color argument just
% passes user input into the plot to change the face color, which will 
% accept any format accepted by the patch() Name-Value argument for
% shape.FaceColor (long name, short name, RGB, hexadecimal).

isofaces = Volume.faces;
isoverts = Volume.all_verts;

% plot isosurface in new figure
patchObj = patch('Faces',isofaces,'Vertices',isoverts);
    % formatting
    patchObj.FaceColor = Color; patchObj.EdgeColor = "none";
    patchObj.FaceAlpha = 0.25; patchObj.EdgeAlpha = patchObj.FaceAlpha;
end