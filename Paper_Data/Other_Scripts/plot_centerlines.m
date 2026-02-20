%% plot centerlines
% author: Jade Lariviere
    % date: Feb. 20, 2026
% to make a graphic for presentation purposes! the VMTK, ORTHO, and
% wire_mask variables should already be in the workspace when running this
% code. ===================================================================

title_label = set_name + ID{n} + geometry(m);

% convert original volume mask to isosurface for elegant display
[x,y,z] = ndgrid(1:size(wire_mask,1),1:size(wire_mask,2),1:size(wire_mask,3));
[isoVol.faces,isoVol.all_verts] = isosurface(x,y,z,wire_mask,0);

% assemble figure from ORTHO, VMTK, and TRUTH
figure;
    patchObj = makePatch(isoVol,[0.5 0.6 0.5]); % change color here
    view(3); axis equal; grid minor; axis ij;
    hold on;
    if exist("TRUTH", 'var')
        plot3(TRUTH(:,1),TRUTH(:,2),TRUTH(:,3),"-k","LineWidth",3);
    end
    plot3(ORTHO(:,1),ORTHO(:,2),ORTHO(:,3),".-r","LineWidth",1.5);
    plot3(VMTK(:,1),VMTK(:,2),VMTK(:,3),".-b","LineWidth",1.5);
    hold off;
    if exist("TRUTH",'var')
        legend("","TRUTH","ORTHO","VMTK","Location","best");
    else
        legend("","ORTHO","VMTK","Location","best");
    end
    xlabel("X"); ylabel("Y"); zlabel("Z");
    title(title_label);