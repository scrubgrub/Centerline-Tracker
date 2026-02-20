%% creates a path with bifurcation at regular intervals
% by Warren Boschen
clear; clc; close all

length = 5; segments = 500;
t = linspace(1, length, segments);
x = t; y = t; z = t;
all_points = [x', y', z'];

figure;
plot3(x,y,z,'b','LineWidth',2);
hold on; grid on; axis equal;
xlabel('X'); ylabel('Y'); zlabel('Z');

branches1 = [100,200,300,400];
x_rel1 = -0.2*t;
y_rel1 = 0.2*t;
z_rel1 = 0.2*t;

for k=1:size(branches1,2)
    idx = branches1(k);

    x0 = x(idx); y0 = y(idx); z0 = z(idx);
    x_branch = x_rel1 + x0 - x_rel1(end);
    y_branch = y_rel1 + y0 - y_rel1(end);
    z_branch = z_rel1 + z0 - z_rel1(end);

    all_points = [all_points; [x_branch', y_branch', z_branch']];
    plot3(x_branch, y_branch, z_branch, 'r', 'LineWidth', 1.5);
    plot3(x0, y0, z0, 'ko', 'MarkerSize', 6, 'MarkerFaceColor', 'k');
end

branches2 = [50,150,250,350,450];
x_rel2 = 0.2*t;
y_rel2 = -0.2*t;
z_rel2 = 0.2*t;

for k=1:size(branches2,2)
    idx = branches2(k);

    x0 = x(idx); y0 = y(idx); z0 = z(idx);
    x_branch = x_rel2 + x0 - x_rel2(end);
    y_branch = y_rel2 + y0 - y_rel2(end);
    z_branch = z_rel2 + z0 - z_rel2(end);

    all_points = [all_points; [x_branch', y_branch', z_branch']];
    plot3(x_branch, y_branch, z_branch, 'g', 'LineWidth', 1.5);
    plot3(x0, y0, z0, 'ko', 'MarkerSize', 6, 'MarkerFaceColor', 'k');
end

branches3 = [25,125,225,325,425];
x_rel3 = -0.2*t;
y_rel3 = 0.2*t;
z_rel3 = -0.2*t;

for k=1:size(branches3,2)
    idx = branches3(k);

    x0 = x(idx); y0 = y(idx); z0 = z(idx);
    x_branch = x_rel3 + x0 - x_rel3(end);
    y_branch = y_rel3 + y0 - y_rel3(end);
    z_branch = z_rel3 + z0 - z_rel3(end);

    all_points = [all_points; [x_branch', y_branch', z_branch']];
    plot3(x_branch, y_branch, z_branch, 'r', 'LineWidth', 1.5);
    plot3(x0, y0, z0, 'ko', 'MarkerSize', 6, 'MarkerFaceColor', 'k');
end

branches4 = [75,175,275,375,475];
x_rel4 = 0.2*t;
y_rel4 = -0.2*t;
z_rel4 = -0.2*t;

for k=1:size(branches4,2)
    idx = branches4(k);

    x0 = x(idx); y0 = y(idx); z0 = z(idx);
    x_branch = x_rel4 + x0 - x_rel4(end);
    y_branch = y_rel4 + y0 - y_rel4(end);
    z_branch = z_rel4 + z0 - z_rel4(end);

    all_points = [all_points; [x_branch', y_branch', z_branch']];
    plot3(x_branch, y_branch, z_branch, 'g', 'LineWidth', 1.5);
    plot3(x0, y0, z0, 'ko', 'MarkerSize', 6, 'MarkerFaceColor', 'k');
end

min_coords = min(all_points, [], 1); max_coords = max(all_points, [], 1);
scaled_points = (all_points - min_coords)./(max_coords - min_coords);
scaled_points = 223*scaled_points + 1;

vol_size = 224;
volume = zeros(vol_size);

indices = round(scaled_points);
indices = max(min(indices, vol_size), 1);

for i =1:size(indices,1)
    volume(indices(i,1), indices(i,2), indices(i,3)) = 1;
end

SE = strel('sphere', 3);
volume = imdilate(volume, SE); volshow(volume);
