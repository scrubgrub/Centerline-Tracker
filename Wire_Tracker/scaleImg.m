%% helper function to upscale image to higher resolution
% Jade Lariviere | last modified Nov. 10, 2025

function [ScaledImage,ScaledDims] = scaleImg(Image, scaleFactor, interpMethod)
% the function wraps the MATLAB imresize function to change the Image's 
% size and performing other processing steps for refinement purposes.
% scaleFactor determines the proportional size of the ScaledImage to the
% input Image, and the argument interpMethod is directly passed to the
% "method" argument of imresize. Returns dimensions ScaledDims. ===========

filtSigma = 0.5; % standard deviation in gaussian filter; higher = more blur

% scale and process image
fin_img = imresize(Image,scaleFactor,interpMethod);
    fin_img = imgaussfilt(fin_img,filtSigma,"Padding","circular");

ScaledImage = fin_img;
ScaledDims = size(fin_img);

% plotting for debugging ==================================================
%{
img = Image; fin_img = ScaledImage;
max_val = 2*max(fin_img,[],'all');
valRange = [0 max_val*1.1]; % intensity scaling for plotting
% highlight center pixel of images (guess coordinate)
img_dim = size(img);
    img(ceil(img_dim(1)/2),ceil(img_dim(2)/2)) = max_val;
scaled_dim = size(fin_img);
    fin_img(ceil(scaled_dim(1)/2),ceil(scaled_dim(2)/2)) = max_val;
    
figure(21); % transpose images for plotting (MATLAB hates being normal)
    subplot(2,1,1); imshow(img',valRange); set(gca,"YDir","normal");
        title('Original Image (Center Highlighted)');
    subplot(2,1,2); imshow(fin_img',valRange); set(gca,"YDir","normal");
        title('Upscaled Image (Center Highlighted)');
%}
end