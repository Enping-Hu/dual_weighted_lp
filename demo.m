clear all;
close all;
clc;


filename='D:\our\dataset\LOL\1.png';
outputDir='D:\our\our_after\final\check';
im = im2double(imread(filename)); 

% Estimate reflectance and illumination layers using our model
addpath(genpath('Third_codes'));
tic;
[L,R,N]=dual_weighted_lp(im);
elapsedTime = toc;
fprintf('The function reweighted_after took %f seconds to execute.\n', elapsedTime);


% Enhanced images
gamma=2.2;
hsv=rgb2hsv(im);
L_gamma=L.^(1/gamma);
S_gamma=R .* L_gamma;
hsv(:,:,3)=real(S_gamma);
enhance = hsv2rgb(hsv);

imwrite(enhance,fullfile(outputDir,'E.png'));
