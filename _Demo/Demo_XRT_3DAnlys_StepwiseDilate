%------------------------------------------------------------------------------------------------
% Code written by Frederik Doerr, Aug 2020 (MATLAB R2019b)
% Application: Demo for XRT_3DAnlys_StepwiseDilate
% 
% frederik.doerr(at)strath.ac.uk | CMAC (http://www.cmac.ac.uk/)
% https://github.com/frederik-d/XRT_3DAnlys


clear all
close all
clc

%% Create an example V
center1 = 40;
center2 = center1 + 40;
dist = sqrt(1*(2*center1)^2);
radius = dist/4 * 1.4;
lims = [0 ceil(center2+1.8*radius)];
[x,y,z] = meshgrid(lims(1):lims(2));
bw1 = sqrt((x-center1).^2 + (y-center1).^2 + ...
   (z-center1).^2) <= radius;
bw2 = sqrt((x-center2).^2 + (y-center2).^2 + ...
   (z-center2).^2) <= radius;

V = bw1 | bw2;

figure, isosurface(x,y,z,V,0.5), axis equal, title('BW')
xlabel x, ylabel y, zlabel z
xlim(lims), ylim(lims), zlim(lims)
view(3), camlight, lighting gouraud

% Growth within full volume
V_ROI_1 = ones(size(V));

% Growth within specified volume
idx_c_Row = round(size(V,1)/2);
idx_c_Img = round(size(V,3)/2);
V_ROI_2 = V_ROI_1;
V_ROI_2(idx_c_Row-2:idx_c_Row+2,:,:) = 0;


figure;imshow(V(:,:,idx_c_Img))
figure;imshow(V_ROI_1(:,:,idx_c_Img))
figure;imshow(V_ROI_2(:,:,idx_c_Img))


%% Stepwise 3D dilation (3d_stepwisedilate)
Opt.numDilationSteps = 10;
Opt.method = '3d_stepwisedilate';
Opt.ExpShorthand = 'Spheres (3d_stepwisedilate)';
Opt.Print_Ctrl_On = false; 

[V_dil,~] = XRT_3DAnlys_StepwiseDilate(V,V_ROI_2,Opt);


figure;imshow(V_dil(:,:,idx_c_Img))

figure, isosurface(x,y,z,V_dil,0.5), axis equal, title('BW')
xlabel x, ylabel y, zlabel z
xlim(lims), ylim(lims), zlim(lims)
view(3), camlight, lighting gouraud


%% Stepwise dilation (3d_bwmorph_thicken)
Opt.numDilationSteps = 10;
Opt.method = '3d_bwmorph_thicken';
Opt.ExpShorthand = 'Spheres (3d_bwmorph_thicken)';
Opt.Print_Ctrl_On = true;

[V_thi,~] = XRT_3DAnlys_StepwiseDilate(V,V_ROI_1,Opt);

figure;imshow(V_thi(:,:,idx_c_Img))

figure, isosurface(x,y,z,V_thi,0.5), axis equal, title('BW')
xlabel x, ylabel y, zlabel z
xlim(lims), ylim(lims), zlim(lims)
view(3), camlight, lighting gouraud
