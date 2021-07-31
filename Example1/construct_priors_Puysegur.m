%                       BayGrav3D, Version 1.0 
%       © 2021 Erin Hightower, California Institute of Technology
%   ----------------------------------------------------------------------
%   Supported in part by the US National Science Foundation
%   (http://www.nsf.gov).
%
%   Free for non-commercial academic use only.
%
%   ----------------------------------------------------------------------
%   This file [construct_prior_Puysegur.m] is part of BayGrav3D.

%   BayGrav3D is free software: you can redistribute it and/or modify it
%   under the terms of the GNU Lesser General Public License as published
%   by the Free Software Foundation, version 3 of the License.
%
%   BayGrav3D is distributed in the hope that it will be useful, but
%   WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTIBILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser
%   General Public License for more details.
%
%   You should have received a copy of the GNU Lesser General Public
%   License along with BayGrav3D. If not, see http://www.gnu.org/licenses/
%
%   ----------------------------------------------------------------------
%   Copyright July 2021, Erin Hightower, California Institute of Technology
%   ALL RIGHTS RESERVED. United States Government Sponsoring Acknowledged.
%
%   Please cite:
%   
%       AND
%   Hightower et al. (2020), A Bayesian 3-D linear gravity inversion for
%   complex density distributions: application to the Puysegur subduction
%   system: Geophysical Journal International, v. 223, p.1899-1918, 
%   doi:10.1093/gji/ggaa425.

% ======================================================================= %
% construct_priors_Puysegur.m
% ---------------------------
% This script constructs a low resolution version of the prior model for 
% the gravity inversion of the Puysegur subduction system. 
% 
% *** This is an example for how a prior could be constructed. Every prior
% is very specific to the particular case study at hand, so design whatever
% prior you need or have data for. 
%
% The required OUTPUTs that is needed for the gravity inversion are:
%   mu_abs:   vector of the prior absolute density of each prism
%   mu_pr:    vector of the prior relative density of each prism
%   sigma_pr: vector of the prior standard deviation of each prism
%   Cm:       prior covariance matrix 
%   layer_mean: vector of the mean density that was removed from each prism
%
% Other OUTPUT:
%   layers: cell array of layers used in the prior model
%   D:      vector of nearest neighbor distances of prisms to data points
%   prism_density_stdev: prior prism standard deviations before being
%                        weighted by the distance from data
% ======================================================================= %

filepath = 'BayGrav3D/Example1/';
modelname = 'Example1';

% Specify mesh and data to be used in constructing the prior
input_mesh      = [filepath modelname '_mesh.mat'];
input_bathy     = [filepath modelname '_data_bathy_preprocessed.txt'];
input_structure = [filepath modelname '_data_structure_preprocessed.txt'];
input_density   = [filepath modelname '_data_density_preprocessed.txt'];

% ======================================================================= %
% Load mesh and data
% ==================
load(input_mesh)

data_bathy     = readtable(input_bathy,'ReadVariableNames',1);
data_structure = readtable(input_structure,'ReadVariableNames',1);
data_density   = readtable(input_density,'ReadVariableNames',1);

data_bathy     = table2array(data_bathy);
data_structure = table2array(data_structure);
data_density   = table2array(data_density);

data_horizons = [data_structure; data_bathy];

data_horizons(data_horizons(:,1) > Xmax,:) = [];
data_horizons(data_horizons(:,2) > Ymax,:) = [];
data_density(data_horizons(:,1) > Xmax,:)  = [];
data_density(data_horizons(:,2) > Ymax,:)  = [];

% ======================================================================= %
% Interpolate horizon interfaces from structural data and map to prisms
% =====================================================================
layerID = data_horizons(:,4);

% Define layer interfaces using the unique layerID's of the input data
h1 = [data_horizons(layerID == 0,1:3); data_horizons(layerID == 1,1:3); data_horizons(layerID == 2,1:3)]; % seafloor
h2 = [data_horizons(layerID == 7,1:3); data_horizons(layerID == 9,1:3)]; % AU basement, including decollement
h3 = [data_horizons(layerID == 8,1:3)]; % PA basement
h4 = [data_horizons(layerID == 10,1:3)]; % AU moho
h5 = [data_horizons(layerID == 11,1:3)]; % PA moho, including decollement

% Create evenly spaced grid on which to interpolate the scattered data into a surface
grid_res = 70;
xlin = linspace(min(data_horizons(:,1)),max(data_horizons(:,1)),grid_res);
ylin = linspace(min(data_horizons(:,2)),max(data_horizons(:,2)),grid_res);
zlin = linspace(0,Zmax,grid_res);

[X2d,Y2d] = meshgrid(xlin,ylin);
[X3d,Y3d,Z3d] = meshgrid(xlin,ylin,zlin);

% Set different xlimits for different sides of the plate boundary
% For AU plate: minX:trench, based on seafloor depth
grid_res = 100;
xlin_au = linspace(min(data_horizons(:,1)), data_horizons((data_horizons(:,3)==max(data_horizons(layerID==1,3))),1), grid_res);
ylin_au = linspace(min(data_horizons(:,2)), max(data_horizons(:,2)), grid_res);
[Xau,Yau] = meshgrid(xlin_au,ylin_au);

% For PA plate: trench:maxX
xlin_pa = linspace(data_horizons((data_horizons(:,3)==max(data_horizons(layerID==1,3))),1), max(data_horizons(:,1)), grid_res);
ylin_pa = linspace(min(data_horizons(:,2)), max(data_horizons(:,2)), grid_res);
[Xpa,Ypa] = meshgrid(xlin_pa,ylin_pa);

% For AU moho: minX:end_of_the_slab
xlin_sl = linspace(min(data_horizons(:,1)), max(data_horizons((data_horizons(:,3)==max(data_horizons(layerID==10,3))),1)), grid_res);
ylin_sl = linspace(min(data_horizons(:,2)), max(data_horizons(:,2)), grid_res);
[Xsl,Ysl] = meshgrid(xlin_sl,ylin_sl);

% Interpolate surfaces from the scattered data for each interface
f1 = scatteredInterpolant(h1(:,1),h1(:,2),h1(:,3),'linear','nearest'); % seafloor surface
Zb1 = f1(X2d,Y2d);
f2 = scatteredInterpolant(h2(:,1),h2(:,2),h2(:,3),'linear','nearest'); % AU basement surface
Zb2 = f2(Xau,Yau);
f3 = scatteredInterpolant(h3(:,1),h3(:,2),h3(:,3),'linear','nearest'); % PA basement surface
Zb3 = f3(Xpa,Ypa);
f4 = scatteredInterpolant(h4(:,1),h4(:,2),h4(:,3),'linear','nearest'); % AU moho surface
Zb4 = f4(Xsl,Ysl);
f5 = scatteredInterpolant(h5(:,1),h5(:,2),h5(:,3),'linear','nearest'); % PA moho surface
Zb5 = f5(Xpa,Ypa);

% Evaluate the surface functions at the x,y coordinates of the prism centroids
layer1surf = interp2(X2d,Y2d,Zb1,prism_centroids(:,1),prism_centroids(:,2)); % seafloor at prism centroids
layer2surf = interp2(Xau,Yau,Zb2,prism_centroids(:,1),prism_centroids(:,2)); % AU basement at prism centroids
layer3surf = interp2(Xpa,Ypa,Zb3,prism_centroids(:,1),prism_centroids(:,2)); % PA basement at prism centroids
layer4surf = interp2(Xsl,Ysl,Zb4,prism_centroids(:,1),prism_centroids(:,2)); % AU moho at prism centroids
layer5surf = interp2(Xpa,Ypa,Zb5,prism_centroids(:,1),prism_centroids(:,2)); % PA moho at prism centroids

% Obtain indices of prisms lying above and below specified surface to define layers with logical vectors
% Include transitional layers between seafloor and rock (~1 prism thick) to prevent jumps in stdev at the seafloor
layer1 = prism_centroids(:,3) < layer1surf; % ocean layer
layer1_trans1 = (prism_centroids(:,3) > layer1surf) & (prism_centroids(:,3) < layer1surf+(Zmax/length(z_c)));
layer1_trans2 = (prism_centroids(:,3) > layer1surf+(Zmax/length(z_c))) & (prism_centroids(:,3) < layer1surf+(2*Zmax/length(z_c)));
layer1_trans3 = (prism_centroids(:,3) > layer1surf+(2*Zmax/length(z_c))) & (prism_centroids(:,3) < layer1surf+(3*Zmax/length(z_c)));

layer2 = (prism_centroids(:,3) > layer1surf) & (prism_centroids(:,3) < layer2surf); % AU plate sediments
layer3 = (prism_centroids(:,3) > layer1surf) & (prism_centroids(:,3) < layer3surf); % PA plate sediment
layer4 = (prism_centroids(:,3) > layer2surf) & (prism_centroids(:,3) < layer4surf); % AU plate basement rock
layer5 = (prism_centroids(:,3) > layer3surf) & (prism_centroids(:,3) < layer5surf); % PA plate basement rock
layer6 = (prism_centroids(:,3) < layer4surf) & (prism_centroids(:,3) > layer5surf); % slab
layer7 = (prism_centroids(:,3) < layer5surf) & (prism_centroids(:,3) > layer1surf) & ...
         (prism_centroids(:,1) > 1.589e5)    & (prism_centroids(:,1) < 1.85e5) & ...
         (prism_centroids(:,2) > 2.15e5)     & (prism_centroids(:,2) < 3.062e5);    % Western half of Puysegur Ridge
layer8 = (prism_centroids(:,3) < layer5surf) & (prism_centroids(:,3) > layer1surf) & ...
         (prism_centroids(:,1) > 1.85e5)     & (prism_centroids(:,1) < 2.185e5) & ...
         (prism_centroids(:,2) > 2.15e5)     & (prism_centroids(:,2) < 3.062e5);    % Eastern half of Pusegur Ridge
layer9 = (prism_centroids(:,3) > layer5surf) & (prism_centroids(:,3) > 2e5);        % extra mantle constraint

layers = {layer1; layer2; layer3; layer4; layer5; layer6; layer7; layer8; layer9};

% ======================================================================= %
% 3D interpolate density (from velocity modeled density) and map to prisms
% ========================================================================
% 3D interpolate the prior density values and errors onto the 3D grid created above 
frho = scatteredInterpolant(data_density(:,1),data_density(:,2),data_density(:,3),data_density(:,4),'linear','nearest');
density = frho(X3d,Y3d,Z3d);

data_density_stdev = data_density;
data_density_stdev(:,4) = [];
data_density_stdev(any(isnan(data_density_stdev),2),:) = [];

frho_stdev = scatteredInterpolant(data_density_stdev(:,1),data_density_stdev(:,2),data_density_stdev(:,3),data_density_stdev(:,4),'linear','nearest');
density_stdev = frho_stdev(X3d,Y3d,Z3d);

% Interpolate the prior values of the 3D grid at the locations of the prism centroids
prism_density = interp3(X3d,Y3d,Z3d,density,prism_centroids(:,1),prism_centroids(:,2),prism_centroids(:,3));
prism_density_stdev = interp3(X3d,Y3d,Z3d,density_stdev,prism_centroids(:,1),prism_centroids(:,2),prism_centroids(:,3));
prism_density_stdev(isnan(prism_density_stdev)) = max(prism_density_stdev);

% ======================================================================= %
% Construct the prior vector using densities above 
% ================================================
% Set prior densities of different model regions defined by the layer vectors above
mu_pr = prism_density;
mu_pr(layer1) = 1027; % fix the density of ocean water
mu_pr(layer4) = 2900; % fix the density of the oceanic plate
mu_pr(mu_pr > 3300) = 3300; % fix the maximum density of the mantle

% Force boundary prism priors to be the same density as their adjacent prisms in the main model domain
mu_pr(prism_centroids(:,2) == min(y_c)) = mu_pr(prism_centroids(:,2) == y_c(2));
mu_pr(prism_centroids(:,2) == max(y_c)) = mu_pr(prism_centroids(:,2) == y_c(end-1));
mu_pr(prism_centroids(:,1) == min(x_c)) = mu_pr(prism_centroids(:,1) == x_c(2));
mu_pr(prism_centroids(:,1) == max(x_c)) = mu_pr(prism_centroids(:,1) == x_c(end-1));

% Remove the mean of the prior densities for each layer to create relative drho
layer_mean = zeros(length(mu_pr),1);
for i = 1:size(prism_centroids,1)
    layer_mu_pr = mu_pr(i:length(z_c):end);
    layer_mean(i:length(z_c):end) = mean(layer_mu_pr);
end

mu_abs = mu_pr;
mu_pr = mu_abs - layer_mean;

% ======================================================================= %
% Construct prior covariance based on nearest-neighbor data
% =========================================================
% Perform knnsearch between prism_centroids and seismic horizon data
[id,D] = knnsearch(data_structure(:,1:3),prism_centroids,'K',1);

% Define the sigma0 values at the prism_centroids
sigma0_prisms = prism_density_stdev;

% Calculate sigma as a function of distance away from the lines given a maximum far-field sigma value
sigma_max = 500;
minD = 100; % distance below which all sigma values take on sigma0
k = (-1/minD)*log(1-(1/(sigma_max-max(sigma0_prisms))));
sigma = (sigma_max-max(sigma0_prisms)).*(1-exp(-k.*D))+sigma0_prisms;

% Reassign sigma values to the ocean and transition layers and boundary prisms 
sigma_raw = sigma;
sigma(layer1) = 5;
sigma(layer1_trans1) = 20;
sigma(layer1_trans2) = 50;
sigma(layer1_trans3) = 100;

sigma(prism_centroids(:,2) == min(y_c)) = sigma(prism_centroids(:,2) == y_c(2));
sigma(prism_centroids(:,2) == max(y_c)) = sigma(prism_centroids(:,2) == y_c(end-1));
sigma(prism_centroids(:,1) == min(x_c)) = sigma(prism_centroids(:,1) == x_c(2));
sigma(prism_centroids(:,1) == max(x_c)) = sigma(prism_centroids(:,1) == x_c(end-1));

sigma_pr = sigma;

% Construct the prior covariance operator
Cm = diag(sigma_pr.^2);
Cm = sparse(Cm);

% Save the prior density vector and covariance matrix to load to inversion
save([filepath modelname '_Prior_PuysegurGJI.mat'],'mu_abs','mu_pr','sigma_pr','Cm','layer_mean','layers','D','prism_density_stdev')

% ======================================================================= %
% Plot prior densities and their standard deviations
% ==================================================
Xgrid = repmat(x_c,length(y_c),1);
Xgrid = repmat(Xgrid,[1 1 length(z_c)]);
Ygrid = repmat(y_c',1,length(x_c));
Ygrid = repmat(Ygrid,[1 1 length(z_c)]);

[~,~,grid_dist] = vec2grid(3,D,x_c,y_c,z_c,[]);
[~,~,grid_sigma] = vec2grid(3,sigma,x_c,y_c,z_c,[]);
[~,~,grid_density] = vec2grid(3,prism_density,x_c,y_c,z_c,[]);
[~,~,grid_density_stdev] = vec2grid(3,prism_density_stdev,x_c,y_c,z_c,[]);

% Contour plot the nearest neighbor distances at various depth slices
zlevels = [1 12 18];
fig1 = contour_priors(zlevels,z_c,grid_dist,Xgrid,Ygrid,'Nearest Neighbor Distance (m)');

% Contour plot the densities at various depths
zlevels = [1 7 12 18];
fig2 = contour_priors(zlevels,z_c,grid_density,Xgrid,Ygrid,'Prism Density (kg/m^3)');

% Contour plot the raw density stdev values at various depths
zlevels = [1 7 12 18];
fig3 = contour_priors(zlevels,z_c,grid_density_stdev,Xgrid,Ygrid,'Prism Density Stdev (kg/m^3)');

% Contour plot the distance weighted stdev at various depth slices: sigma
zlevels = [1 7 12 18];
fig4 = contour_priors(zlevels,z_c,grid_sigma,Xgrid,Ygrid,'Prior StDev');

% Plot prior standard deviation as a function of distance from seismic lines
fig5 = figure;
plot(D./1000,sigma,'bo')
xlim([0 350])
xlabel('Nearest Neighbor Distance (km)')
ylabel('Prior Standard Deviation (kg/m^3)')

% ======================================================================= %
% Plot interpolated surfaces, seismic lines, and prism_centroids
% ==============================================================
olivine_green = [0.4660 0.6740 0.1880];
royal_purple = [0.4940 0.1840 0.5560];

fig6 = figure;

% plot bathymetry as a colored shaded surface
surf(X2d,Y2d,Zb1); hold on
cb = colorbar('eastoutside');
caxis([min(h1(:,3)) max(h1(:,3))])
colormap(flipud(parula))
ylabel(cb,'Bathymetric Depth')

% plot data points and prism centroids for the seafloor
scatter3(prism_centroids(layer1,1),prism_centroids(layer1,2),prism_centroids(layer1,3),'b','filled')
hold on
scatter3(h1(:,1),h1(:,2),h1(:,3),'c','filled')
hold on

% plot surface for the basement and prism centroids for the sediments
surf(Xau,Yau,Zb2,'FaceColor',olivine_green); hold on
scatter3(prism_centroids(layer2,1),prism_centroids(layer2,2),prism_centroids(layer2,3),'m','filled')
hold on
surf(Xpa,Ypa,Zb3,'FaceColor',olivine_green); hold on
scatter3(prism_centroids(layer3,1),prism_centroids(layer3,2),prism_centroids(layer3,3),'m','filled')
hold on

% plot surface for the subducting slab and prism centroids for the slab
surf(Xsl,Ysl,Zb4,'FaceColor',royal_purple); hold on
scatter3(prism_centroids(layer4,1),prism_centroids(layer4,2),prism_centroids(layer4,3),'g','filled')
hold on
surf(Xpa,Ypa,Zb5,'FaceColor',royal_purple); hold on
scatter3(prism_centroids(layer5,1),prism_centroids(layer5,2),prism_centroids(layer5,3),'c','filled')
hold on

% plot prism centroids for the deeper portion of the slab
scatter3(prism_centroids(layer6,1),prism_centroids(layer6,2),prism_centroids(layer6,3),'k','filled')
hold on

% plot prism centroids for the west and east halves of Puysegur Ridge
scatter3(prism_centroids(layer7,1),prism_centroids(layer7,2),prism_centroids(layer7,3),'r','filled')
hold on
scatter3(prism_centroids(layer8,1),prism_centroids(layer8,2),prism_centroids(layer8,3),'r','filled')
hold on

% plot the basement horizon along the seismic lines
plot3(h3(:,1),h3(:,2),h3(:,3),'r*','MarkerSize',15)

set(gca,'ZDir','reverse')
xlabel('X (m)')
ylabel('Y (m)')
zlabel('Depth (m)')

% ======================================================================= %
% Plot slices through 3D density distribution of the prior
% ========================================================
% Grid prior into 3D space for plotting
[~,~,mu_pr3D] = vec2grid(3,mu_pr,x_c,y_c,z_c,[]);
[~,~,mu_abs3D] = vec2grid(3,mu_abs,x_c,y_c,z_c,[]);

xslice = [100e3];
yslice = [50e3];
zslice = [5e3 12e3];

viewangle = [-35 30];

% Plot 3D absolute density of the prior
fig7 = figure;
plot_density_model(x,y,z,xslice,yslice,zslice,[],mu_abs3D)
set(gca,'ZDir','reverse')
caxis([min(mu_abs3D(:)) max(mu_abs3D(:))])
cb = colorbar;
ylabel(cb,'Density (kg/m^3)')
xlabel('X')
ylabel('Y')
zlabel('Depth')
xlim([x(2) x(end-1)])
ylim([y(2) y(end-1)])
view(viewangle)

% Plot 3D relative density of the prior
fig8 = figure;
plot_density_model(x,y,z,xslice,yslice,zslice,[],mu_pr3D)
set(gca,'ZDir','reverse')
caxis([min(mu_pr3D(:)) max(mu_pr3D(:))])
cb = colorbar;
ylabel(cb,'Differential Density (kg/m^3)')
xlabel('X')
ylabel('Y')
zlabel('Depth')
xlim([x(2) x(end-1)])
ylim([y(2) y(end-1)])
view(viewangle)


function [fig] = contour_priors(layers,zrange,var,X,Y,label)
    fig = figure;
    for i = 1:length(layers)
        [~,cf] = contourf(X(:,:,layers(i)),Y(:,:,layers(i)),var(:,:,layers(i)),40);
        cf.ContourZLevel = zrange(layers(i));
        hold on
    end
    view(3)
    set(gca,'ZDir','reverse')
    cb = colorbar;
    ylabel(cb,label)
    xlabel('X (m)')
    ylabel('Y (m)')
    zlabel('Depth (m)')
    title(label)
end




