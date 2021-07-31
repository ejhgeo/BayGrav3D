%                       BayGrav3D, Version 1.0 
%       © 2021 Erin Hightower, California Institute of Technology
%   ----------------------------------------------------------------------
%   Supported in part by the US National Science Foundation
%   (http://www.nsf.gov).
%
%   Free for non-commercial academic use only.
%
%   ----------------------------------------------------------------------
%   This file [plot_Prior.m] is part of BayGrav3D.

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
% plot_Prior.m
% --------------
% This script plots the prior gravity and density, as well as the prior
% variance in 3D associated with a given inversion.
%
% *** This is an example: by all means modify as needed or write your own 
%     script to make your plots from the inversion output.
%
% INPUT:
%   filepath:    path to directory where this a given model is stored
%   modelname:   name of the model/inversion to be plotted
%   input_data:  .mat file of results produced by grav_inversion.m
%   input_mesh:  .mat file of mesh parameters used in the inversion
%   input_model: .txt file of posterior model produced by grav_inversion.m
%   drop_boundary: 'Y' to drop the boundary prisms 
%                  'N' to leave the boundary prisms in the plots
%
% OUTPUT: 
%   plots of prior variance in 3D and prior gravity and density
% ======================================================================= %

filepath = 'BayGrav3D/Example1/';
modelname = 'Example1';

input_data = [filepath modelname '_InversionResults.mat'];
input_mesh = [filepath modelname '_mesh.mat'];
input_model = [filepath modelname '_PosteriorModel.txt'];
drop_boundary = 'Y';

load(input_data)
load(input_mesh)
model = readtable(input_model,'ReadVariableNames',1);
model = table2array(model);

if strcmp(drop_boundary,'Y') == 1
    [model,layers,layer_mean,prism_centroids,m,m_abs,mu_pr,mu_abs,sigma_pr,post_stdev,post_res,x,y,x_c,y_c] ...
        = drop_boundary_prisms(model,layers,layer_mean,x,x_c,y,y_c);
end

Xmin = min(x);
Xmax = max(x);
Ymin = min(y);
Ymax = max(y);

% ======================================================================= %
% Grid the results for plotting
% =============================
% Grid coordinates into 3D space to plot in 3D
[~,~,X3D] = vec2grid(3,prism_centroids(:,1),x_c,y_c,z_c,[]);
[~,~,Y3D] = vec2grid(3,prism_centroids(:,2),x_c,y_c,z_c,[]);
[~,~,Z3D] = vec2grid(3,prism_centroids(:,3),x_c,y_c,z_c,[]);

% Reshape gravity into 2D grid for plotting
grid_res = 300;
Xdata = grav_data(:,1);
Ydata = grav_data(:,2);
[X,Y,grav_prior2D] = vec2grid(2,grav_prior,Xdata,Ydata,[],grid_res);

% Reshape prior density and standard deviation values into 3D grid 
[~,~,mu_pr3D] = vec2grid(3,mu_pr,x_c,y_c,z_c,[]);
[~,~,mu_abs3D] = vec2grid(3,mu_abs,x_c,y_c,z_c,[]);
[~,~,sigma_pr3D] = vec2grid(3,sigma_pr,x_c,y_c,z_c,[]);

% ======================================================================= %
% Plot prior density distribution and forward gravity
% ===================================================
xslice = [100e3];
yslice = [50e3];
zslice = [5e3 12e3];

contours = 30;
contour_color = [50/255 50/255 50/255];
contour_linewidth = 0.1;

viewangle = [-35 30];

fig1 = figure;

% Contour plot the gravity from the prior
ax1 = subplot(2,1,1);
contourf(X,Y,grav_prior2D,contours,'LineColor',contour_color,'LineWidth',contour_linewidth);
axis([Xmin Xmax Ymin Ymax 0 0.1])
set(gca,'YDir','normal')
ax1.ZAxis.Visible = 'off';
ax1.Color = 'none';
ax1.Box = 'off';
view(viewangle)

% Plot the 3D density distribution of the prior
ax2 = subplot(2,1,2);
plot_density_model(x,y,z,xslice,yslice,zslice,[],mu_abs3D)
set(gca,'ZDir','reverse')
caxis([min(mu_abs3D(:)) max(mu_abs3D(:))])
xlabel('X')
ylabel('Y')
zlabel('Depth')

cb1 = colorbar(ax1);
ylabel(cb1,'Gravity (mGal)')
cb2 = colorbar(ax2);
ylabel(cb2,'Density (kg/m^3)')


% ======================================================================= %
% Plot the Prior Standard Deviation in 3D
% =======================================
zlevels = [1 round(0.5*length(z_c)) round(0.75*length(z_c)) length(z_c)];
contours = 15;
contour_color = [50/255 50/255 50/255];
contour_linewidth = 0.1;

% Contour standard deviation in 3D at the different zlevels
fig2 = figure;
for i = 1:length(zlevels)
    [~,cf] = contourf(X3D(:,:,zlevels(i)),Y3D(:,:,zlevels(i)),sigma_pr3D(:,:,zlevels(i)),contours);
    set(cf,'LineColor',contour_color,'LineWidth',contour_linewidth)
    cf.ContourZLevel = z_c(zlevels(i));
    hold on
end
view([-35 15])
set(gca,'ZDir','reverse')
xlabel('X')
ylabel('Y')
zlabel('Depth')

cb = colorbar;
ylabel(cb,{'Prior \sigma','(kg/m^3)'})



