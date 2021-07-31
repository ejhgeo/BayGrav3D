%                       BayGrav3D, Version 1.0 
%       © 2021 Erin Hightower, California Institute of Technology
%   ----------------------------------------------------------------------
%   Supported in part by the US National Science Foundation
%   (http://www.nsf.gov).
%
%   Free for non-commercial academic use only.
%
%   ----------------------------------------------------------------------
%   This file [plot_GravityAndDensityProfiles.m] is part of BayGrav3D.

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
% plot_GravityAndDensityProfiles.m
% --------------------------------
% This script plots profiles of the observed, predicted, and prior gravity
% associated with a given inversion above the corresponding slice through
% the 3D density domain of the predicted density solved for in the
% inversion, the prior density used to constrain the inversion, and the
% posterior standard deviation of the estimated density. 
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
%   plots of gravity profiles over their corresponding density profiles 
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
grav_obs = grav_data(:,3);
grid_res = 300;
Xdata = grav_data(:,1);
Ydata = grav_data(:,2);
grav_error = abs(grav_obs - grav_pred);

[X,Y,grav_obs2D]  = vec2grid(2,grav_obs,Xdata,Ydata,[],grid_res);
[~,~,grav_prior2D] = vec2grid(2,grav_prior,Xdata,Ydata,[],grid_res);
[~,~,grav_pred2D]  = vec2grid(2,grav_pred,Xdata,Ydata,[],grid_res);
[~,~,grav_error2D] = vec2grid(2,grav_error,Xdata,Ydata,[],grid_res);

% Reshape prior density and standard deviation values into 3D grid 
[~,~,mu_pr3D] = vec2grid(3,mu_pr,x_c,y_c,z_c,[]);
[~,~,mu_abs3D] = vec2grid(3,mu_abs,x_c,y_c,z_c,[]);

% Reshape predicted density and standard deviation values into 3D grid
[~,~,m3D] = vec2grid(3,m,x_c,y_c,z_c,[]);
[~,~,m_abs3D] = vec2grid(3,m_abs,x_c,y_c,z_c,[]);
[~,~,post_stdev3D] = vec2grid(3,post_stdev,x_c,y_c,z_c,[]);

% ======================================================================= %
% Create gravity profiles
% =======================
slice_orientation = 'x';
slice_position = 150e3;

% Extract gravity profile at the slice location using interp2
[profile_axis,grav_obs_profile] = data_profile(X,Y,grav_obs2D,slice_position,slice_orientation);
[~,grav_pred_profile] = data_profile(X,Y,grav_pred2D,slice_position,slice_orientation);
[~,grav_prior_profile] = data_profile(X,Y,grav_prior2D,slice_position,slice_orientation);

% ======================================================================= %
% Plot the gravity transect over the density profile
% ==================================================
clims_grav = [min(grav_obs_profile) max(grav_obs_profile)];
clims_density = [min(mu_abs) max(mu_abs)];

xlims = [0 Xmax];
ylims = [0 Ymax];
zlims = [0 Zmax];

fig = figure;

% Plot gravity profiles for observed, predicted and prior gravity
ax1 = subplot(4,1,1);
p1 = plot(profile_axis,grav_obs_profile,'k','LineWidth',1.5);
hold on
p2 = plot(profile_axis,grav_pred_profile,'--k','LineWidth',1.5);
hold on
p3 = plot(profile_axis,grav_prior_profile,'-.k','LineWidth',1.5);
xlim([0 max(profile_axis)])
ylim(clims_grav)
xticklabels([])
ylabel('Gravity (mGal)')

% Plot the posterior density distribution profile from the inversion
ax2 = subplot(4,1,2);
[view_dir] = plot_density_model(x,y,z,slice_position,slice_position,slice_position,slice_orientation,m_abs3D);
set(gca,'ZDir','reverse')
caxis(clims_density)
view(view_dir)
rotate3d off
xlim(xlims)
ylim(ylims)
zlim(zlims)
xticklabels([])
yticklabels([])
zlabel('Depth')

% Plot the prior density distribution profile
ax3 = subplot(4,1,3);
[view_dir] = plot_density_model(x,y,z,slice_position,slice_position,slice_position,slice_orientation,mu_abs3D);
set(gca,'ZDir','reverse')
caxis(clims_density)
view(view_dir)
rotate3d off
xlim(xlims)
ylim(ylims)
zlim(zlims)
xticklabels([])
yticklabels([])
zlabel('Depth')

% Plot posterior standard deviation of density from inversion 
ax4 = subplot(4,1,4);
[view_dir] = plot_density_model(x,y,z,slice_position,slice_position,slice_position,slice_orientation,post_stdev3D);
set(gca,'ZDir','reverse')
colormap(ax4,flipud(pink))
caxis([min(post_stdev3D(:)) max(post_stdev3D(:))])
view(view_dir)
rotate3d off
xlim(xlims)
ylim(ylims)
zlim(zlims)
xlabel('Distance')
ylabel('Distance')
zlabel('Depth')

% Plot colorbars for density profiles and legend for gravity profiles
ax1.OuterPosition(3) = 0.9;
ax2.OuterPosition(3) = 0.9;
ax3.OuterPosition(3) = 0.9;
ax4.OuterPosition(3) = 0.9;

cb2 = colorbar(ax2);
cb2.Position = [ax2.OuterPosition(1)+ax2.OuterPosition(3)-0.07 ax3.OuterPosition(2) cb2.Position(3) 2*ax2.OuterPosition(4)];
ylabel(cb2,'Density (kg/m^3)')

cb4 = colorbar(ax4);
cb4.Position = [ax4.OuterPosition(1)+ax4.OuterPosition(3)-0.07 ax4.Position(2) cb4.Position(3) ax4.Position(4)];
ylabel(cb4,'Posterior \sigma (kg/m^3)')

lg = legend(ax1,[p1 p2 p3],{'Obs. Gravity','Predicted Gravity','Prior Gravity'},'FontSize',6);
lg.Position = [ax1.Position(1)+0.015 ax1.Position(2)+0.015 lg.Position(3) lg.Position(4)];








