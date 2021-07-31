%                       BayGrav3D, Version 1.0 
%       © 2021 Erin Hightower, California Institute of Technology
%   ----------------------------------------------------------------------
%   Supported in part by the US National Science Foundation
%   (http://www.nsf.gov).
%
%   Free for non-commercial academic use only.
%
%   ----------------------------------------------------------------------
%   This file [plot_Gravity.m] is part of BayGrav3D.

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
% plot_Gravity.m
% --------------
% This script plots the observed, predicted, and residual gravity
% associated with a given inversion.
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
%   plots of observed, predicted, and residual gravity 
% ======================================================================= %

filepath = 'BayGrav3D/Example1/';
modelname = 'Example1';

input_data = [filepath modelname '_InversionResults.mat'];
input_mesh = [filepath modelname '_mesh.mat'];
input_model = [filepath modelname '_PosteriorModel.txt'];
drop_boundary = 'Y';

% ======================================================================= %
% Load the data
% =============
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

% ======================================================================= %
% Plot prior density distribution and forward gravity
% ===================================================
contours = 30;
contour_color = [50/255 50/255 50/255];
contour_linewidth = 0.05;
colorlims = [min(grav_obs) max(grav_obs)];

fig = figure;

% Plot the observed gravity contour map
ax1 = subplot(3,1,1);
contourf(X,Y,grav_obs2D,contours,'LineColor',contour_color,'LineWidth',contour_linewidth)
hold on
caxis(colorlims)
ylabel('Y')
xticklabels([])

% Plot predicted gravity contour map
ax2 = subplot(3,1,2);
contourf(X,Y,grav_pred2D,contours,'LineColor',contour_color,'LineWidth',contour_linewidth)
hold on
caxis(colorlims)
ylabel('Y')
xticklabels([])

% Plot residual between the observed and predicted gravity
ax3 = subplot(3,1,3);
contourf(X,Y,grav_error2D,contours,'LineColor',contour_color,'LineWidth',contour_linewidth)
caxis(colorlims)
xlabel('X')
ylabel('Y')

ax1.OuterPosition(3) = 0.5;
ax2.OuterPosition(3) = 0.5;
ax3.OuterPosition(3) = 0.5;

cb = colorbar;
cb.Position = [ax1.OuterPosition(1)+ax1.OuterPosition(3)-0.03 ax1.Position(2) cb.Position(3) ax1.Position(4)];
ylabel(cb,'Gravity (mGal)')

title(ax1,'Observed Gravity')
title(ax2,'Predicted Gravity')
title(ax3,'Residual Gravity')









