%                       BayGrav3D, Version 1.0 
%       © 2021 Erin Hightower, California Institute of Technology
%   ----------------------------------------------------------------------
%   Supported in part by the US National Science Foundation
%   (http://www.nsf.gov).
%
%   Free for non-commercial academic use only.
%
%   ----------------------------------------------------------------------
%   This file [plot_VarianceAndResolution.m] is part of BayGrav3D.

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
% plot_VarianceAndResolution.m
% --------------
% This script plots the posterior resolution and standard deviation in 3D
% space for the results of a given inversion. 
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
%   plots of 3D resolution and standard deviation 
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

% Grid posterior standard deviation and resolution into 3D space
[~,~,post_stdev3D] = vec2grid(3,post_stdev,x_c,y_c,z_c,[]);
[~,~,post_res3D] = vec2grid(3,post_res,x_c,y_c,z_c,[]);

% ======================================================================= %
% Plot the posterior standard deviation and resolution in 3D
% ==========================================================
zslice = [0.5 9 17 24.9].*1e3;
view_angle = [-35 15];

xlims = [0 Xmax];
ylims = [0 Ymax];
zlims = [0 Zmax];

fig = figure;

% Plot 3D resolution
ax1 = subplot(2,1,1); 
plot_density_model(x,y,z,[],[],zslice,[],post_res3D)
view(view_angle)
set(gca,'ZDir','reverse')
caxis([min(post_res) 1])
xlim(xlims)
ylim(ylims)
zlim(zlims)
xlabel('X')
ylabel('Y')
zlabel('Z')

% Plot 3D standard deviation
ax2 = subplot(2,1,2);
plot_density_model(x,y,z,[],[],zslice,[],post_stdev3D)
view(view_angle)
set(gca,'ZDir','reverse')
caxis([min(post_stdev) max(post_stdev)])
xlim(xlims)
ylim(ylims)
zlim(zlims)
xlabel('X')
ylabel('Y')
zlabel('Depth')

% Plot colorbars for standard deviation and resolution
ax1.OuterPosition(3) = 0.5;
ax2.OuterPosition(3) = 0.5;

cb1 = colorbar(ax1);
cb1.Position = [ax1.OuterPosition(1)+ax1.OuterPosition(3)-0.03 ax1.Position(2) cb1.Position(3) ax1.Position(4)];
ylabel(cb1,'Resolution')

cb2 = colorbar(ax2);
cb2.Position = [ax2.OuterPosition(1)+ax2.OuterPosition(3)-0.03 ax2.Position(2) cb2.Position(3) ax2.Position(4)];
ylabel(cb2,'Posterior \sigma (kg/m^3)')











