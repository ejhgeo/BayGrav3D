%                       BayGrav3D, Version 1.0 
%       © 2021 Erin Hightower, California Institute of Technology
%   ----------------------------------------------------------------------
%   Supported in part by the US National Science Foundation
%   (http://www.nsf.gov).
%
%   Free for non-commercial academic use only.
%
%   ----------------------------------------------------------------------
%   This file [drop_boundary_prisms.m] is part of BayGrav3D.

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
% drop_boundary_prisms.m
% ----------------------
% Drops the boundary prisms from the x and y dimensions of the 3D
% subsurface mesh for plotting purposes. Should be called at the beginning
% of a plotting script if desired. 
%
% * This has a lot of input/output and it's a bit hard coded, but I only
% use it for plotting, so that you don't have to plot the boundary prisms.
% *** Modify as needed to suit your model's needs. 
%
% INPUT:
%   model:      posterior_model array saved from grav_inversion.m
%   layers:     the cell array of layers saved from the prior construction
%   layer_mean: mean that was substracted from each prism centroid, vector
%               of equivalent length to m, the model parameters
%   x:   nodal x coordinates of the mesh to modify
%   x_c: centroid x coordinates of the mesh to modify
%   y:   nodal y coordinates of the mesh to modify
%   y_c: centroid y coordinates of the mesh to modify
%
% OUTPUT:
%   Input variables to the function but with values on the boundary prisms 
%   dropped.
%   All mesh variables, as saved from grav_inversion.m, but with values on
%   the boundary prisms dropped.
% ======================================================================= %

function [model,layers,layer_mean,prism_centroids,m,m_abs,mu_pr,mu_abs,sigma_pr,post_stdev,post_res,x,y,x_c,y_c] ...
          = drop_boundary_prisms(model,layers,layer_mean,x,x_c,y,y_c)
    % Find the indices of the boundary prisms to drop
    ind1 = find(model(:,1) == max(model(:,1)));
    ind2 = find(model(:,1) == min(model(:,1)));
    ind3 = find(model(:,2) == max(model(:,2)));
    ind4 = find(model(:,2) == min(model(:,2)));
    ind = [ind1; ind2; ind3; ind4];

    % Drop boundary prisms from the combined posterior model data
    model(ind,:) = [];
    
    % Drop boundary prisms from layers and layer_mean
    layer_out = cell(length(layers),1);
    for i = 1:length(layers)
        layer = layers{i};
        layer(ind) = [];
        layer_out{i} = layer;
    end
    layers = layer_out;
    layer_mean(ind,:) = [];
    
    prism_centroids = model(:,1:3);
    m = model(:,4);
    m_abs = model(:,5);
    mu_pr = model(:,6);
    mu_abs = model(:,7);
    sigma_pr = model(:,8);
    post_stdev = model(:,9);
    post_res = model(:,10);
    
    x = x(2:end-1);
    y = y(2:end-1);
    x_c = x_c(2:end-1);
    y_c = y_c(2:end-1);
end
    
    
    