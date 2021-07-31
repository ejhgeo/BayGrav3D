%                       BayGrav3D, Version 1.0 
%       © 2021 Erin Hightower, California Institute of Technology
%   ----------------------------------------------------------------------
%   Supported in part by the US National Science Foundation
%   (http://www.nsf.gov).
%
%   Free for non-commercial academic use only.
%
%   ----------------------------------------------------------------------
%   This file [plot_density_model.m] is part of BayGrav3D.

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
% plot_density_model.m
% --------------------
% Plots either a slice through the density model in 2D or multiple slices
% in 3D depending on the input options. 
%
% Plots prisms with x,y,z not x_c, y_c, z_c, otherwise matlab draws and 
% colors in the prisms starting with the centroid position, not the edge of
% the prism. To make them the same length, you have to only plot 1:end-1, 
% then plot the end prism, using end-1:end separately, so there are many
% combinations of ranges of x,y,z that have to be plotted
%
% To plot a fully 3D slice plot with multiple slices, set
% profile_orientation = [], in which case the function will not zero out 
% other directions or set a view orientation.
%
% INPUT:
%   xpl: x coordinates of the mesh (e.g. x, not x_c)
%   ypl: y coordinates of the mesh (e.g. y, not y_c)
%   zpl: z coordinates of the mesh (e.g. z, not z_c)
%   xsl: x values of N-S slices that intersect the x-axis, can be []
%   ysl: y values of E-W slices that intersect the y-axis, can be []
%   xsl: z values of horizontal slices that intersect the z-axis, can be []
%   profile_orientation: for 2D slice view only: 
%                        'x' for slices parallel to x-axis (intersecting y)
%                        'y' for slices parallel to y-axis (intersecting x)
%   rho: 3D gridded density model to slice through
%
% OUTPUT:
%   plots the density slice when called within a figure
%   [view_dir]: direction to look straight on at the slice in 2D, only
%               output if profile_orientation is set; pass to view() in
%               your figure
% ======================================================================= %

function [view_dir] = plot_density_model(xpl,ypl,zpl,xsl,ysl,zsl,profile_orientation,rho)
    if (strcmp(profile_orientation,'x') == 1)
        xsl = [];
        zsl = [];
        view_dir = [0 0];
    elseif (strcmp(profile_orientation,'y') == 1)
        ysl = [];
        zsl =[];
        view_dir = [90 0];
    end
    
    s1 = slice(xpl(1:end-1),ypl(1:end-1),zpl(1:end-1),rho,xsl,ysl,zsl);
    set(s1,'EdgeColor','none')
    hold on
    s2 = slice(xpl(1:end-1),ypl(1:end-1),zpl(end-1:end),rho(:,:,end-1:end),xsl,ysl,zsl);
    set(s2,'EdgeColor','none')
    hold on
    s3 = slice(xpl(1:end-1),ypl(end-1:end),zpl(1:end-1),rho(end-1:end,:,:),xsl,ysl,zsl);
    set(s3,'EdgeColor','none')
    hold on
    s4 = slice(xpl(1:end-1),ypl(end-1:end),zpl(end-1:end),rho(end-1:end,:,end-1:end),xsl,ysl,zsl);
    set(s4,'EdgeColor','none')
    hold on
    s5 = slice(xpl(end-1:end),ypl(1:end-1),zpl(1:end-1),rho(:,end-1:end,:),xsl,ysl,zsl);
    set(s5,'EdgeColor','none')
    hold on
    s6 = slice(xpl(end-1:end),ypl(end-1:end),zpl(1:end-1),rho(end-1:end,end-1:end,:),xsl,ysl,zsl);
    set(s6,'EdgeColor','none')
    hold on
    s7 = slice(xpl(end-1:end),ypl(1:end-1),zpl(end-1:end),rho(:,end-1:end,end-1:end),xsl,ysl,zsl);
    set(s7,'EdgeColor','none')
    hold on
    s8 = slice(xpl(end-1:end),ypl(end-1:end),zpl(end-1:end),rho(end-1:end,end-1:end,end-1:end),xsl,ysl,zsl);
    set(s8,'EdgeColor','none')
end