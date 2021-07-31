%                       BayGrav3D, Version 1.0 
%       © 2021 Erin Hightower, California Institute of Technology
%   ----------------------------------------------------------------------
%   Supported in part by the US National Science Foundation
%   (http://www.nsf.gov).
%
%   Free for non-commercial academic use only.
%
%   ----------------------------------------------------------------------
%   This file [construct_Gmatrix.m] is part of BayGrav3D.

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
% construct_Gmatrix.m
% -------------------
% Computes the G matrix, using the function gravprism.m, that is used in 
% the gravity inversion. Also, computes the gravitational attraction of the
% top 1 m of ocean water to use as a correction to the gravity. 
%
% INPUT: 
%   x:         x-coordinates of the mesh nodes, a vector
%   y:         y-coordinates of the mesh nodes, a vector
%   z:         z-coordinates of the mesh nodes, a vector
%   x_c:       x-coordinates of the prism centroids, a vector
%   y_c:       y-coordinates of the prism centroids, a vector
%   z_c:       z-coordinates of the prism centroids, a vector
%   data:      data = [x,y,gravity]
%              Data array with which to calculate the prism to observation
%              point geometry.
%
% OUTPUT:
%   G:         G matrix used in the inversion
%   grav_corr: gravity for the top 1 m of ocean used as a correction
% ======================================================================= %

function [G,grav_corr] = construct_Gmatrix(x,y,z,x_c,y_c,z_c,data)    
    num_param = length(x_c)*length(y_c)*length(z_c);
    
   % Replicate prism edge coordinates to create a vector with edge 
   % coordinate values for each prism.
   % Prisms are ordered in the direction of increasing z, then y, then x
    x1 = x(1:end-1);
    x1 = repmat(x1,length(z_c)*length(y_c),1);
    x1 = reshape(x1,[num_param 1]);
    
    x2 = x(2:end);
    x2 = repmat(x2,length(z_c)*length(y_c),1);
    x2 = reshape(x2,[num_param 1]);
    
    y1 = y(1:end-1);
    y1 = repmat(y1,length(z_c),1);
    y1 = reshape(y1,[length(z_c)*length(y_c) 1]);
    y1 = repmat(y1,1,length(x_c));
    y1 = reshape(y1,[num_param 1]);
    
    y2 = y(2:end);
    y2 = repmat(y2,length(z_c),1);
    y2 = reshape(y2,[length(z_c)*length(y_c) 1]);
    y2 = repmat(y2,1,length(x_c));
    y2 = reshape(y2,[num_param 1]);
    
    z1 = z(1:end-1);
    z1 = repmat(z1,length(y_c)*length(x_c),1);
    z1 = z1';
    z1 = reshape(z1,[num_param 1]);
    
    z2 = z(2:end);
    z2 = repmat(z2,length(y_c)*length(x_c),1);
    z2 = z2';
    z2 = reshape(z2,[num_param 1]);
    
    edge1 = [x1 y1 z1];
    edge2 = [x2 y2 z2];
    
    % Calculate prism edge to observation point geometry per the equation 
    % in Turcotte and Schubert, 2014, Geodynamics, 3rd Edition, p. 527.
    % Uses function gravprism.m and produces G matrix used for inversion
    geom = zeros(num_param,size(data,1));
    for i = 1:num_param
        dx1 = edge1(i,1) - data(:,1);
        dx2 = edge2(i,1) - data(:,1);
        dy1 = edge1(i,2) - data(:,2);
        dy2 = edge2(i,2) - data(:,2);
        dz1 = edge1(i,3) - zeros(size(data,1),1);
        dz2 = edge2(i,3) - zeros(size(data,1),1); 
        
        [geometry] = gravprism(dx1,dx2,dy1,dy2,dz1,dz2);
        geom(i,:) = geometry;
    end

   % Determine the geometric constant of the prism that is the top 1 m of
   % ocean water for the correction to the gravity data
    dx1_corr = x(1)-data(:,1);
    dx2_corr = x(end)-data(:,1);
    dy1_corr = y(1)-data(:,2);
    dy2_corr = y(end)-data(:,2);
    dz1_corr = 0 - zeros(size(data,1),1);
    dz2_corr = z(1) - zeros(size(data,1),1);
    
   % Calculate the G matrix used in the gravity inversion  
    const = (6.6732e-11)*1e5;
    G = const*geom';
    
    rho_seawater = 1027;
    geom_corr = gravprism(dx1_corr,dx2_corr,dy1_corr,dy2_corr,dz1_corr,dz2_corr);
    grav_corr = const.*geom_corr.*rho_seawater;
    
end