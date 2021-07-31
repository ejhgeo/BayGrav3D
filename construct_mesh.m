%                       BayGrav3D, Version 1.0 
%       © 2021 Erin Hightower, California Institute of Technology
%   ----------------------------------------------------------------------
%   Supported in part by the US National Science Foundation
%   (http://www.nsf.gov).
%
%   Free for non-commercial academic use only.
%
%   ----------------------------------------------------------------------
%   This file [construct_mesh.m] is part of BayGrav3D.

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
% construct_mesh.m
% ----------------
% Constructs the mesh for the 3D subsurface domain, pads the domain with 
% 'infinite edge' boundary prisms, and exports node and centroid coords.
%
% INPUT:
%   Xmin:   minimum X limit of the model domain 
%   Xmax:   maximum X limit of the model domain
%   Ymin:   minimum Y limit of the model domain
%   Ymax:   maximum Y limit of the model domain
%   Zmin:   minimum Z limit of the model domain
%   Zmax:   maximum Z limit of the model domain
%           X,Y,Z limits are usually taken from the extent of the prior
%           data or the bathymetry grid; should all be in meters.
%
%   style:  'var' OR 'equal'
%           'var':   creates a mesh with decreasing resolution (increasing
%                    prism height) in the vertical direction with depth
%           'equal': creates a mesh with equal sized prisms everywhere
%
%   np_x:   number of prisms in the x-direction
%   np_y:   number of prisms in the y-direction
%   np_z:   number of prisms in the z-direction
%           This is the total number of prisms EXCLUDING the boundary
%           prisms. The boundary prisms get added in this script when the
%           mesh is constructed. The input for number of prisms to the
%           grav_inversion.m script should be the total number of prisms in
%           each direction overall. That script takes care of exlcuding the
%           boundary prisms from the count for input to this function.
%
%   res_z:  [resZ_min resZ_max]
%           Min and max resolution in vertical direction (prism height (m))
%
% OUTPUT: 
%   x:         x-coordinates of the mesh nodes, a vector
%   y:         y-coordinates of the mesh nodes, a vector
%   z:         z-coordinates of the mesh nodes, a vector
%   x_c:       x-coordinates of the prism centroids, a vector
%   y_c:       y-coordinates of the prism centroids, a vector
%   z_c:       z-coordinates of the prism centroids, a vector
%   prism_centroids: centroid coordinate for each prism in the mesh, matrix
% ======================================================================= %

function [x,y,z,x_c,y_c,z_c,prism_centroids] = construct_mesh(Xmin,Xmax,Ymin,Ymax,Zmin,Zmax,style,np_x,np_y,np_z,res_z,filename)
   if strcmp(style,'var') == 1
    % Create a grid with equal discretization in the x,y-directions and
    % increasing height in the z-direction
    zrange = 0:(1/(np_z-1)):1;
    Rz_min = res_z(1);
    Rz_max = res_z(2);
    zrange = (Zmax-Zmin).*zrange + Zmin;

    dz_var = ((Rz_max-Rz_min)/(Zmax-Zmin))*zrange + Rz_min;
    scale_z = (Zmax-Zmin)/sum(dz_var);
    
    dx = ones(1,np_x)*(Xmax-Xmin)/(np_x);
    dy = ones(1,np_y)*(Ymax-Ymin)/(np_y);
    dz = scale_z*dz_var;
    
   elseif strcmp(style,'equal') == 1
    % Create a grid with equal discretization in all the x,y,z-directions
    dx = ones(1,np_x)*(Xmax-Xmin)/(np_x);
    dy = ones(1,np_y)*(Ymax-Ymin)/(np_y);
    dz = ones(1,np_z)*(Zmax-Zmin)/(np_z);
   end
   
   % Calculate prism edge coordinates in each direction 
   % Add large edge prisms to create "infinite edge" boundary condition
    bound = 1000000;
    x = [Xmin-bound Xmin (Xmin+cumsum(dx)) (Xmin+sum(dx)+bound)];
    y = [Ymin-bound Ymin (Ymin+cumsum(dy)) (Ymin+sum(dy)+bound)];
    z = [Zmin Zmin+cumsum(dz)];
    
   % Calculate prism centroid coordinates and total number of parameters
    x_c = zeros(1,length(x)-1);
        for i = 1:length(x)-1
            x_c(i) = (x(i)+x(i+1))/2;
        end 
    y_c = zeros(1,length(y)-1);
        for j = 1:length(y)-1
            y_c(j) = (y(j)+y(j+1))/2;
        end
    z_c = zeros(1,length(z)-1);
        for k = 1:length(z)-1
            z_c(k) = (z(k)+z(k+1))/2;
        end
        
    prism_centroids = allcomb(x_c,y_c,z_c);
    
    save([filename '_mesh.mat'], 'Xmin','Xmax','Ymin','Ymax','Zmin','Zmax','np_x','np_y','np_z','res_z',...
                                 'x','y','z','x_c','y_c','z_c','prism_centroids')
    fid = fopen([filename '_prism_centroids.txt'],'w');
    fprintf(fid,'%s %s %s\n','x','y','z');
    fprintf(fid,'%.6f %.6f %.6f\n',prism_centroids');
    fclose(fid);
end