%                       BayGrav3D, Version 1.0 
%       © 2021 Erin Hightower, California Institute of Technology
%   ----------------------------------------------------------------------
%   Supported in part by the US National Science Foundation
%   (http://www.nsf.gov).
%
%   Free for non-commercial academic use only.
%
%   ----------------------------------------------------------------------
%   This file [gravprism.m] is part of BayGrav3D.

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
% gravprism.m
% -----------
% Computes the geometrical coefficient needed to describe the gravitational
% attraction of a buried prism relative to an observation point on the
% surface, using the equation from Turcotte & Schubert, 3rd Edition (2014)
%
% INPUT:
%   dx1:    [x_1 - x_p] 
%           [x-coord of prism's 1st edge] - [x-coord of observation point]
%   dx2:    [x_2 - x_p]
%           [x-coord of prism's 2nd edge] - [x-coord of observation point]
%
%   dy1:    [y_1 - y_p]
%           [y-coord of prism's 1st edge] - [y-coord of observation point]
%   dy2:    [y_2 - y_p]
%           [y-coord of prism's 2nd edge] - [y-coord of observation point]
%
%   dz1:    [z_1 - z_p]
%           [z-coord of prism's top edge] - [z-coord of obs. point = 0]
%   dz2:    [z_2 - z_p]
%           [z-coord of prism's bottom edge] - [z-coord of obs. point = 0]
% 
%   For a single prism:     1 x N vector (1 prism  x # obs. points)
%   For an array of prisms: M x N matrix (# prisms x # obs. points)
%
% OUTPUT:
%   G:      M x N coefficient matrix of prism to obs. point geometry
% ======================================================================= %

function [G] = gravprism(dx1,dx2,dy1,dy2,dz1,dz2)
R111 = sqrt((dx1.^2) + (dy1.^2) + (dz1.^2));
R112 = sqrt((dx1.^2) + (dy1.^2) + (dz2.^2));
R121 = sqrt((dx1.^2) + (dy2.^2) + (dz1.^2));
R122 = sqrt((dx1.^2) + (dy2.^2) + (dz2.^2));
R211 = sqrt((dx2.^2) + (dy1.^2) + (dz1.^2));
R212 = sqrt((dx2.^2) + (dy1.^2) + (dz2.^2));
R221 = sqrt((dx2.^2) + (dy2.^2) + (dz1.^2));
R222 = sqrt((dx2.^2) + (dy2.^2) + (dz2.^2));

g111 = -(dz1.*atan((dx1.*dy1)./(dz1.*R111)) - dx1.*log(R111+dy1) - dy1.*log(R111+dx1));
g112 = (dz2.*atan((dx1.*dy1)./(dz2.*R112)) - dx1.*log(R112+dy1) - dy1.*log(R112+dx1));
g121 = (dz1.*atan((dx1.*dy2)./(dz1.*R121)) - dx1.*log(R121+dy2) - dy2.*log(R121+dx1));
g122 = -(dz2.*atan((dx1.*dy2)./(dz2.*R122)) - dx1.*log(R122+dy2) - dy2.*log(R122+dx1));

g211 = (dz1.*atan((dx2.*dy1)./(dz1.*R211)) - dx2.*log(R211+dy1) - dy1.*log(R211+dx2));
g212 = -(dz2.*atan((dx2.*dy1)./(dz2.*R212)) - dx2.*log(R212+dy1) - dy1.*log(R212+dx2));
g221 = -(dz1.*atan((dx2.*dy2)./(dz1.*R221)) - dx2.*log(R221+dy2) - dy2.*log(R221+dx2));
g222 = (dz2.*atan((dx2.*dy2)./(dz2.*R222)) - dx2.*log(R222+dy2) - dy2.*log(R222+dx2));

G = (g111 + g112 + g121 + g122 + g211 + g212 + g221 + g222);
end