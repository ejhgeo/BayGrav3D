%                       BayGrav3D, Version 1.0 
%       © 2021 Erin Hightower, California Institute of Technology
%   ----------------------------------------------------------------------
%   Supported in part by the US National Science Foundation
%   (http://www.nsf.gov).
%
%   Free for non-commercial academic use only.
%
%   ----------------------------------------------------------------------
%   This file [data_profile.m] is part of BayGrav3D.

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
% data_profile.m
% --------------
% Samples a 2D grid (i.e. gravity) along a given transect and outputs a
% profile along that transect.
%
% INPUT:
%   X:    X coordinates of grid to sample, a 2D grid
%   Y:    Y coordinates of grid to sample, a 2D grid
%   grid: grid to sample (i.e. gravity), a 2D grid
%   slice_position: position (in distance) where the profile intersects the
%                   chosen axis, a scalar
%   slice_orientation: 'x' or 'y', orienation of the resulting profile, 'x' 
%                                  being parallel to the x axis (E-W
%                                  profile) and 'y' being parallel to the y
%                                  axis (N-S profile)
%
% OUPUT:
%   profile_axis: distance along the profile, a vector
%   profile:      data [i.e. gravity] values along the profile, a vector
% ======================================================================= %

function [profile_axis,profile] = data_profile(X,Y,grid,slice_position,slice_orientation)
    if (strcmp(slice_orientation, 'x') == 1)
        profile = interp2(X,Y,grid,X,slice_position);
        profile_axis = X;
    elseif (strcmp(slice_orientation, 'y') == 1)
        profile = interp2(X,Y,grid,slice_position,Y);
        profile_axis = Y;
    end
end
