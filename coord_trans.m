%                       BayGrav3D, Version 1.0 
%       © 2021 Erin Hightower, California Institute of Technology
%   ----------------------------------------------------------------------
%   Supported in part by the US National Science Foundation
%   (http://www.nsf.gov).
%
%   Free for non-commercial academic use only.
%
%   ----------------------------------------------------------------------
%   This file [coord_trans.m] is part of BayGrav3D.
%
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
% coord_trans.m
% -------------
% Rotates data coordinates so that the gravity grid is aligned with x and
% y axes (i.e. with due N and E axis) and hence with the axes of the 
% ensuing model domain, using two reference points. This ensures that the 
% model domain is entirely covered by the gravity grid and aligned in a 
% cartesian geometry. Depending on your data, this step is optional and 
% should be applied only in the pre- and post-processing. 

% *** All data, including both gravity and priors, need to be projected and
% rotated by the same amount so that everything is in the same reference
% coordinate system ***
%
% INPUT: 
%   data: x,y,[g,z,etc.]: array to be rotated
%   refpoints: x,y: arrray of at least 2 points forming a line parallel to 
%                   an edge of the grid
%   rot: 'ccw' or 'cw': direction in which to rotate the coordinates, 
%                       counterclockwise or clockwise
%
% OUTPUT:
%   newdata: x,y,[g,z,etc.] array of the data in the rotated coordinates
% ======================================================================= %

function [newdata] = coord_trans(data,refpoints,rot)
    data_xy = [data(:,1) data(:,2)]';
    minX = min(refpoints(:,1));
    maxX = max(refpoints(:,1));
    minY = min(refpoints(:,2));
    maxY = max(refpoints(:,2));
    
    if (strcmp(rot,'ccw') == 1)
        % Angle from horizontal by which to rotate coords counterclockwise
        theta = atan((maxY-minY)/(maxX-minX));
    elseif (strcmp(rot,'cw') == 1)
        % OR angle from horizontal by which to rotate coords clockwise
        theta = -atan((maxY-minY)/(maxX-minX));
    end
    
    % Rotate the coordinates with the rotation matrix
    R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
    rotated_coord = (R*data_xy)';
    
    % Assign new coordinates to the data array and return new data
    newdata = data;
    newdata(:,1) = rotated_coord(:,1);
    newdata(:,2) = rotated_coord(:,2);
end
    
    