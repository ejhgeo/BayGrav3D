%                       BayGrav3D, Version 1.0 
%       © 2021 Erin Hightower, California Institute of Technology
%   ----------------------------------------------------------------------
%   Supported in part by the US National Science Foundation
%   (http://www.nsf.gov).
%
%   Free for non-commercial academic use only.
%
%   ----------------------------------------------------------------------
%   This file [vec2grid.m] is part of BayGrav3D.

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
% vec2grid.m
% ----------
% Reshapes vectors into 2D or 3D grids as necessary to be plotted. 
%
% INPUT
%   dim: 2 or 3, for 2D or 3D grid, respectively (2D for gravity, 3D for
%        density, or other mesh values)
%   vec: vector to reshape
%   x:   x coordinates of the resulting grid (e.g. Xdata or x_c)
%   y:   y coordinates of the resulting grid (e.g. Ydata or y_c)
%   z:   z coordinates of the resulting grid (e.g. z_c)
%   res: resolution of desired 2D grid (only for gravity grids)
%
% OUTPUT
%   X:    2D grid of X values (for gravity grids)
%   Y:    2D grid of Y values (for gravity grids)
%   grid: 2D or 3D resulting grid
% ======================================================================= %

function [X,Y,grid] = vec2grid(dim,vec,x,y,z,res)
    if dim == 2
        xlin = linspace(min(x),max(x),res);
        ylin = linspace(min(y),max(y),res);
        [Xgrid,Ygrid] = meshgrid(xlin,ylin);
        X = Xgrid(1,:);
        Y = Ygrid(:,1)';
        
        fgrid = scatteredInterpolant(x,y,vec);
        grid = fgrid(Xgrid,Ygrid);
        
    elseif dim == 3
        grid = permute(reshape(vec,length(z),length(y),length(x)), [3 2 1]);
        grid = permute(grid,[2 1 3]); 
        X = [];
        Y = [];
    end
end