%                       BayGrav3D, Version 1.0 
%       © 2021 Erin Hightower, California Institute of Technology
%   ----------------------------------------------------------------------
%   Supported in part by the US National Science Foundation
%   (http://www.nsf.gov).
%
%   Free for non-commercial academic use only.
%
%   ----------------------------------------------------------------------
%   This file [preprocessing.m] is part of BayGrav3D.

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
% preprocessing.m
% ---------------
% This is a preprocessing script to:
%   1. rotate and transform all data files into the same reference frame, 
%      with the gravity grid aligned with the x,y axes and the x and y 
%      extents of the grid shifted to a 0,0 origin. 
%   2. Generates the mesh using the spatial extent of the gravity grid as
%      reference. 
%   3. Prints a file with all the input data file names and mesh parameters
%      used to generate this model
% 
% ***MODIFY SCRIPT AS NEEDED TO SUITE YOUR NEEDS
% 
% Necessary output for use with the gravity inversion:
%   *** All coordinates must be projected to cartesian x, y distances in m
% 
%   data_gravity:    [x,y,g] text file of gravity values to be inverted
%   prism_centroids: [x,y,z] text file of centroid of each prism in mesh
%   _mesh.mat:       .mat file produced by construct_mesh.m that contains 
%                    the nodal and centroid coordinate vectors of the mesh:
%                    [x,y,z,x_c,y_c,z_c,prism_centroids]
%   input_data_file: text file listing the input files and mesh paramters
%                    used to construct a certain prior and the mesh for 
%                    this model
% ======================================================================= %

filepath = 'BayGrav3D/Example1/';
modelname = 'Example1';

% Reference points with which to rotate the data coordinates
input_refpoints = [filepath 'rotational_ref_line_SISIEline02.txt'];

% Gravity data to import
input_gravity = [filepath 'data_gravity_sr3.7_20_EN.csv'];

% Optional data fields to transform 
input_structure = [filepath 'data_priors_horizons_MCSbase_all_w_SedThick.txt'];
input_bathy = [filepath 'data_priors_bathymetry.txt'];
input_density = [filepath 'data_priors_density_all_6-30-20_UTM.txt'];
input_velocity = '';
input_quakes = '';
input_coast = '';

% Define parameters for generating the mesh 
mesh_lim_method = 'fromData'; % either 'fromData' or 'manual'
mesh_lims = [];
mesh_style = 'var';
np_x = 20;
np_y = 20;
np_z = 20;
res_z = [5 100];
input_Zmax = 25e3;

% ======================================================================= %
% Import reference points to use for coordinate transformation
% ============================================================
refpoints = readtable(input_refpoints,'ReadVariableNames',1);
refpoints.formation = []; % comment out depending on your file
refpoints = table2array(refpoints);

% ======================================================================= %
% Import and transform gravity data coordinates 
% ============================================= 
data_gravity = readtable(input_gravity,'ReadVariableNames',1);
data_gravity = table2array(data_gravity);
data_gravity = coord_trans(data_gravity,refpoints,'ccw');

shiftx = min(data_gravity(:,1));
shifty = min(data_gravity(:,2));

data_gravity(:,1) = data_gravity(:,1) - shiftx;
data_gravity(:,2) = data_gravity(:,2) - shifty;

fid = fopen([filepath modelname '_data_gravity_preprocessed.txt'],'w');
fprintf(fid,'%s %s %s\n','x','y','g');
fprintf(fid,'%.6f %.6f %.6f\n',data_gravity');
fclose(fid);

% ======================================================================= %
% Import and transform structural horizons from seismic/well data for prior
% ======================================================================= %
if strcmp(input_structure,'') == 1
    disp('No structural data to transform')
else
    data_structure = readtable(input_structure,'ReadVariableNames',1);
    data_structure.formation = [];
    data_structure = table2array(data_structure);

    data_structure = coord_trans(data_structure,refpoints,'ccw');
    data_structure(:,1) = data_structure(:,1) - shiftx;
    data_structure(:,2) = data_structure(:,2) - shifty;
    
    fid = fopen([filepath modelname '_data_structure_preprocessed.txt'],'w');
    fprintf(fid,'%s %s %s %s\n','x','y','z','horizonID');
    fprintf(fid,'%.6f %.6f %.6f %2d\n', data_structure');
    fclose(fid);
end

% ======================================================================= %
% Import and transform bathymetry data for prior 
% ==============================================
if strcmp(input_bathy,'') == 1
    disp('No bathymetry data to transform')
else
    data_bathy = readtable(input_bathy,'ReadVariableNames',1);
    data_bathy.formation = [];
    data_bathy = table2array(data_bathy);

    data_bathy = coord_trans(data_bathy,refpoints,'ccw');
    data_bathy(:,1) = data_bathy(:,1) - shiftx;
    data_bathy(:,2) = data_bathy(:,2) - shifty;
    
    fid = fopen([filepath modelname '_data_bathy_preprocessed.txt'],'w');
    fprintf(fid,'%s %s %s %s\n','x','y','z','horizonID');
    fprintf(fid,'%.6f %.6f %.6f %2d\n',data_bathy');
    fclose(fid);
end

% ======================================================================= %
% Import densities (from seismic velocities) for the prior 
% ========================================================
if strcmp(input_density,'') == 1
    disp('No density data to transform')
else
    data_density = readtable(input_density,'ReadVariableNames',1);
    data_density = table2array(data_density);
    data_density(:,4) = data_density(:,4).*1e3; % convert to kg/m^3
    data_density(:,5) = data_density(:,5).*1e3;

    data_density = coord_trans(data_density,refpoints,'ccw');
    data_density(:,1) = data_density(:,1) - shiftx;
    data_density(:,2) = data_density(:,2) - shifty; 
    
    fid = fopen([filepath modelname '_data_density_preprocessed.txt'],'w');
    fprintf(fid,'%s %s %s %s %s\n','x','y','z','density','stdev');
    fprintf(fid,'%.6f %.6f %.6f %.4f %.4f\n',data_density');
    fclose(fid);
end

% ======================================================================= %
% Import earthquake data to use in the prior 
% ==========================================
if strcmp(input_quakes,'') == 1
    disp('No earthquake location data to transform')
else
    data_quakes = readtable(input_quakes,'ReadVariableNames',1);
    data_quakes = table2array(data_quakes);
    data_quakes(:,3) = data_quakes(:,3)*1e3; % convert to m depth

    data_quakes = coord_trans(data_quakes,refpoints,'ccw');
    data_quakes(:,1) = data_quakes(:,1) - shiftx;
    data_quakes(:,2) = data_quakes(:,2) - shifty;
    
    fid = fopen([filepath modelname '_data_earthquakes_preprocessed.txt'],'w');
    fprintf(fid,'%s %s %s %s\n','x','y','z','mag');
    fprintf(fid,'%.6f %.6f %.6f %.2f\n',data_quakes');
    fclose(fid);
end

% ======================================================================= %
% Import coastline data to use in the prior 
% =========================================
if strcmp(input_coast,'') == 1
    disp('No coastline data to transform')
else
    data_coast = readtable(input_coast);
    data_coast.Var3 = [];
    data_coast = table2array(data_coast);

    data_coast = coord_trans(data_coast,refpoints,'ccw');
    data_coast(:,1) = data_coast(:,1) - shiftx;
    data_coast(:,2) = data_coast(:,2) - shifty;
    
    fid = fopen([filepath modelname '_data_coastline_preprocessed.txt'],'w');
    fprintf(fid,'%s %s\n','x','y');
    fprintf(fid,'%.6f %.6f\n',data_coast');
    fclose(fid);
end

% ======================================================================= %
% Set spatial limits of the subsurface domain and generate mesh
% =============================================================

% set spatial limits of the mesh to be generated
% Zmin must be > 0 to stabilize G under the possibility of nearly coincident data and prism corner points
if strcmp(mesh_lim_method,'fromData') == 1
    Xmin = min(data_gravity(:,1));
    Xmax = max(data_gravity(:,1));
    Ymin = min(data_gravity(:,2));
    Ymax = max(data_gravity(:,2));
    Zmin = 1;
    Zmax = input_Zmax;
elseif strcmp(mesh_lim_method,'manual') == 1
    Xmin = mesh_lims(1);
    Xmax = mesh_lims(2);
    Ymin = mesh_lims(3);
    Ymax = mesh_lims(4);
    Zmin = 1;
    Zmax = input_Zmax;
end

% generate the mesh node and centroid coordinates using construct_mesh.m
[x,y,z,x_c,y_c,z_c,prism_centroids] = construct_mesh(Xmin,Xmax,Ymin,Ymax,Zmin,Zmax,mesh_style,np_x-2,np_y-2,np_z,res_z,[filepath modelname]);

% ======================================================================= %
% Save prior data file names and mesh parameters to file for reference
% ====================================================================
input_data_file = [filepath modelname '_InputData-and-MeshParameters.txt'];
fid = fopen(input_data_file,'w');
fprintf(fid,'%s\n', 'INPUT DATA FILES FOR GRAVITY:');
fprintf(fid,'%s\n\t%s%s\n', 'gravity', '= ', input_gravity);
fprintf(fid,'\n');
fprintf(fid,'%s\n', 'INPUT DATA FILES FOR CONSTRUCTING THE PRIOR:');
fprintf(fid,'%s\n\t%s%s\n', 'structure', '= ', input_structure);
fprintf(fid,'%s\n\t%s%s\n', 'bathymetry', '= ', input_bathy);
fprintf(fid,'%s\n\t%s%s\n', 'density', '= ', input_density);
fprintf(fid,'%s\n\t%s%s\n', 'velocity', '= ', input_velocity);
fprintf(fid,'%s\n\t%s%s\n', 'earthquakes', '= ', input_quakes);
fprintf(fid,'%s\n\t%s%s\n', 'coastline', '= ', input_coast);
fprintf(fid,'\n');
fprintf(fid,'%s\n', 'INPUT REFERENCE POINTS AND SHIFT VALUE FOR TRANSFORMING DATA');
fprintf(fid,'%s\n\t%s%s\n', 'coord rotation refpoints', '= ', input_refpoints);
fprintf(fid,'%s\n\t%s%s%s\n', 'data shiftx, shifty', '= ', [num2str(shiftx) ','],num2str(shifty));
fprintf(fid,'\n');
fprintf(fid,'%s\n', 'INPUT MESH PARAMETERS FOR CONSTRUCTING THE GRID');
fprintf(fid,'%s\n\t%s%s%s%s\n', 'mesh x,y,z number of prisms', '= ', [num2str(np_x) ', '],[num2str(np_y) ', '],num2str(np_z));
fprintf(fid,'%s\n\t%s%s\n', 'mesh min,max vertical resolution', '= ', num2str(res_z));
fprintf(fid,'%s\n\t%s%s\n', 'mesh style', '= ', mesh_style);
fprintf(fid,'%s\n\t%s%s%s%s%s%s%s\n', 'mesh limits', '= ', ...
    [num2str(Xmin) ', '],[num2str(Xmax) ', '],[num2str(Ymin) ', '],[num2str(Ymax) ', '],[num2str(Zmin) ', '],num2str(Zmax));
fclose(fid);
