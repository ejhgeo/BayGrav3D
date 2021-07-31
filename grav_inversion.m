%                       BayGrav3D, Version 1.0 
%       © 2021 Erin Hightower, California Institute of Technology
%   ----------------------------------------------------------------------
%   Supported in part by the US National Science Foundation
%   (http://www.nsf.gov).
%
%   Free for non-commercial academic use only.
%
%   ----------------------------------------------------------------------
%   This file [grav_inversion.m] is part of BayGrav3D.

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
% grav_inversion.m
% ----------------
% This is the main script in BayGrav3D that runs the full gravity
% inversion, provided the chosen input files and regularization parameters.
%
% INPUT:
%   filepath: path to directory where this model is being stored
%   modelname: name for the current model being inverted
%   input_mesh: .mat file storing the mesh to use, generated with
%               construct_mesh.m
%   input_parameter_file: .txt file storing the details of the data, mesh,
%                         and priors used for this model, generated with 
%                         preprocessing.m
%   inversion_runfile: name of the .txt file to store the input parameters 
%                      for this inversion, appended to input_parameter_file
%   input_gravity: file of the (preprocessed) gravity data to invert
%   input_data_stdev: chosen standard deviation of the data, scalar
%                     * This assumes a constant variance on the data,
%                     however, this can be a vector. If so, you will have
%                     to change the data covariance matrix setup
%                     accordingly. 
%   input_prior: .mat file storing the variables for the prior [mu_pr,
%                mu_abs, sigma_pr, Cm, prism_density_stdev, D, layer_mean,
%                layers] to use in the inversion
%   input_reg_type: chosen level of Tikhonov regularization,
%                   0: Zeroth Order Tikhonov Regularization
%                   111: First order reg. in all directions
%                   222: Second order reg. in all directions
%                   221: Second order reg. in x,y; first order reg. in z
%   input_TikW_X: Tikhonov regularization weight in x-direction
%   input_TikW_Y: Tikhonov regularization weight in y-direction
%   input_TikW_Z: Tikhonov regularization weight in z-direction
%
% OUTPUT:
%   ..._Inversion_RunFile.txt: file storing data and parameters used in
%                             this particular run of the inversion
%   ..._PosteriorModel.txt: array of the resulting model results for all
%                          mesh values:
%           [x_c y_c z_c m m_abs mu_pr mu_abs sigma_pr post_stdev post_res]
%   ..._InversionResults.mat: .mat file storing all results of the
%                             inversion, including model parameter 
%                             solution, prior values, posterior and prior 
%                             covariances, misfits, mesh coordinates, and 
%                             observed, predicted, and prior gravity
% ======================================================================= %

filepath = 'BayGrav3D/Example1/';
modelname = 'Example1';

% Specify mesh and model parameter files to read from
input_mesh = [filepath modelname '_mesh.mat'];
input_parameter_file = [filepath modelname '_InputData-and-MeshParameters.txt'];
inversion_runfile = [filepath modelname '_Inversion_RunFile.txt'];

% Specify gravity data file to import and gravity standard deviation
input_gravity = [filepath modelname '_data_gravity_preprocessed.txt'];
input_data_stdev = 1.7;

% Specify prior data to import
input_prior = [filepath modelname '_Prior_PuysegurGJI.mat'];

% Specify regularization type and the weight in each dimension
input_reg_type = 221;
input_TikW_X = 5e6;
input_TikW_Y = 5e6;
input_TikW_Z = 1e0;

% ======================================================================= %
% Create the inversion runfile
% ============================
% Print the time and date of the current inversion to write to file
inv_datetime = datestr(now,'YYYY-mm-dd_HH-MM');

% Create inversion runfile with date-time, chosen inversion parameters, and
% append input parameter file created during preprocessing
fid = fopen(inversion_runfile, 'wt');
fprintf(fid,'%s\n', inv_datetime);
fprintf(fid,'\n');
fprintf(fid,'%s\n\t%s%s\n', 'prior model', '= ', input_prior);
fprintf(fid,'%s\n\t%s%s\n', 'data stdev', '= ', num2str(input_data_stdev));
fprintf(fid,'%s\n\t%s%s\n', 'regularization type', '= ', num2str(input_reg_type));
fprintf(fid,'%s\n\t%s%s%s%s\n', 'Tikhonov weight in X,Y,Z', '= ', ...
        [num2str(input_TikW_X,'%1.0e') ', '], [num2str(input_TikW_Y,'%1.0e') ', '], num2str(input_TikW_Z,'%1.0e'));
fprintf(fid,'\n');
fclose(fid);

fr = fopen(input_parameter_file, 'rt');
fw = fopen(inversion_runfile, 'at');
while feof(fr) == 0 
    tline = fgetl(fr);
    fwrite(fw, sprintf('%s\n',tline));
end
fclose(fr);
fclose(fw);

% ======================================================================= %
% Load mesh, prior model, and gravity data
% ========================================
% Load the mesh and prior from mat file
load(input_mesh)
load(input_prior)

% Load input data
grav_data = readtable(input_gravity,'ReadVariableNames',1);
grav_data = table2array(grav_data);
Xdata = grav_data(:,1);
Ydata = grav_data(:,2);

% ======================================================================= %
% Create G matrix, calculate forward gravity on prior, and setup data
% ===================================================================
% Calculate G matrix (prism to observation point geometry)
disp('Calculating G matrix of prism to data point geometry')
[G,grav_corr] = construct_Gmatrix(x,y,z,x_c,y_c,z_c,grav_data);

% Calculate forward gravity of only the prior
grav_prior = G*mu_pr;

% Remove the effect of top 1 m of seawater from the gravity data
disp('Constructing d vector and data covariance matrix')
d = grav_data(:,3) - grav_corr;

% Construct the data covariance operator
stdev = (input_data_stdev^2).*ones(1,length(d));
Cd = spdiags(stdev',0,length(d),length(d));

% ======================================================================= %
% Construct the Tikhonov regularization matrix
% ============================================
% Assign weights for the Tikhonov regularization and calculate L matrix
disp('Calculating Tikhonov Regularization matrix')
A = 1e8; % weight for the boundary prisms in x-direction
B = 1e8; % weight for the boundary prisms in y-direction
a = input_TikW_X;
b = input_TikW_Y;
c = input_TikW_Z;

L = tikhonov_reg(input_reg_type,x_c,y_c,z_c,a,b,c,A,B);

% ======================================================================= %
% Perform linear least square inversion
% =====================================
disp('Performing Inversion')

% Calculate model parameter covariance matrix: inverse of the Hessian
C = (G'*(Cd^(-1))*G + L + (Cm^(-1)))^(-1);

% Calculate the gradient
Grad = G'*(Cd^(-1))*d + (Cm^(-1))*mu_pr;

% Solve for the model parameters and predicted gravity
m = C*Grad;
grav_pred = G*m;
clear G

% ======================================================================= %
% Calculate the posterior misfit and resolution
% =============================================
disp('Calculating model misfit and resolution')

% Calculate posterior standard deviation from C and resolution matrix R
post_stdev = sqrt(diag(C));
R = speye(length(m)) - C*(Cm^(-1));
post_res = diag(R);
clear R

% Calculate total model misfit from the complete misfit equation
F = 0.5 * (1/input_data_stdev^2)*((d-grav_pred)'*(d-grav_pred)) + ...
          m'*L*m + 0.5*sum(((mu_pr-m).^2)./diag(Cm));
      
% Calculate LS misfit on the gravity and parameters relative to the prior
Fg_LS = 0.5 * (1/input_data_stdev^2)*((d-grav_pred)'*(d-grav_pred));
Fm_LS = 0.5 * sum(((mu_pr - m).^2)./diag(Cm));

% Calculate mean absolute error on gravity and parameters relative to prior
Fg_MAE = (1/length(d)) * sum(abs(d-grav_pred));
Fm_MAE = (1/length(m)) * sum(abs(mu_pr - m));

disp('F = '); disp(F)
disp('Fg_MAE = '); disp(Fg_MAE)
disp('Fm_MAE = '); disp(Fm_MAE)
disp('Mean Posterior StDev = '); disp(mean(post_stdev))
disp('Mean Resolution = '); disp(mean(post_res));

% Convert differential density back into absolute density
m_abs = m + layer_mean;

% ======================================================================= %
% Save inversion results to .txt file and .mat file
% ==========================================

% Combine prism coordinates and model parameter values (includes boundary prisms)
posterior_model = [prism_centroids m m_abs mu_pr mu_abs sigma_pr post_stdev post_res]; 

% Save the output of the model to .txt file and .mat file
save([filepath modelname '_InversionResults.mat'],'d','grav_data','grav_prior','grav_pred',...
                                                  'prism_centroids','x','y','z','x_c','y_c','z_c',...
                                                  'mu_pr','mu_abs','sigma_pr','layers','layer_mean',...
                                                  'm','m_abs','post_stdev','post_res',...
                                                  'F','Fg_LS','Fm_LS','Fg_MAE','Fm_MAE');
                                              
fid = fopen([filepath modelname '_PosteriorModel.txt'],'w');
fprintf(fid,'%s %s %s %s %s %s %s %s %s %s\n','x_c','y_c','z_c','m','m_abs','mu_pr','mu_abs','sigma_pr','post_stdev','post_res');
fprintf(fid,'%.6f %.6f %.6f %.4f %.4f %.4f %.4f %.4f %.4f %.4f\n',posterior_model');
fclose(fid);




