%                       BayGrav3D, Version 1.0 
%       © 2021 Erin Hightower, California Institute of Technology
%   ----------------------------------------------------------------------
%   Supported in part by the US National Science Foundation
%   (http://www.nsf.gov).
%
%   Free for non-commercial academic use only.
%
%   ----------------------------------------------------------------------
%   This file [tikhonov_reg.m] is part of BayGrav3D.

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
% tikhonov_reg.m
% --------------
% Constructs the complete tikhonov regularization matrix used to stabilize
% the gravity inversion. Computes the finite difference approximation to
% the first or second derivative of the model parameters, based on the grid
% geometry, for each direction (x, y, z) separately. Also constructs the
% first order boundary condition regularization matrices applied only to
% the edge parameters. The final output L matrix is a linear combination of
% all the uni-directional and boundary L matrices. 
%
% INPUT:
%   num_param: total number of model parameters
%   Reg:       regularization type, one of the following:
%              0:   Zeroth Order
%              111: 1st order regularization in all directions
%              222: 2nd order regularization in all directions
%              221: 2nd order regularization in horz; 1st order in vertical
%
%   x_c:       x-coordinates of prism centroids
%   y_c:       y-coordinates of prism centroids
%   z_c:       z-coordinates of prism centroids
%
%   a:         Tikhonov regularization weight in the x-direction
%   b:         Tikhonov regularization weight in the y-direction
%   c:         Tikhonov regularization weight in the z-direction
%   A:         Boundary condition regularization weight in the x-direction
%   B:         Boundary condition regularization weight in the y-direction
%              Boundary condition weight must be very large: => 1e8
% OUTPUT:
%   L:         M x M regularization matrix used in gravity inversion
% ======================================================================= %

function [L] = tikhonov_reg(Reg,x_c,y_c,z_c,a,b,c,A,B)
    num_param = length(x_c)*length(y_c)*length(z_c);
    
    % Calculate distance between centroids in x,y,z to normalize Tikhonov
    % D[x,y,z]1 refers to first order delta [x,y,z] 
    % D[x,y,z]2 refers to second order delta [x,y,z]
    
    % X-DIRECTION
    xdiff1 = zeros(length(x_c)-1,1);
        for i = 1:length(x_c)-1
            xdiff1(i) = x_c(i+1)-x_c(i);
        end
    xdiff2 = zeros(length(xdiff1)-1,1);
        for i = 1:length(xdiff1)-1
            xdiff2(i) = 0.5*(xdiff1(i+1)+xdiff1(i));
        end
    xdiff1 = repmat(xdiff1,1,length(z_c)*length(y_c));
    xdiff2 = repmat(xdiff2,1,length(z_c)*length(y_c));
    Dx1 = reshape(xdiff1',1,[]);
    Dx2 = reshape(xdiff2',1,[]);

    % Y-DIRECTION
    ydiff1 = zeros(length(y_c)-1,1);
        for j = 1:length(y_c)-1
            ydiff1(j) = y_c(j+1)-y_c(j);
        end
    ydiff2 = zeros(length(ydiff1)-1,1);
        for j = 1:length(ydiff1)-1
            ydiff2(j) = 0.5*(ydiff1(j+1)+ydiff1(j));
        end
    ydiff1 = repmat(ydiff1,1,length(z_c));        % repeat same ydiff as we loop over z
    ydiff1 = vertcat(ydiff1,ones(1,length(z_c))); % add one row of ones to divide from the row of L that we zero out
    ydiff1 = reshape(ydiff1',1,[]);               % reshape into a vector 
    ydiff1 = repmat(ydiff1,length(x_c),1);        % replicate this vector to repeat the sequence for each new x-value

    ydiff2 = repmat(ydiff2,1,length(z_c));
    ydiff2 = vertcat(ydiff2,ones(2,length(z_c))); % add two rows of ones to divide from the rows of L that we zero out
    ydiff2 = reshape(ydiff2',1,[]); 
    ydiff2 = repmat(ydiff2,length(x_c),1);

    Dy1 = reshape(ydiff1',1,[]);        % reshape back into a vector
    Dy1(end-length(z_c)+1:end) = [];    % truncate the last length(z_c) rows because we need vector to be num_param-length(z_c) long 

    Dy2 = reshape(ydiff2',1,[]);        % reshape back into a vector
    Dy2(end-2*length(z_c)+1:end) = [];  % truncate the last 2*length(z_c) rows because vector needs to be num_param-2*length(z_c) long

    % Z-DIRECTION
    zdiff1 = zeros(length(z_c)-1,1);
        for k = 1:length(z_c)-1
            zdiff1(k) = z_c(k+1)-z_c(k);
        end
    zdiff2 = zeros(length(zdiff1)-1,1);
        for k = 1:length(zdiff1)-1
            zdiff2(k) = 0.5*(zdiff1(k+1)+zdiff1(k));
        end
    zdiff1 = horzcat(zdiff1',1);
    zdiff1 = repmat(zdiff1,length(x_c)*length(y_c)*length(z_c)/length(zdiff1),1);

    zdiff2 = horzcat(zdiff2',ones(1,2));
    zdiff2 = repmat(zdiff2,(length(y_c)*length(x_c)),1);

    Dz1 = reshape(zdiff1',1,[]);
    Dz1(end) = [];

    Dz2 = reshape(zdiff2',1,[]);
    Dz2(end-1:end) = [];
    
    % Create 1st order Tikhonov Reg. matrices and horizontal boundary reg. matrices
    % FOR X-DIRECTION
        % The second diagonal is placed on the length(z_c)*length(y_c)+1 
        % column because the model parameters up to this element all occur 
        % at the same x-value
        % Each length(z_c)*length(y_c) block of rows is divided by its 
        % corresponding delta x value, which are stored in the vector Dx1
        % L1x is an M-length(z_c)*length(y_c) X M matrix
        
        l1x = diag(-1*ones(1,num_param),0);
        l2x = diag(ones(1,num_param-(length(z_c))*(length(y_c))),(length(z_c))*(length(y_c)));
        L1x = l1x + l2x;
        L1x(num_param-(length(z_c))*(length(y_c))+1:end,:) = [];
        L1x = sparse(L1x);
        L1x = spdiags(1./Dx1',0,num_param-length(z_c)*length(y_c),num_param-length(z_c)*length(y_c))*L1x;

        % Split L1x to apply boundary condition in the x-direction
        Bx = full(L1x);
        Bx(length(z_c)*length(y_c)+1:end-length(z_c)*length(y_c),:) = 0;
        Bx = sparse(Bx);
    
    % FOR Y-DIRECTION
        % The second diagonal is placed on the length(z_c)+1 column because
        % model parameters up to this element all occur at the same y-value
        % Every length(z_c)*(length(y_c)-1) rows we skip length(z_c) rows 
        % (make zero) so that non-adjacent parameters are not regularized
        % Each length(z_c)*(length(y_c)-1) rows are divided by thier 
        % corresponding delta y value, which are stored in the vector Dy1
        % L1y is an M-length(z_c) X M matrix

        l1y = diag(-1*ones(1,num_param),0);
        l2y = diag(ones(1,num_param-length(z_c)),length(z_c));
        L1y = l1y + l2y;
        L1y(num_param-length(z_c)+1:end,:) = [];
        for i = 1:length(z_c)
            L1y(length(z_c)*(length(y_c)-1)+i:length(z_c)*(length(y_c)-1)+length(z_c):end,:) = 0;
        end
        L1y = sparse(L1y);
        L1y = spdiags(1./Dy1',0,num_param-length(z_c),num_param-length(z_c))*L1y;

        % Boundary condition in the y-direction
        By = full(L1y);
        for i = 1:length(z_c)*(length(y_c)-3)
            By(length(z_c)+i:length(z_c)*(length(y_c)):end,:) = 0;
        end
        By = sparse(By);
    
    % FOR Z-DIRECTION
        % The second diagonal is placed on the second column because the 
        % parameters are adjacent to one another in the z-direction
        % Every length(z_c)-1 rows, we skip (zero) 1 row so that we don't 
        % regularize parameters at the top and bottom with each other
        % Each non-zero block of length(z_c)-1 rows is divided by its 
        % corresponding delta z value, which are stored in the vector Dz1
        % L1z is an M-1 X M matrix
        
        l1z = diag(-1*ones(1,num_param),0);
        l2z = diag(ones(1,num_param-1),1);
        L1z = l1z + l2z;
        L1z(end,:) = [];
        L1z(length(z_c):length(z_c):end,:) = 0;
        L1z = sparse(L1z);
        L1z = spdiags(1./Dz1',0,num_param-1,num_param-1)*L1z;
    
    % Choose regularization type/combination and return L matrix
    if Reg == 0
        % Zeroth Order Regularization 
        Lm = speye(num_param,num_param); 
        L = (a^2)*(Lm'*Lm) + (A^2)*(Bx'*Bx) + (B^2)*(By'*By);
        
    elseif Reg == 111 
        % First Order Regularization in All Directions
        % Zero out regularization between boundary and adjacent prisms in 
        % the x-direction: Bx is used for those prisms instead
        L1x = full(L1x);
        L1x(1:length(z_c)*length(y_c),:) = 0;
        L1x(end-length(z_c)*length(y_c)+1:end,:) = 0;
        L1x = sparse(L1x);
        
        % Zero out regularization between boundary and adjacent prisms in 
        % the y-direction: By is used for those prisms instead
        % L1y is already length(z_c)*(length(y_c)-1), so to zero out 
        % boundary prisms start with length(z_c)*(length(y_c)-2)
        % Loop over 3*length(z_c) to elimnate 1) North boundary reg, 2) reg
        % between North and South (already zeroed out above), and 3) South
        % boundary reg
        % Eliminate boundary reg on first length(z_c) pairs of parameters
        L1y = full(L1y);
        for i = 1:3*length(z_c)
            L1y(length(z_c)*(length(y_c)-2)+i:length(z_c)*(length(y_c)-2)+(2*length(z_c)):end,:) = 0;
        end
        L1y(1:length(z_c),:) = 0;
        L1y = sparse(L1y);
        
        % Linearly combine x,y,z reg and boundary reg into single matrix L
        L = (a^2)*(L1x'*L1x) + (b^2)*(L1y'*L1y) + (c^2)*(L1z'*L1z) + (A^2)*(Bx'*Bx) + (B^2)*(By'*By);
        
    elseif Reg == 222  
        % Second Order Regularization in All Directions
        % FOR X-DIRECTION
            % The second diagonal is on the length(z_c)*length(y_c)+1 column
            % The third diagonal is on the 2*length(z_c)*length(y_c)+1 column
            % Each length(z_c)*length(y_c) block of rows in L1x is subtracted 
            % from the subsequent block of length(z_c)*length(y_c) rows in L1x 
            % Zero out regularization between boundary and adjacent prisms 
            % in the x-direction: Bx is used for those prisms  
            % Each of the resulting length(z_c)*length(y_c) blocks of rows 
            % are divided by their corresponding del x value, which are 
            % stored in the vector Dx2
            % L2x is an M-2*length(z_c)*length(x_c) X M matrix
            L1x = full(L1x);
            idx = 0:length(x_c)-2;
            L2x = zeros(num_param-2*length(z_c)*length(y_c),num_param);
            for i = 2:length(x_c)-1
                L2idx = L1x(((1 + idx(i)*length(z_c)*length(y_c)):(length(z_c)*length(y_c) + idx(i)*length(z_c)*length(y_c))),:) - ...
                    L1x(((1 + (idx(i)-1)*length(z_c)*length(y_c)):(length(z_c)*length(y_c) + (idx(i)-1)*length(z_c)*length(y_c))),:);
                L2x(((1 + (idx(i)-1)*length(z_c)*length(y_c)):(length(z_c)*length(y_c) + (idx(i)-1)*length(z_c)*length(y_c))),:) = L2idx;
            end
            clear L1x

            L2x(1:length(z_c)*length(y_c),:) = 0;
            L2x(end-length(z_c)*length(y_c)+1:end,:) = 0;

            L2x = sparse(L2x);
            L2x = spdiags(1./Dx2',0,num_param-2*length(z_c)*length(y_c),num_param-2*length(z_c)*length(y_c))*L2x;
        
        % FOR Y-DIRECTION 
            % The second diagonal is placed on the length(z_c)+1 column
            % The third diagonal is placed on the 2*length(z_c)+1 column
            % Each block of length(z_c) rows is subtracted from the 
            % subsequent block of length(z_c) rows
            % Every length(z_c)*(length(y_c)-2) rows we skip (zero out) 2*length(z_c) rows
            % Zero out regularzation between boundary and adjacent prisms 
            % in the y-direction: By is used for those prisms
            % Each of the resulting non-zero block of length(z_c) rows is 
            % divided by its corresponding del y value, which are stored in
            % the diagonal matrix Dy2
            % L2y is an M-2*length(z_c) X M matrix 
            L1y = full(L1y);
            idy = 0:length(x_c)*length(y_c)-2;
            L2y = zeros(num_param-2*length(z_c),num_param);
            for j = 2:(length(y_c)*length(x_c))-1
                L2idy = L1y(((1 + idy(j)*length(z_c)):(length(z_c) + idy(j)*length(z_c))),:) - ...
                    L1y(((1 + (idy(j)-1)*length(z_c)):(length(z_c) + (idy(j)-1)*length(z_c))),:);
                L2y(((1 + (idy(j)-1)*length(z_c)):(length(z_c) + (idy(j)-1)*length(z_c))),:) = L2idy;
            end
            clear L1y

            for j = 1:2*length(z_c)
                L2y(length(z_c)*(length(y_c)-2)+j:length(z_c)*(length(y_c)-2)+2*length(z_c):end,:) = 0;
            end

            % L2y is already length(z_c)*(length(y_c)-2), so to zero out 
            % boundary prisms start with length(z_c)*(length(y_c)-3)
            % Loop over 4*length(z_c) to elimnate 1) North boundary reg, 
            % 2) 2 regs between North and South (already zeroed out above),
            % and 3) South boundary reg
            % Eliminate boundary reg on the first length(z_c) triplets of parameters
            for i = 1:4*length(z_c)
            L2y(length(z_c)*(length(y_c)-3)+i:length(z_c)*(length(y_c)-3)+(3*length(z_c)):end,:) = 0;
            end
            L2y(1:length(z_c),:) = 0;

            L2y = sparse(L2y);
            L2y = spdiags(1./Dy2',0,num_param-2*length(z_c),num_param-2*length(z_c))*L2y;
        
        % FOR Z-DIRECTION
            % The second diagonal is placed in the second column
            % The third diagonal is placed on the third column 
            % Each row is subtracted from the subsequent row 
            % Every length(z_c)-2 rows we skip (zero out) two rows 
            % Each of the resulting non-zero blocks of length(z_c)-2 rows 
            % is divided by its corresponding del z value, which are stored
            % in the diagonal matrix Dz2
            % L2z is an M-2 X M matrix
            L1z = full(L1z);
            L1z_add = L1z(2:end,:);
            L1z_add = vertcat(L1z_add,zeros(1,size(L1z_add,2)));
            L2z = L1z_add - L1z;
            for k = 1:2
                L2z((length(z_c)-2)+k:(length(z_c)-2)+2:end,:) = 0;
            end
            L2z(end,:) = [];
            L2z = sparse(L2z);
            L2z = spdiags(1./Dz2',0,num_param-2,num_param-2)*L2z;
        
        L = (a^2)*(L2x'*L2x) + (b^2)*(L2y'*L2y) + (c^2)*(L2z'*L2z) + (A^2)*(Bx'*Bx) + (B^2)*(By'*By);
        
    elseif Reg == 221
        % Second Order Regularization in Horizontal; First Order in Vertical
        % FOR X-DIRECTION
            % The second diagonal is on the length(z_c)*length(y_c)+1 column
            % The third diagonal is on the 2*length(z_c)*length(y_c)+1 column
            % Each length(z_c)*length(y_c) block of rows in L1x is subtracted 
            % from the subsequent block of length(z_c)*length(y_c) rows in L1x 
            % Zero out regularization between boundary and adjacent prisms 
            % in the x-direction: Bx is used for those prisms  
            % Each of the resulting length(z_c)*length(y_c) blocks of rows 
            % are divided by their corresponding del x value, which are 
            % stored in the vector Dx2
            % L2x is an M-2*length(z_c)*length(x_c) X M matrix
            L1x = full(L1x);
            idx = 0:length(x_c)-2;
            L2x = zeros(num_param-2*length(z_c)*length(y_c),num_param);
            for i = 2:length(x_c)-1
                L2idx = L1x(((1 + idx(i)*length(z_c)*length(y_c)):(length(z_c)*length(y_c) + idx(i)*length(z_c)*length(y_c))),:) - ...
                    L1x(((1 + (idx(i)-1)*length(z_c)*length(y_c)):(length(z_c)*length(y_c) + (idx(i)-1)*length(z_c)*length(y_c))),:);
                L2x(((1 + (idx(i)-1)*length(z_c)*length(y_c)):(length(z_c)*length(y_c) + (idx(i)-1)*length(z_c)*length(y_c))),:) = L2idx;
            end
            clear L1x

            L2x(1:length(z_c)*length(y_c),:) = 0;
            L2x(end-length(z_c)*length(y_c)+1:end,:) = 0;

            L2x = sparse(L2x);
            L2x = spdiags(1./Dx2',0,num_param-2*length(z_c)*length(y_c),num_param-2*length(z_c)*length(y_c))*L2x;
        
        % FOR Y-DIRECTION 
            % The second diagonal is placed on the length(z_c)+1 column
            % The third diagonal is placed on the 2*length(z_c)+1 column
            % Each block of length(z_c) rows is subtracted from the 
            % subsequent block of length(z_c) rows
            % Every length(z_c)*(length(y_c)-2) rows we skip (zero) 2*length(z_c) rows
            % Zero out regularzation between boundary and adjacent prisms 
            % in the y-direction: By is used for those prisms
            % Each of the resulting non-zero block of length(z_c) rows is 
            % divided by its corresponding del y value, which are stored 
            % in the diagonal matrix Dy2
            % L2y is an M-2*length(z_c) X M matrix 
            L1y = full(L1y);
            idy = 0:length(x_c)*length(y_c)-2;
            L2y = zeros(num_param-2*length(z_c),num_param);
            for j = 2:(length(y_c)*length(x_c))-1
                L2idy = L1y(((1 + idy(j)*length(z_c)):(length(z_c) + idy(j)*length(z_c))),:) - ...
                    L1y(((1 + (idy(j)-1)*length(z_c)):(length(z_c) + (idy(j)-1)*length(z_c))),:);
                L2y(((1 + (idy(j)-1)*length(z_c)):(length(z_c) + (idy(j)-1)*length(z_c))),:) = L2idy;
            end
            clear L1y

            for j = 1:2*length(z_c)
                L2y(length(z_c)*(length(y_c)-2)+j:length(z_c)*(length(y_c)-2)+2*length(z_c):end,:) = 0;
            end

            % L2y is already length(z_c)*(length(y_c)-2), so to zero out 
            % boundary prisms start with length(z_c)*(length(y_c)-3)
            % Loop over 4*length(z_c) to elimnate 1) North boundary reg, 
            % 2) 2 regs between North and South (already zeroed out above),
            % and 3) South boundary reg
            % Eliminate boundary reg on the first length(z_c) triplets of parameters
            for i = 1:4*length(z_c)
            L2y(length(z_c)*(length(y_c)-3)+i:length(z_c)*(length(y_c)-3)+(3*length(z_c)):end,:) = 0;
            end
            L2y(1:length(z_c),:) = 0;

            L2y = sparse(L2y);
            L2y = spdiags(1./Dy2',0,num_param-2*length(z_c),num_param-2*length(z_c))*L2y;
        
        % FOR Z-DIRECTION: use L1z already computed above
        % Combine component-wise matrices into a single reg. matrix
        L = (a^2)*(L2x'*L2x) + (b^2)*(L2y'*L2y) + (c^2)*(L1z'*L1z) + (A^2)*(Bx'*Bx) + (B^2)*(By'*By);
    end
end
