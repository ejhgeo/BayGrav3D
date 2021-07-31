======================================
Gravitational Anomaly over a 3D Prism
======================================

The gravitational anomaly of a rectangular prism is given as in Turcotte and Schubert (2014) as:

.. math::

  \Delta g = \Gamma \Delta \rho \sum_{i=1}^{2} \sum_{j=1}^{2} \sum_{k=1}^{2} \mu_{ijk}[\Delta z_k arctan(\frac{\Delta x_i \Delta y_i}{\Delta z_k R_{ijk}}) - \Delta x_i ln(R_{ijk} + \Delta y_j) - \Delta y_j ln(R_{ijk} + \Delta x_i)]

where :math:`\Delta x_i = (x_i - x_p)`, :math:`\Delta y_j = (y_j - y_p)`, :math:`\Delta z_k = (z_k - z_p)`, and :math:`\mu_{ijk} = (-1)^i (-1)^j (-1)^k`. :math:`\Delta \rho` is the density contrast of the prism, and :math:`\Gamma` is the universal gravitational constant. :math:`x_p,  y_p,` and :math:`z_p` are the coordinates of the measurement point P, and :math:`x_i, y_j,` and :math:`z_k` are the coordinates of the corners of the prism, where :math:`(i, j, k) = (1, 2)`. :math:`R_{ijk}` is the distance from the measurement point to a corner at :math:`x_i, y_j, z_k` and is given by

.. math::

  R_{ijk} = (\Delta x_{i}^{2} + \Delta y_{j}^{2} + \Delta z_{k}^{2})^{1/2}

For computational purposes, this equation can be expanded, demonstrating that the gravitational attraction of the prism is essentially the sum of the gravitational attraction of all its sides. Expanding over the first sum gives:

.. math::

  \begin{align}
    \Delta g = \Gamma \Delta \rho \sum_{i=1}^{2} \sum_{j=1}^{2} [\mu_{ij1}[\Delta z_1 tan^{-1}(\frac{\Delta x_i \Delta y_j}{\Delta z_1 R_{ij1}}) - \Delta x_i ln(R_{ij1} + \Delta y_j) - \Delta y_j ln(R_{ij1} + \Delta x_i)] \nonumber \\ 
    + \mu_{ij2}[\Delta z_2 tan^{-1} (\frac{\Delta x_i \Delta y_j}{\Delta z_2 R_{ij2}}) - \Delta x_i ln(R_{ij2} + \Delta y_j) - \Delta y_j ln(R_{ij2} + \Delta x_i)]] \nonumber
  \end{align}

Expanding over the second sum gives:

.. math::
  \begin{align}
    \Delta g = \Gamma \Delta \rho \sum_{i=1}^{2}[\mu_{i11}[\Delta z_1 tan^{-1}(\frac{\Delta x_i \Delta y_1}{\Delta z_1 R_{i11}}) - \Delta x_i ln(R_{i11} + \Delta y_1) - \Delta y_1 ln(R_{i11} + \Delta x_i)] \nonumber \\
    + \mu_{i12}[\Delta z_2 tan^{-1}(\frac{\Delta x_i \Delta y_1}{\Delta z_2 R_{i12}}) - \Delta x_i ln(R_{i12} + \Delta y_1) - \Delta y_1 ln(R_{i12} + \Delta x_i)] \nonumber \\
    + \mu_{i21}[\Delta z_1 tan^{-1}(\frac{\Delta x_i \Delta y_2}{\Delta z_1 R_{i21}}) - \Delta x_i ln(R_{i21} + \Delta y_2) - \Delta y_2 ln(R_{i21} + \Delta x_i)] \nonumber \\
    + \mu_{i22}[\Delta z_2 tan^{-1}(\frac{\Delta x_i \Delta y_2}{\Delta z_2 R_{i22}}) - \Delta x_i ln(R_{i22} + \Delta y_2) - \Delta y_2 ln(R_{i22} + \Delta x_i)]] \nonumber
  \end{align}

Expanding over the final sum gives:

.. math::
  \begin{align}
    \Delta g = \Gamma \Delta \rho [\mu_{111}[\Delta z_1 tan^{-1}(\frac{\Delta x_1 \Delta y_1}{\Delta z_1 R_{111}}) - \Delta x_1 ln(R_{111}+\Delta y_1) - \Delta y_1 ln(R_{111}+\Delta x_1)] \nonumber \\
    + \mu_{112}[\Delta z_2 tan^{-1}(\frac{\Delta x_1 \Delta y_1}{\Delta z_2 R_{112}}) - \Delta x_1 ln(R_{112}+\Delta y_1) -\Delta y_1 ln(R_{112}+\Delta x_1)] \nonumber \\
    + \mu_{121}[\Delta z_1 tan^{-1}(\frac{\Delta x_1 \Delta y_2}{\Delta z_1 R_{121}}) - \Delta x_1 ln(R_{121} + \Delta y_2) - \Delta y_2 ln(R_{121}+\Delta x_1)] \nonumber \\
    + \mu_{122}[\Delta z_2 tan^{-1}(\frac{\Delta x_1 \Delta y_2}{\Delta z_2 R_{122}}) - \Delta x_1 ln(R_{122}+\Delta y_2) - \Delta y_2 ln(R_{122}+\Delta x_1)] \nonumber \\
    + \mu_{211}[\Delta z_1 tan^{-1}(\frac{\Delta x_2 \Delta y_1}{\Delta z_1 R_{211}}) - \Delta x_2 ln(R_{211}+\Delta y_1) - \Delta y_1 ln(R_{211}+\Delta x_2)] \nonumber \\
    + \mu_{212}[\Delta z_2 tan^{-1}(\frac{\Delta x_2 \Delta y_1}{\Delta z_2 R_{212}}) - \Delta x_2 ln(R_{212}+\Delta y_1) - \Delta y_1 ln(R_{212}+\Delta x_2)] \nonumber \\
    + \mu_{221}[\Delta z_1 tan^{-1}(\frac{\Delta x_2 \Delta y_2}{\Delta z_1 R_{221}}) - \Delta x_2 ln(R_{221}+\Delta y_2) - \Delta y_2 ln(R_{221}+\Delta x_2)] \nonumber \\ 
    + \mu_{222}[\Delta z_2 tan^{-1}(\frac{\Delta x_2 \Delta y_2}{\Delta z_2 R_{222}}) - \Delta x_2 ln(R_{222}+\Delta y_2) - \Delta y_2 ln(R_{222}+\Delta x_2)]] \nonumber
  \end{align}

The expanded sums together define the geometry of the prism relative to the observation point, which is a constant in the following inversion. In other words, each prism in the domain will have a single geometry coefficient for each observation point on the grid. 

To set-up this gravity formula for the inversion, we have to extend it to the case of multiple prisms. In the inversion, we wish to invert satellite gravity data of N observation points to obtain the best estimate of the density of M subsurface prisms (and hence M model parameters) according to the model defined above. 

First, define each :math:`\Delta x_1` as the difference between :math:`x_{p}` and :math:`x_{i}` such that :math:`\Delta x_1` is a 1 x N vector of differences between the prism's :math:`x_1` edge and each observation point. Similarly, define each :math:`\Delta x_2` as the difference between the next edge of the prism in the x direction and each observation point, such that it is also a 1 x N vector. This can be done iteratively for each prism such that :math:`\Delta x_1` and :math:`\Delta x_2` are each an M x N matrix, where each row corresponds to a prism and each column an observation point, the elements of which describe the difference in the x-direction between that prism's edge (1 or 2) and the observation point. For example, :math:`\Delta x_1(1,1)` is the difference between the first prism's first edge and the first observation point. Likewise, the same procedure is repeated for :math:`\Delta y_1, \Delta y_2, \Delta z_1, \Delta z_2`.

Once this geometry is defined, the operations in the expanded sum can then be carried out as element-wise matrix operations, such that there is an final M x N matrix that describes the set geometry of each prism relative to each observation point. We will call this matrix **Z**. The gravity anomaly at any observation point due to the combined effect of all the prisms is then the matrix product of **Z** and :math:`\mathbf{\Delta \rho}`, where :math:`\mathbf{\Delta \rho}` is an M x 1 vector containing the differential density of each prism. 

.. math::
  
  \Delta \mathbf{g} = \Gamma \mathbf{Z}^T \Delta \mathbf{\rho}

