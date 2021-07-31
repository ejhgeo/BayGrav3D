==============================
Linear Least Squares Inversion
==============================

The Linear Least Squares Solution
----------------------------------

The method for linear and nonlinear least squares inversion is covered in many texts, but we will be following the notation and convention of Aster et al. (2013) *Parameter Estimation and Inverse Problems*. Before implementing a full gravity inversion for the Puysegur region, we perform gravity inversions on synthetic gravity models for both a simple buried object and for a simplified subduction system, to set-up and test the appropriate equations and the performance of the inversion.

For N data points and M model parameters, where :math:`g_i(\textbf{m})` is the model prediction of the :math:`i^{th}` datum, the least squares misfit is defined as:

.. math::

  F(\textbf{m})=\frac{1}{2} \sum_{i=1}^{N} (d_i - g_i(\textbf{m}))^2

The Gauss-Newton algorithm for the solution of the model parameters that minimizes this least squares misfit is:

.. math::

    \textbf{m} = \textbf{m}_o + (\textbf{G}^T\textbf{G})^{-1} \textbf{G}^T (\textbf{d} - \textbf{g}(\textbf{m}_o))

Here, :math:`\textbf{G}= \frac{\partial g_i}{\partial m_k}` is an N x M matrix of the partial derivatives of **g**. For a case where the predicted data are a linear combination of the model parameters (i.e. the case where the derivative of the model with respect to the model parameters is not dependent on the model parameters), our prediction can be written directly as :math:`\mathbf{Gm}`, in which case the linear least squares solution reduces to:

.. math::
    \textbf{m} = (\textbf{G}^T\textbf{G})^{-1} \textbf{G}^T (\textbf{d})
    
We can also write this equation as the product of the Hessian, the second derivative of the misfit with respect to the model parameters (:math:`\mathbf{H} = \frac{\partial^2 F}{\partial m_k \partial m_k}`, and the gradient, :math:`\mathbf{\gamma}=\frac{\partial F}{\partial m_k}`, the first derivative of the misfit with respect to the model parameters. 

.. math::
    \mathbf{m} = \mathbf{H}^{-1}\mathbf{\gamma}
    
From section 1, we know that our model is, in component form:
    
.. math::
    \Delta g_i = \Gamma Z_{ik}^T \Delta \rho_k
    
For each observation $g_i$, the derivative with respect to the model parameter :math:`\Delta \rho_k` is simply :math:`\Gamma Z_{ik}^T`. This means that :math:`\frac{\partial g_1}{\partial m_1}` is :math:`\Gamma Z_{11}^T`, :math:`\frac{\partial g_1}{\partial m_2}` is :math:`\Gamma Z_{12}^T`, and so on. Continuing this and populating our **G** matrix as defined above, we can see that the **G** matrix is in fact just :math:`\Gamma \mathbf{Z}^T`. This means that our linear model can be written:
  
.. math::
    \Delta \mathbf{g} = \mathbf{G} \Delta \mathbf{\rho}
    
Once **G** is computed, it remains a constant, and the only unknowns are the :math:`\Delta \rho`, which we will refer to as **m**, the model parameters to be estimated, and the model is in fact linear. The data **d** are the observed gravity anomaly values from the Sandwell et al. (2014) global gravity grid and the model parameters to be estimated are the differential densities of each discretized block in the subsurface :math:`\Delta \rho_i`. 


Bayes Theorem and Incorporating Data Errors
--------------------------------------------

To accommodate both data errors and prior constraints on the model parameters, we take a Bayesian approach and calculate the probability that the model parameters take on certain values given the observed data and the prior information. Bayes Theorem states that the probability of the model parameters given the data is proportional to the probability of the data given the model times the probability of the model:

.. math::

    P(\mathbf{m}|\mathbf{d}) \propto P(\mathbf{d}|\mathbf{m})P(\mathbf{m})

In this case, P(**m**) is a prior that we can use to constrain the model parameters to certain values given our existing geophysical knowledge. 

First, however, we will consider the case where we have a constant prior, such that :math:`P(\textbf{m}) = 1` and use Bayes Theorem to incorporate the effect of the data errors on our model. Thus, Bayes Theorem simplifies to:

.. math::
    P(\mathbf{m}|\mathbf{d}) \propto P(\mathbf{d}|\mathbf{m})=P(d_1|\mathbf{m})P(d_2|\mathbf{m})...P(d_n|\mathbf{m})

.. note::
    This makes the key assumption that the data are independent. In the case of gravity, we are dealing with the relative attraction of both adjacent and distal blocks of mass and with gridded gravity data that has undergone some degree of interpolation. Thus, the data are not actually independent. That is, the probability that the gravity is high in one location, say over the Puysegur Ridge, given the model, is not entirely independent of the fact that it is low in another, like the trench, given the model, because that gravity high is the result of a model that simultaneously requires an adjacent gravity low. This can be seen in forward modeling of trench-perpendicular gravity profiles that show trade off in the fit of the gravity between the two sides of the plate boundary, depending on the chosen density distribution. Thus, the probability of two data points for instance is in reality :math:`P(d1,d2)=P(d1)P(d2|d1)` and so on and so forth for n data points. How :math:`P(d2|d1)` would differ from simply :math:`P(d2)` defined by a Gaussian for the point d2 with some known :math:`\sigma_d` is unknown and difficult to quantify given the interdependence of all the data points. If we assume that the perturbation to the probability of d2 by the interdependence is small and that the actual data points we use are distinct gravity measurements, then we can make the simplifying assumption that the data are independent. 

Thus, we can incorporate the data error and priors with the following formulation. Assuming each data point can be represented by a Gaussian distribution with known error, we can take the product of all these probability distributions. The product of a set of exponentials (i.e. of a set of Gaussians) becomes the exponential of the sum of the Gaussian exponents:

.. math::
    P(\mathbf{d}|\mathbf{m})=\frac{1}{\sigma_{d}^{N} (2\pi)^{N/2}} exp[-\frac{1}{2\sigma_{d}^{2}} \sum_{i=1}^{N} (d_i - \mu_i)^2]
    
where :math:`\sigma_d` is the standard deviation of each data point with mean :math:`\mu`.
    
We can see that the exponent term is similar to the definition of the misfit defined previously, so we can define the new misfit as: 
    
.. math::
    F(\mathbf{\mu})=\frac{1}{2\sigma_{d}^{2}} \sum_{i=1}^{N} (d_i - \mu_i)^2
    
That is, we are no minimizing the difference between the known and predicted gravity, given the error in the gravity data. This is for the case when our known $\sigma_d$ is the same for all the data points. For the case where each data point has a different standard deviation, the :math:`\sigma_{d}^{2}` moves inside the sum. 
    
Now, we can write:
    
.. math::
    P(\textbf{d}|\textbf{m})=\frac{1}{\sigma_{d}^{N} (2\pi)^{N/2}} exp(-F(\mathbf{\mu}))
    
and dropping the coefficient, we can make this a proportionality and write this with the simplified Bayes Theorem:
   
.. math::
    P(\textbf{m}|\textbf{d}) \propto P(\textbf{d}|\textbf{m}) \propto exp(-F(\mathbf{\mu}))
    
Now, for the general case where each of the $d_i$ are predicted by the model $g_i(\textbf{m})$, $\mu=g_i(\textbf{m})$ and we have:
    
.. math::
    P(\textbf{m}|\textbf{d}) \propto exp(-F(g_i(\textbf{m})))
    F(\textbf{m}) = \frac{1}{2\sigma_{d}^{2}} \sum_{i=1}^{N} (d_i - g_i(\textbf{m}))^2
    
Thus, minimizing the new misfit F(**m**) is equivalent to maximizing :math:`P(\textbf{m}|\textbf{d})`: finding the best estimate of **m** that maximizes the posterior probability :math:`P(\textbf{m}|\textbf{d})`.
    
We can now define a diagonal and symmetric weight matrix **W** and use a variable substitution to derive the least squares solution that is weighted by the data errors. Let :math:`d_i'=d_i/\sigma_i` and :math:`g_i'=g_i/\sigma_i`. One can construct a matrix **W** such that :math:`\mathbf{d'}=\textbf{W}\textbf{d}` and :math:`\mathbf{g'}=\textbf{W}\textbf{g}`. Starting with our misfit expression, we can make this substitution and arrive at a general expression for the nonlinear LSS in terms of :math:`\mathbf{d'}` and :math:`\mathbf{g'}`.
  
.. math::

    \begin{align}
        F(\textbf{m}) &= \frac{1}{2} \sum_{i=1}^{N} \bigg(\frac{d_i - g_i(\textbf{m})}{\sigma_i}\bigg)^2 \nonumber \\
        &= \frac{1}{2} \sum_{i=1}^{N} \bigg(\frac{d_i}{\sigma_i} - \frac{g_i}{\sigma_i}\bigg)^2 \nonumber \\ 
        &= \frac{1}{2} \sum_{i=1}^{N} (d'_i - g'_i)^2 \nonumber
    \end{align}
    
Now this is just the form of the original "constant :math:`\sigma`" case of the misfit and we can apply the Gauss-Newton Algorithm for the nonlinear least squares solution directly for a case of :math:`\mathbf{G'}`.

.. math::
    F(\textbf{m})=\frac{1}{2} (\mathbf{d'} - \mathbf{g'})^T (\mathbf{d'} - \mathbf{g'})
    \textbf{m} = \textbf{m}_o + (\mathbf{G'}^T \mathbf{G'})^{-1} \mathbf{G'}^T (\mathbf{d'} - \mathbf{g'}(\mathbf{m}_o))
    
Taking the original definition of :math:`\mathbf{d'}` and :math:`\mathbf{g'}` above, we can substitute back in to find the weighted nonlinear least squares solution. We also know that :math:`\mathbf{G'} = \nabla_m \mathbf{g'}`, just by definition of the **G** matrix in general, remembering that **W** does not depend on the model parameters and is thus a matrix of constants in this case. 
  
.. math::
    \begin{align}
        \mathbf{G'} &= \nabla_m \mathbf{g'} \nonumber \\
        &= \nabla_m \mathbf{Wg} \nonumber \\
        &= \mathbf{W} \nabla_m \mathbf{g} \nonumber \\
        &= \mathbf{WG} \nonumber
    \end{align}

.. math::
    \begin{align}
        \textbf{m} &= \textbf{m}_o + ((\textbf{WG})^T (\textbf{WG}))^{-1} (\textbf{WG})^T (\textbf{W}(\textbf{d} - \textbf{g}(\mathbf{m}_o))) \nonumber \\
        &= \textbf{m}_o + (\textbf{G}^T \textbf{W}^T \textbf{WG})^{-1} \textbf{G}^T \textbf{W}^T \textbf{W} (\textbf{d} - \textbf{g}(\mathbf{m}_o)) \nonumber \\
        &= \textbf{m}_o + (\textbf{G}^T \textbf{W}^2 \textbf{G})^{-1} \textbf{G}^T \textbf{W}^2 (\textbf{d} - \textbf{g}(\mathbf{m}_o)) \nonumber
    \end{align}
    
In the linear case, :math:`\mathbf{g}(\mathbf{m}_o) = \mathbf{Gm}_o` and the solution reduces to:
    
.. math::
    \textbf{m} = (\textbf{G}^T \textbf{W}^2 \textbf{G})^{-1} \textbf{G}^T \textbf{W}^2 \textbf{d}
    
In the case of constant priors and with constant data error :math:`\sigma_d` for all points, this equation will actually collapse back to its original unweighted form. However, for the case of non-constant priors, the weights do effect the final estimate of **m**. We will now derive a similar weighted expression that regularizes the solution and stabilizes the inversion. 
    
Tikhonov Regularization
------------------------

Linear least squares, even when using the generalized inverse or the truncated generalized inverse to handle small singular values, is often not sufficient for many inverse problems. Thus, another form of regularization must be applied; a common one is Tikhonov regularization. Zeroth order Tikhonov regularization favors models that are small. In other words, it is identical to applying a prior that is a Gaussian with a mean of zero and thus minimizes the square of the model parameters. First order Tikhonov minimizes the square of the first derivative of the spatially discretized model parameters (i.e. the slopes), which are spatially discretized and thus serves as a flatness criterion. Second order Tikhonov minimizes the square of the second derivative of the model parameters with respect to space (i.e. the gradient) and thus serves as a smoothness criterion. 

Because Tikhonov regularization is essentially applying a prior of either zero, flatness, or smoothness, we can derive the regularized least squares solution by defining a new misfit, and because Tikhonov minimizes the square of the zeroth, first, or second derivative, the misfit equation remains exactly quadratic and the solution linear. Above, we defined a new misfit that incorporated the data error. Likewise, we adjust the misfit equation to reflect the additional minimization of the model parameters or their first or second derivatives. As our model is linear, we will use the fact that :math:`\mathbf{g} = \mathbf{Gm}` and derive the solution using augmented matrices. **L** is either the identity matrix, the first derivative finite difference operator, or the second derivative finite difference operator for zeroth, first, and second order Tikhonov, respectively, and :math:`\alpha` is a constant. Using the weighted form of **d** and **G** defined previously, the new misfit that minimizes both the least squares difference between the true and predicted gravity and the differences between model parameters is:

.. math::
    \begin{align}
        F(\mathbf{m}) &= \frac{1}{2}(\mathbf{Wd} - \mathbf{WG}\mathbf{m})^{T}(\mathbf{Wd}-\mathbf{WG}\mathbf{m}) + \alpha^2 (\mathbf{L}\mathbf{m})^T (\mathbf{L}\mathbf{m}) \\
        F(\mathbf{m}) &= \frac{1}{2} \bigg( (\mathbf{Wd} - \mathbf{WGm})^T (\mathbf{Wd} - \mathbf{WGm}) + 2\alpha^2 (\mathbf{Lm})^T (\mathbf{Lm})\bigg)
    \end{align}

Because alpha is just a constant, we can just absorb the 2 into alpha and write this equation with augmented matrices.

.. math::
    \begin{align}
        F(\textbf{m}) &= \frac{1}{2}\bigg( \begin{bmatrix} \mathbf{Wd} \\ \mathbf{0} \end{bmatrix}
            - \begin{bmatrix} \mathbf{WG} \\ \alpha \mathbf{L} \end{bmatrix} \mathbf{m} \bigg)^T 
            \bigg(\begin{bmatrix} \mathbf{Wd} \\ \mathbf{0} \end{bmatrix}
            - \begin{bmatrix} \mathbf{WG} \\ \alpha \mathbf{L} \end{bmatrix} \mathbf{m}\bigg) \\
        F(\mathbf{m}) &= \frac{1}{2} (\mathbf{d}_{aug} - \mathbf{G}_{aug} \mathbf{m})^T (\mathbf{d}_{aug} - \mathbf{G}_{aug} \mathbf{m})
    \end{align}

Now this is just the original form of the misfit equation and we can use the Gauss-Newton solution for the non-linear least squares with these new variables:

.. math::
    \begin{align}
        \mathbf{m} &= \mathbf{m}_o + (\mathbf{G}_{aug}^{T} \mathbf{G}_{aug})^{-1} \mathbf{G}_{aug}^{T} (\mathbf{d}_{aug} - \mathbf{G}_{aug} \mathbf{m}_o) \nonumber \\
        &= \mathbf{m}_o + (\mathbf{G}_{aug}^{T} \mathbf{G}_{aug})^{-1} \mathbf{G}_{aug}^{T} \mathbf{d}_{aug} - (\mathbf{G}_{aug}^{T} \mathbf{G}_{aug})^{-1} \mathbf{G}_{aug}^{T} \mathbf{G}_{aug} \mathbf{m}_o \nonumber \\
        &= \mathbf{m}_o + (\mathbf{G}_{aug}^{T} \mathbf{G}_{aug})^{-1} \mathbf{G}_{aug}^{T} \mathbf{d}_{aug} - \mathbf{m}_o \nonumber \\
        &= (\mathbf{G}_{aug}^{T} \mathbf{G}_{aug})^{-1} \mathbf{G}_{aug}^{T} \mathbf{d}_{aug} \nonumber
    \end{align}

Because our problem is linear, the nonlinear LS solution reduces back to the linear form as expected. Now substituting back in our definitions of :math:`\mathbf{G}_{aug}` and :math:`\mathbf{d}_{aug}`, we can arrive at the final weighted, Tikhonov regularized solution. 

.. math::
    \begin{align} 
        \mathbf{G}_{aug}^{T} \mathbf{G}_{aug} &= \begin{bmatrix} (\mathbf{WG})^{T} \alpha \mathbf{L}^{T} \end{bmatrix} \begin{bmatrix} \mathbf{WG} \\ \alpha \mathbf{L} \end{bmatrix} = (\mathbf{WG})^{T} \mathbf{WG} + \alpha^2 \mathbf{L}^{T} \mathbf{L} \nonumber \\
        &= \mathbf{G}^T \mathbf{W}^T \mathbf{WG} + \alpha^2 \mathbf{L}^T \mathbf{L} \nonumber \\
        &= \mathbf{G}^T \mathbf{W}^2 \mathbf{G} + \alpha^2 \mathbf{L}^T \mathbf{L} \nonumber \\
        \mathbf{G}_{aug}^{T} \mathbf{d}_{aug} &= \begin{bmatrix} (\mathbf{WG})^{T} \alpha \mathbf{L}^{T} \end{bmatrix} \begin{bmatrix} \mathbf{Wd} \\ \mathbf{0} \end{bmatrix} = (\mathbf{WG})^{T} \mathbf{Wd} \nonumber \\
        &= \mathbf{G}^T \mathbf{W}^2 \mathbf{d} \nonumber \\
        & \nonumber \\
        \mathbf{m} &= (\mathbf{G}^{T} \mathbf{W}^2 \mathbf{G} + \alpha^2 \mathbf{L}^{T} \mathbf{L})^{-1} \mathbf{G}^{T} \mathbf{W}^2 \mathbf{d} \nonumber
    \end{align}

Knowing that :math:`\mathbf{m} = \mathbf{H}^{-1} \mathbf{\gamma}` we can also define the weighted Tikhonov regularized Hessian and gradient: 

.. math::
    \mathbf{H} = (\mathbf{G}^{T} \mathbf{W}^2 \mathbf{G} + \alpha^2 \mathbf{L}^{T} \mathbf{L})
    \mathbf{\gamma} = (\mathbf{G}^{T} \mathbf{W}^2 \mathbf{d})

For three dimensional models, first and second order Tikhonov regularization are implemented using a finite-difference approximation to the Laplacian operator of the appropriate dimensionality, and is thus the first or second derivative in the x-direction plus the first or second derivative in the y-direction plus the first or second derivative in the z-direction. Because the discretization of the grid will be different in the x, y, and z directions, we need to apply three different regularizations, with associated constants :math:`\alpha` for the x-direction, :math:`\beta` for the y-direction, and :math:`\zeta` for the z-direction. Following the derivation of the weighted regularized solution above, we can see by inspection that the appropriate weighted Tikhonov regularized solution for the 3-dimensional case is:

.. math::
    \mathbf{m}=(\mathbf{G}^{T}\mathbf{W}^2\mathbf{G} + \alpha^2\mathbf{L_x}^{T}\mathbf{L_x} + \beta^2\mathbf{L_y}^{T}\mathbf{L_y} + \zeta^2\mathbf{L_z}^{T}\mathbf{L_z})^{-1} (\mathbf{G}^{T}\mathbf{W}^2\mathbf{d})

First Order Tikhonov Regularization with Variable Grid Spacing
---------------------------------------------------------------

Determining the structure of the operator **L** for either first or second order Tikhonov is dependent on the grid structure. In the general, one-dimensional case, the finite difference approximation to the first derivative is:

.. math::
    \frac{\partial m_k}{\partial x} = \frac{1}{\Delta x}(-m_{k} + m_{k+1})

This equation can be represented in the form of a tri-diagonal matrix operator, **L**, acting on a vector of the spatially discretized model parameters. **L1** is an M-1 x M matrix.  

.. math::
    \begin{align}
        \frac{\partial m_k}{\partial x} &= \frac{1}{\Delta x} 
        \begin{bmatrix}
        -1 & 1 & 0 & \cdots & 0 \\
        0 & -1 & 1 & \cdots & 0 \\
        \vdots & \ddots & \ddots & \ddots & 0 \\
        0 & \cdots & 0 & -1 & 1 
        \end{bmatrix} 
        \begin{bmatrix}
        m_1 \\ m_2 \\ \vdots \\ m_k
        \end{bmatrix} \nonumber \\
        &= \frac{1}{\Delta x} \mathbf{L1m} \nonumber
    \end{align} 

Because the discretization is variable even within the x, y, and z directions, the **L** matrices are unique for each of those directions and :math:`\Delta x, \Delta y, ` or :math:`\Delta z` is different for each pair of model parameters being regularized, so the :math:`1/\Delta x` term must be brought inside the **L** matrix. The discretization is such that the model parameters are numbered first in the z-direction, then in the y-direction, then in the x-direction. This means that in the x direction we have to regularize parameters that are spaced length(z)*length(y) apart. In other words, model parameter 1 and model parameter 2 are not adjacent to one another in the x-direction because they lie at the same x-value. This means that we place the second diagonal of the matrix at $1+length(z)*length(y)$, and the **L1** matrix for the x-direction takes the form:

.. math::
    \begin{align}\mathbf{L1}_x = \frac{\partial m_k}{\partial x} &=  
        \begin{bmatrix}
        \frac{-1}{\Delta x_1} & 0 & \cdots & \frac{1}{\Delta x_1} & 0 & \cdots & & & & 0 \\
        0 & \frac{-1}{\Delta x_1} & 0 & \cdots & \frac{1}{\Delta x_1} & 0 & \cdots & & & 0 \\
        \vdots & \ddots & \ddots & \ddots & & \ddots & \ddots & & & \vdots \\
        0 & \cdots & 0 & \frac{-1}{\Delta x_p} & 0 & \cdots & \frac{1}{\Delta x_p} & 0 & \cdots & 0 \\
        0 & & \cdots & 0 & \frac{-1}{\Delta x_{p+1}} & 0 & \cdots & \frac{1}{\Delta x_{p+1}} & 0 & \cdots \\
        \vdots & & & & \ddots & \ddots & \ddots & & \ddots & \ddots\\
        0 & & & & \cdots & 0 & \frac{-1}{\Delta x_{k-l}} & 0 & \cdots & \frac{1}{\Delta x_{k-l}}
        \end{bmatrix}
        \begin{bmatrix}
        m_1 \\ m_2 \\ \vdots \\ m_p \\ m_{p+1} \\ \vdots \\ m_k
        \end{bmatrix} \nonumber \\ \nonumber
    \end{align}

where the length of the :math:`\Delta x` vector that is applied to the **L1** matrix is :math:`length(x)-1`. Every :math:`length(z)*length(y)` rows are divided by :math:`\Delta x_i` because for each of these sets of rows we are regularizing the difference between model parameters for all prisms in y and in z that have a :math:`\Delta x` of :math:`x_{p+1} - x_p`. After that set of rows, we move to the next pair of model parameters in the x-direction, which have a new :math:`\Delta x` value. The :math:`\mathbf{L1}_x` matrix is an :math:`(M-length(z)*length(y))` x :math:`M` matrix. 

We can likewise create a similar matrix for the y-direction. However, because we move along the y-direction second in the grid, we place the second non-zero diagonal at :math:`1+length(z)`, and in order to avoid regularizing nonadjacent parameters at the edges of the grid space (i.e. at the beginning and end of y where the numbering jumps to the next x-value), we also have to skip (zero out) a set of :math:`length(z)` rows after each set of :math:`length(z)*(length(y)-1)` rows. Each non-null set of rows consists of :math:`length(y)-1` smaller sets of rows, each of :math:`length(z)`, which are divided by their corresponding :math:`\Delta y`. After the null rows, the sequence is repeated. The :math:`\mathbf{L1}_y` matrix is an :math:`(M-length(z))` x :math:`M` matrix.

.. math::
    \begin{align}\mathbf{L1}_y = \frac{\partial m_k}{\partial y} &=  
        \begin{bmatrix}
        \frac{-1}{\Delta y_1} & 0 & \cdots & \frac{1}{\Delta y_1} & 0 & \cdots & & & & 0 \\
        0 & \frac{-1}{\Delta y_1} & 0 & \cdots & \frac{1}{\Delta y_1} & 0 & \cdots & & & 0 \\
        \vdots & \ddots & \ddots & \ddots & & \ddots & \ddots & & & \vdots \\
        0 & \cdots & 0 & \frac{-1}{\Delta y_p} & 0 & \cdots & \frac{1}{\Delta y_p} & 0 & \cdots & 0 \\
        0 & & \cdots & 0 & \frac{-1}{\Delta y_{p+1}} & 0 & \cdots & \frac{1}{\Delta y_{p+1}} & 0 & \cdots \\
        \vdots & & & & \ddots & \ddots & \ddots & & \ddots & \ddots\\
        0 & & & & \cdots & 0 & \frac{-1}{\Delta y_{k-l}} & 0 & \cdots & \frac{1}{\Delta y_{k-l}}
        \end{bmatrix}
        \begin{bmatrix}
        m_1 \\ m_2 \\ \vdots \\ m_p \\ m_{p+1} \\ \vdots \\ m_k
        \end{bmatrix} \nonumber
    \end{align}

The **L** matrix for the z-direction takes the form of the standard first derivative finite difference operator, except that every :math:`length(z)-1` rows are skipped in order to avoid regularizing parameters located at the top and bottom of the grid and not adjacent to one another. Each row in a set of :math:`length(z)-1` rows is divided by its corresponding :math:`\Delta z` value. Then a row is skipped (zeroed out) and the sequence of :math:`\Delta z` is repeated. This makes the :math:`\mathbf{L1}_z` matrix and :math:`(M-1)` x :math:`M` matrix.

.. math::
    \begin{align}\mathbf{L1}_z = \frac{\partial m_k}{\partial z} &=  
        \begin{bmatrix}
        \frac{-1}{\Delta z_1} & \frac{1}{\Delta z_1} & 0 & \cdots & & & \cdots & 0 \\
        0 & \frac{-1}{\Delta z_1} & \frac{1}{\Delta z_1} & 0 & \cdots & & \cdots & 0 \\
        \vdots & \ddots & \ddots & \ddots & & & & \vdots \\
        0 & \cdots & 0 & 0 & 0 & 0 & \cdots & 0 \\
        0 & & \cdots & 0 & \frac{-1}{\Delta z_{p1}} & \frac{1}{\Delta z_{p}} & 0 & \vdots \\
        \vdots & & & & \ddots & \ddots & \ddots & & \\
        0 & \cdots & & & \cdots & 0 & \frac{-1}{\Delta z_{k-l}} & \frac{1}{\Delta z_{k-l}} 
        \end{bmatrix}
        \begin{bmatrix}
        m_1 \\ m_2 \\ \vdots \\ m_p \\ m_{p+1} \\ \vdots \\ m_k
        \end{bmatrix} \nonumber
    \end{align}
    
Even though :math:`\mathbf{L1}_x`, :math:`\mathbf{L1}_y`, and :math:`\mathbf{L1}_z` are each different sizes, the product of their transpose with themselves still results in an M x M matrix that is compatible with the least squares solution **m**. To be more explicit about the first order Tikhonov solution for variable grid refinement, we can write the solution as:

.. math::
    \mathbf{m} = (\mathbf{G}^T \mathbf{W}^2 \mathbf{G} + \alpha^2 \mathbf{L1}_x^T \mathbf{L1}_x + \beta^2 \mathbf{L1}_y^T \mathbf{L1}_y + \zeta^2 \mathbf{L1}_z^T \mathbf{L1}_z)^{-1} (\mathbf{G}^T \mathbf{W}^2 \mathbf{d})

Second Order Tikhonov Regularization for Variable Grid Spacing
---------------------------------------------------------------

In the general, one-dimensional case, the finite difference approximation to the second derivative is:

.. math::
    \frac{\partial^2 m_k}{\partial x^2} = \frac{1}{\Delta x^2}(m_{k-1} - 2m_k + m_{k+1})

This equation can be represented in the form of a tri-diagonal matrix operator, **L**, acting on a vector of the spatially discretized model parameters. **L2** is an :math:`(M-2)` x :math:`M` matrix. 

.. math::
    \begin{align}
        \frac{\partial^2 m_k}{\partial x^2} &= \frac{1}{\Delta x^2} 
        \begin{bmatrix}
        1 & -2 & 1 & 0 & \cdots & 0 \\
        0 & 1 & -2 & 1 & \cdots & 0 \\
        \vdots & \ddots & \ddots & \ddots & \ddots & 0 \\
        0 & \cdots & 0 & 1 & -2 & 1 
        \end{bmatrix}
        \begin{bmatrix}
        m_1 \\ m_2 \\ \vdots \\ m_k
        \end{bmatrix} \nonumber \\
        &= \frac{1}{\Delta x^2} \mathbf{L2m} \nonumber
    \end{align}

However, because the discretization is variable, we cannot use the straightforward formulation above and we have to consider that for each row, there are two different :math:`\Delta x` involved. This being said, the simplest way to calculate the appropriate :math:`\mathbf{L2}` matrix is as a function of the :math:`\mathbf{L1}` matrix. The first derivative finite difference approximation is the approximation to the slope of a curve at the midpoint between two points :math:`m_{k-1}` and :math:`m_k`, separated by a distance :math:`\Delta x_i`. We have already accounted for the variable :math:`\Delta x_i` in the :math:`\mathbf{L1}` matrix. The second derivative approximation then is the slope of the slope, i.e. the curvature, so we can do the same finite difference approximation again, by calculating the difference in the slopes between model parameters over the distance :math:`\delta x_i` - the difference between the points at which the slopes were evaluated, which are the midpoints for each :math:`\Delta x_i` segment. Thus, the second derivative finite difference approximation can be written as a function of the first derivative finite difference approximation as:

.. math::
    \frac{\partial^2 m_k}{\partial x^2} = \frac{1}{\delta x_j}\bigg(\frac{\partial m_{k+1}}{\partial x_{i+1}} - \frac{\partial m_k}{\partial x_i}\bigg)

    \frac{\partial^2 m_k}{\partial x^2} = \frac{1}{\delta x_j}\bigg(\frac{-m_k + m_{k+1}}{\Delta x_{i+1}} - \bigg(\frac{-m_{k-1} + m_k}{\Delta x_i}\bigg)\bigg)

We can see that if all the :math:`\Delta x_i`, and likewise :math:`\delta x_j`, are equal, then this equation becomes the standard finite difference approximation to the second derivative above, but with all the :math:`\Delta x_i` being different, we can leave the equation in this form and see that it is in fact the difference between different subsets of the :math:`\mathbf{L1}` matrices, normalized by :math:`\delta x = \frac{1}{2}(\Delta x_{i+1} + \Delta x_i)`. This method can be applied separately to the x, y, and z directions. The sets of rows that are subtracted from one another depend on the ordering of the model parameters and which direction is being regularized. 

Like with first order Tikhonov, the **L** matrices for second order Tikhonov are unique for each direction, x, y, and z. Following the reasoning for the :math:`\mathbf{L1}_x` matrix, the second diagonal is at :math:`1+length(z)*length(y)` and now the third diagonal is at :math:`1+2*length(z)*length(y)`. Each block of :math:`length(z)*length(y)` rows in :math:`\mathbf{L1}_x` is subtracted from the subsequent :math:`length(z)*length(y)` block of rows in :math:`\mathbf{L1}_x`, and that new block is then multiplied by one over the respective distance :math:`\delta x` for the parameters being regularized by those rows. This ultimately results in :math:`\mathbf{L2}_x` being an :math:`(M - 2*length(z)*length(y))` x :math:`M` matrix.

.. math::
    \begin{align}\mathbf{L2}_x =
        \begin{bmatrix}
            \frac{1}{\delta x_1} \\
            \frac{1}{\delta x_1} \\
            \vdots \\
            \frac{1}{\delta x_p} \\
            \frac{1}{\delta x_p} \\
            \vdots \\ 
        \end{bmatrix} 
        .*
        \begin{bmatrix} 
        \frac{1}{\Delta x_1} & \cdots & \frac{-1}{\Delta x_2}-\frac{1}{\Delta x_1} & \cdots & \frac{1}{\Delta x_2} & \cdots & & & & 0 \\
        0 & \frac{1}{\Delta x_1} & \cdots & \frac{-1}{\Delta x_2}-\frac{1}{\Delta x_1} & \cdots & \frac{1}{\Delta x_2} & \cdots & & & 0 \\
        \vdots & & \ddots & & \ddots & & \ddots & & & \vdots\\
        0 & \cdots & 0 & \frac{1}{\Delta x_p} & \cdots & \frac{-1}{\Delta x_{p+1}}-\frac{1}{\Delta x_p} & \cdots & \frac{1}{\Delta x_{p+1}} & & \vdots \\
        0 & \cdots & & 0 & \frac{1}{\Delta x_{p}} & \cdots & \frac{-1}{\Delta x_{p+1}}-\frac{1}{\Delta x_{p}} & \cdots & \frac{1}{\Delta x_{p+1}} & 0\\
        \vdots & & & & \ddots & \ddots & \ddots & & \ddots & \ddots\\
        \end{bmatrix} \nonumber 
    \end{align}

The computation of the :math:`\mathbf{L2_y}` matrix is similar to that of the :math:`\mathbf{L2_x}` matrix in that is it determined by the :math:`\mathbf{L1_y}` matrix. However, now the first diagonal is at :math:`1 + length(z)` and the second diagonal is at :math:`1+ 2*length(z)`. Each block of :math:`length(z)` rows in :math:`\mathbf{L1_y}` is subtracted from the subsequent block of :math:`length(z)` rows of :math:`\mathbf{L1_y}`. To make sure we don't regularize parameters at opposite edges of the domain, we skip (zero out) a set of :math:`2*length(z)` rows every :math:`length(z)*(length(y)-2)` rows. Each of the resulting non-zero blocks of :math:`length(z)` rows is divided element wise by its corresponding :math:`\delta y` value, which result in an :math:`(M-2*length(z))` x :math:`M` matrix.

For :math:`\mathbf{L2_z}`, we do the same thing, except now each row is subtracted from the subsequent row, since the grid loops over the z-direction first. The second diagonal is in the second column and the third diagonal is in the third column. To avoid regularizing the top of the domain with the bottom, we zero two rows every :math:`length(z)-2` rows. Each of the resulting non-zero blocks of :math:`length(z)-2` rows is then divided by its corresponding :math:`\delta z` value, and the resulting :math:`\mathbf{L2_z}` matrix is an :math:`(M-2)` x :math:`M` matrix. 

Boundary Conditions and Variable Order Tikhonov Regularization
---------------------------------------------------------------

To ensure mathematical stability, we need to enforce an infinite boundary condition. By not imposing flatness in the far-field, we get abrupt density changes at the edge of the model domain and the classical gravity edge effect as a result. This does not make physical sense given the actual gravity data, so we suppress them with the infinite edge boundary condition, which lets the gravity smoothly continue off the edges of the model area. This is easily done in the forward model by making sure the edge prisms in the x and y directions are all sufficiently large that they are essentially infinite. The edge prisms of the gridded domain should each be at least long enough to extend the boundaries of the domain beyond the area for which we are calculating the gravity. However, to ensure that this boundary condition is applied in the actual inversion itself, there needs to be a flatness criterion between the edge parameters and their adjacent parameters. This can be done using first order Tikhonov regularization with a very large alpha value that minimizes the difference between the edge parameters and their adjacent values so much that the predicted density values of those two parameters are the same. Doing this requires that we apply different orders and/or strengths of Tikhonov regularization to different sets of model parameters simultaneously. In doing so, we can apply second order smoothing to the main model area while maintaining the flatness on the boundary. 

Like before, variable order Tikhonov can be achieved by redefining the misfit equation:

.. math::
    \begin{align} 
        F(\mathbf{m}) &= \frac{1}{2}(\mathbf{Wd} - \mathbf{WG}\mathbf{m})^{T}(\mathbf{Wd}-\mathbf{WG}\mathbf{m}) \nonumber \\
        &+ \alpha^2 (\mathbf{L_x}\mathbf{m})^T (\mathbf{L_x}\mathbf{m}) + \beta^2(\mathbf{L_y}\mathbf{m})^T (\mathbf{L_y}\mathbf{m}) + \zeta^2(\mathbf{L_z}\mathbf{m})^T (\mathbf{L_z}\mathbf{m}) \nonumber \\
        &+ b^2(\mathbf{B_x}\mathbf{m})^T (\mathbf{B_x}\mathbf{m}) + b^2(\mathbf{B_y}\mathbf{m})^T (\mathbf{B_y}\mathbf{m}) \nonumber
    \end{align}

where b is the weight of the first order Tikhonov regularization applied to the boundary condition. :math:`b=1e8` is usually sufficient to flatten the edges of the model. :math:`\mathbf{B_x}` and :math:`\mathbf{B_y}` are the regularization matrices that apply the boundary conditions in the x and y directions, respectively. To implement the variable order Tikhonov, the structure of the :math:`\mathbf{L_x}` and :math:`\mathbf{L_y}` matrices, either for first or second order, must be adjusted accordingly. The :math:`\mathbf{B_x}` and :math:`\mathbf{B_y}` matrices are constructed from the :math:`\mathbf{L1_x}` and :math:`\mathbf{L1_y}` matrices, respectively. They are essentially equal, except that only the rows that regularize parameters including edge values remain non-null; all rows that regularize parameters within the interior of the domain are zeroed out. The :math:`\mathbf{L_x}` and :math:`\mathbf{L_y}` matrices are then made to be the opposite. All rows that regularize parameters that include edges are zeroed out, so that only the interior parameters are regularized by these matrices. We can then easily apply different regularization weights to the edge versus the interior parameters to enforce the boundary condition. 

The five Tikhonov regularization terms can be combined into a single :math:`\mathbf{L}` matrix such that 

.. math::
    F(\mathbf{m}) = \frac{1}{2}(\mathbf{Wd} - \mathbf{WGm})^T (\mathbf{Wd} - \mathbf{WGm}) + \mathbf{m}^T \mathbf{L} \mathbf{m} 

    \mathbf{m} = (\mathbf{G}^T \mathbf{W}^2 \mathbf{G} + \mathbf{L})^{-1}(\mathbf{G}^T\mathbf{W}^2\mathbf{d})

where :math:`\mathbf{L} = \alpha^2 \mathbf{L_x}^T\mathbf{L_x} + \beta^2 \mathbf{L_y}^T\mathbf{L_y} + \zeta^2 \mathbf{L_z}^T\mathbf{L_z} + b^2 \mathbf{B_x}^T \mathbf{B_x} + b^2 \mathbf{B_y}^T \mathbf{B_y}`. We can do this because we can rearrange the regularization terms in the misfit equation in the following way.

.. math::
    \begin{align}
        F(\mathbf{m}) &= \cdots + \alpha^2 (\mathbf{L_x}\mathbf{m})^T (\mathbf{L_x}\mathbf{m}) + \beta^2(\mathbf{L_y}\mathbf{m})^T (\mathbf{L_y}\mathbf{m}) + \zeta^2(\mathbf{L_z}\mathbf{m})^T (\mathbf{L_z}\mathbf{m}) \nonumber \\
        &+ b^2(\mathbf{B_x}\mathbf{m})^T (\mathbf{B_x}\mathbf{m}) + b^2(\mathbf{B_y}\mathbf{m})^T (\mathbf{B_y}\mathbf{m}) \nonumber \\
        \\
        &= \cdots + \alpha^2 (\mathbf{m}^T \mathbf{L_x}^T)(\mathbf{L_x m}) + \beta^2 (\mathbf{m}^T \mathbf{L_y}^T)(\mathbf{L_y m}) + \zeta^2(\mathbf{m}^T \mathbf{L_z}^T)(\mathbf{L_z m}) \nonumber \\
        &+ b^2(\mathbf{m}^T \mathbf{B_x}^T)(\mathbf{B_x m}) + b^2(\mathbf{m}^T \mathbf{B_y}^T)(\mathbf{B_y m}) \nonumber \\
        \\
        &= \cdots + (\alpha^2 \mathbf{m}^T\mathbf{L_x}^T\mathbf{L_x} + \beta^2 \mathbf{m}^T \mathbf{L_y}^T \mathbf{L_y} + \zeta^2 \mathbf{m}^T \mathbf{L_z}^T \mathbf{L_z} + b^2 \mathbf{m}^T \mathbf{B_x}^T \mathbf{B_x} + b^2 \mathbf{m}^T \mathbf{B_y}^T \mathbf{B_y})\mathbf{m} \nonumber \\
        \\
        &= \cdots + \mathbf{m}^T(\alpha^2 \mathbf{L_x}^T\mathbf{L_x} + \beta^2\mathbf{L_y}^T \mathbf{L_y} + \zeta^2 \mathbf{L_z}^T \mathbf{L_z} + b^2 \mathbf{B_x}^T \mathbf{B_x} + b^2 \mathbf{B_y}^T \mathbf{B_y})\mathbf{m} \nonumber 
    \end{align}



    
