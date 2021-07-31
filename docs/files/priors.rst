=======
Priors
=======

Incorporating Non-Constant Priors
----------------------------------

Having stabilized the solution with Tikhonov Regularization, we can now incorporate further constraints in the form of non-constant priors. Tikhonov regularization is itself a type of prior and acts as a relative equality constraint, meaning it smooths the model. Using seismic data, we can also apply absolute equality constraints to certain model parameters such that they take on values within a specified range, consistent with our knowledge of the region's geology and structure. Different areas of the modeled subsurface (i.e. different groups of model parameters) can have different priors depending on the extent of our knowledge about that area. We will use Gaussian priors, where the values of the chosen model parameters can vary within some specified range (the standard deviation) about some expected value (the mean). Gaussians are the simplest option to use and we don't have a more valid reason to use a different probability density function. Even though some of the priors will be the same for large sections of the model, we will assume that the priors on each model parameter are independent, such that:

.. math::
    P(\textbf{m})=P(m_1)P(m_2)P(m_3)...P(m_M)

We use the example of a buried sphere for the first set of synthetic tests. Lets say we know the sphere has a differential density of :math:`\Delta \rho_{sphere}`. So, the priors on all the parameters (prisms) inside a set of points that define a sphere are such that these :math:`\Delta \rho` values are expected to be around some density value :math:`\mu_{sphere}` and can vary between :math:`\pm \sigma_{sphere}`. Let us assume for now that the priors on all other model parameters (i.e. the prisms outside the spherical region) are uniform (i.e. that the probability of these model parameters being any value is 1). Thus, the prior on each of the sphere parameters will be:

.. math::
    P(m_k)=\frac{1}{\sigma_{sphere} \sqrt{2 \pi}} exp (-\frac{1}{2\sigma_{sphere}^{2}} (m_k - \mu_{sphere})^2)

In the same way as we did above with incorporating the data error, we can use the fact that the product of a set of Gaussians is the exponential of the sum of the Gaussian exponents and get that for the case of the :math:`\mu_{sphere}` prior:

.. math::
    P(\mathbf{m}_{sphere})=\frac{1}{\sigma_{sphere}^{N_{sphere}} (2 \pi)^{N_{sphere}/2}} exp(-\frac{1}{2\sigma_{sphere}^{2}} \sum_{k=1}^{N_{sphere}} (m_k - \mu_{sphere})^2)   

where :math:`N_{sphere}` is the number of model parameters inside the spherical region. Dropping the exponential to make this a proportionality as before and returning to Bayes Theorem, we can write:

.. math::
    \begin{align}
        P(\textbf{m}|\textbf{d}) &\propto e^{-F(\textbf{m})} e^{-\frac{1}{2\sigma_{sphere}^{2}} \sum_{k=1}^{N_{sphere}} (m_k - \mu_{sphere})^2} \nonumber \\
        &\propto e^{-(F(\textbf{m})+\frac{1}{2\sigma_{sphere}^{2}}\sum_{k=1}^{N_{sphere}}(m_k - \mu_{sphere})^2)} \nonumber
    \end{align}

Thus, the new misfit that now incorporates both data errors and non-constant priors is as follows, with :math:`F(\mathbf{m})` being the original Tikhonov regularized misfit equation.

.. math::
    F_{new}(\textbf{m})=F(\textbf{m}) + \frac{1}{2\sigma_{sphere}^{2}}\sum_{k=1}^{N_{sphere}}(m_k - \mu_{sphere})^2

Likewise for different priors applied to different sets of model parameters, one can apply the same method above and the resulting misfit will be the sum of the original misfit and each of the exponential components that result from the different priors. 

There are a couple of ways to go about deriving a new equation for the least squares solution from this misfit. Recall that the least squares misfit can also be written :math:`\textbf{m} = \textbf{m}_o + \textbf{H}^{-1}\mathbf{\gamma}`, where :math:`\textbf{H}_{kk}=\frac{\partial^{2}F}{\partial m_k \partial m_{k'}}` and :math:`\mathbf{\gamma}_k=\frac{\partial F}{\partial m_k}`, and **m** is in terms of the original F(**m**). To update this expression, with our new :math:`F_{new}(\textbf{m})`, one can simply update the gradient and the Hessian term by term by taking the first and second derivative of the new misfit, respectively, such that :math:`\textbf{m} = \textbf{m}_o + \textbf{H}_{new}^{-1} \mathbf{\gamma}_{new}`. However, for applying many different priors to different sets of model parameters, where we have to take care on the bookkeeping so that we are constraining the right parameters in the right places, it is easiest to derive the constrained solution using augmented matrices as we did with the Tikhonov regularization above. In this way, we can write the new misfit as:

.. math::
    \begin{align}
        F(\mathbf{m}) &= \frac{1}{2\sigma_d^2}(\mathbf{d}-\mathbf{Gm})^T(\mathbf{d}-\mathbf{Gm}) + \alpha^2 (\mathbf{Lm})^T (\mathbf{Lm}) + \frac{1}{2\sigma_p^2}(\bm{\mu}_p-\mathbf{m})^T (\bm{\mu}_p-\mathbf{m}) \nonumber \\
        F(\mathbf{m}) &= \frac{1}{2}\sum_{i=1}^N \bigg(\frac{d_i - (\mathbf{Gm})_i}{\sigma_{d_i}}\bigg)^2 + \alpha^2 \sum_{j=1}^M (\mathbf{Lm})_j^2 + \frac{1}{2}\sum_{j=1}^M \bigg(\frac{\mu_{p_j} - m_j}{\sigma_{p_j}}\bigg)^2 \label{misfit_3dTik_pr}
    \end{align}

where for simplicity, we have included only one Tikhonov term here (there are in fact five). As before, let there be a diagonal weight matrix :math:`\mathbf{W}=diag(1/\mathbf{\sigma_d})` with the data errors on the diagonal, and similarly, let there be a diagonal matrix :math:`\mathbf{S} = diag(1/\mathbf{\sigma_p})` with the chosen standard deviations of the priors for each model parameter on the diagonal. As we saw previously the weighted misfit can now be written: 

.. math::
    \begin{align}
        F(\mathbf{m}) &= \frac{1}{2}(\mathbf{Wd}-\mathbf{WGm})^T(\mathbf{Wd}-\mathbf{WGm}) + \alpha^2 (\mathbf{Lm})^T (\mathbf{Lm}) + \frac{1}{2}(\mathbf{S}\bm{\mu}-\mathbf{Sm})^T (\mathbf{S}\bm{\mu}-\mathbf{Sm}) \nonumber \\
        &= \frac{1}{2}\Bigg(\begin{bmatrix} \mathbf{Wd}\\ \mathbf{0}\\ \mathbf{S}\bm{\mu} \end{bmatrix} - \begin{bmatrix} \mathbf{WG}\\ \alpha\mathbf{L}\\ \mathbf{S} \end{bmatrix} \mathbf{m}\Bigg)^T \nonumber \Bigg(\begin{bmatrix} \mathbf{Wd}\\ \mathbf{0}\\ \mathbf{S}\bm{\mu} \end{bmatrix} - \begin{bmatrix} \mathbf{WG}\\ \alpha \mathbf{L}\\ \mathbf{S} \end{bmatrix} \mathbf{m}\Bigg) \nonumber \\
        &= \frac{1}{2} (\mathbf{d}_{aug} - \mathbf{G}_{aug}\mathbf{m})^T (\mathbf{d}_{aug} - \mathbf{G}_{aug}\mathbf{m}) \nonumber
    \end{align}

Now, as before, this is just the original form of our misfit equation, so we can write the Gauss-Newton solution for the least squares problem and substitute back in our definitions of :math:`\mathbf{G}_{aug}` and :math:`\mathbf{d}_{aug}` to arrive at the final Tikhonov regularized solution with prior constraints. 

.. math::
    \begin{align}
        \mathbf{m} &= \mathbf{m}_o + (\mathbf{G}_{aug}^T \mathbf{G}_{aug})^{-1} \mathbf{G}_{aug}^T (\mathbf{d}_{aug} - \mathbf{G}_{aug}\mathbf{m}_o) \nonumber \\
        &= \mathbf{m}_o + (\mathbf{G}_{aug}^T \mathbf{G}_{aug})^{-1} \mathbf{G}_{aug}^T \mathbf{d}_{aug} - (\mathbf{G}_{aug}^T \mathbf{G}_{aug})^{-1} \mathbf{G}_{aug}^T \mathbf{G}_{aug}\mathbf{m}_o \nonumber \\
        &= (\mathbf{G}_{aug}^T \mathbf{G}_{aug})^{-1} \mathbf{G}_{aug}^T \mathbf{d}_{aug} \nonumber \\ 
        \nonumber \\
        \mathbf{G}_{aug}^T \mathbf{G}_{aug} &= \begin{bmatrix} \mathbf{WG}\\ \alpha \mathbf{L}\\ \mathbf{S} \end{bmatrix}^T \begin{bmatrix} \mathbf{WG}\\ \alpha \mathbf{L} \mathbf{S} \end{bmatrix} = \begin{bmatrix} (\mathbf{WG})^T \alpha \mathbf{L}^T \mathbf{S}^T \end{bmatrix} \begin{bmatrix} \mathbf{WG}\\ \alpha \mathbf{L}\\ \mathbf{S} \end{bmatrix} \nonumber \\
        &= (\mathbf{WG})^T \mathbf{WG} + \alpha^2 \mathbf{L}^T\mathbf{L} + \mathbf{S}^T\mathbf{S} \nonumber \\
        &= \mathbf{G}^T\mathbf{W}^2\mathbf{G} + \alpha^2\mathbf{L}^T\mathbf{L} + \mathbf{S}^2 \nonumber \\
        \mathbf{G}_{aug}^T \mathbf{d}_{aug} &= \begin{bmatrix} \mathbf{WG}\\ \alpha \mathbf{L}\\ \mathbf{S} \end{bmatrix}^T \begin{bmatrix} \mathbf{Wd}\\ \mathbf{0}\\ \mathbf{S}\bm{\mu} \end{bmatrix} = \begin{bmatrix} (\mathbf{WG})^T \ \alpha \mathbf{L}^T \mathbf{S}^T \end{bmatrix} \begin{bmatrix} \mathbf{WG}\\ \mathbf{0}\\ \mathbf{S}\bm{\mu} \end{bmatrix} \nonumber \\
        &= (\mathbf{WG}^T\mathbf{Wd} + \mathbf{S}^T\mathbf{S}\bm{\mu} \nonumber \\
        &= \mathbf{G}^T \mathbf{W}^2 \mathbf{d} + \mathbf{S}^2 \bm{\mu} \nonumber \\
        \nonumber \\
        \mathbf{m} &= (\mathbf{G}^T\mathbf{W}^2\mathbf{G} + \alpha^2\mathbf{L}^T\mathbf{L} + \mathbf{S}^2)^{-1} (\mathbf{G}^T\mathbf{W}^2\mathbf{d} + \mathbf{S}^2\bm{\mu}) \nonumber
    \end{align}

Expanding this to the three dimensional case and using the combined **L** defined in the previous section, we get the full solution.

.. math::
    \begin{equation}
        \mathbf{m} = (\mathbf{G}^T\mathbf{W}^2\mathbf{G} + \mathbf{L} + \mathbf{S}^2)^{-1} (\mathbf{G}^T\mathbf{W}^2\mathbf{d} + \mathbf{S}^2\bm{\mu})
    \end{equation}

:math:`\mathbf{S}` is an M x M diagonal matrix where the column and/or row number of the elements along the diagonal correspond to the index number of that model parameter, which has an associated x,y,z coordinate. Likewise, :math:`\mathbf{\mu}` is an M x 1 vector for which each row corresponds to the model parameter of that number. To apply different priors to different model parameters, one simply needs to locate which model parameters are in the region of interest, say inside the sphere, using the coordinates, then use the indices of those model parameters to assign the appropriate standard deviations and means to the elements of :math:`\mathbf{S}` and :math:`\mathbf{\mu}`, respectively, that have the same index. If an element on the diagonal of **S** equal zero, then no prior is applied to that model parameter. 

