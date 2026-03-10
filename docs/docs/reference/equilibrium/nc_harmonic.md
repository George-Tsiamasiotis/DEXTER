Single perturbation harmonic from a netCDF file of the form:

$$
\alpha_{\{m, n\}}(\psi/\psi_p) \cos\big(m\theta - n\zeta + \phi(\psi/\psi_p) \big)
$$

where $\alpha$ and $\phi$ can be expressed as functions of either or both $\psi, \psi_p$, and
are calculated by interpolation over the numerical data.

::: dexter.NcHarmonic
    options:
        inherited_members: true
        docstring_section_style: list
        summary:
            functions: true
