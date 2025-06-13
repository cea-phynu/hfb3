@mainpage

@page Formalism Formalism

# Formalism

## The HFB method
The Hartree-Fock-Bogoliubov method used in this solver relies on the following three main assumptions/principles:

1. The nucleon-nucleon effective interaction is a two-body density-dependent interaction:

\f{equation}{
\hat{H} \equiv \sum_{ab} t_{ab}c_a^\dagger c_b+\frac{1}{4}\sum_{abcd}\tilde{\nu}_{abcd}c_a^\dagger c_b^\dagger c_d c_c.
\f}

2. The stationary solutions are given by the variational principle:

\f{equation}{
\partial E[\psi]=0
\f}
with
\f{equation}{
E[\psi]\equiv\frac{\langle\psi|\hat{H}|\psi\rangle}{\langle\psi|\psi\rangle}.
\f}

3. The wave functions are defined using a Bogoliubov transformation:

\f{equation}{
|\psi\rangle \equiv \prod_b \beta_b|0\rangle
\f}
with
\f{equation}{
\beta_b^\dagger\equiv\sum_a\left(u_{ab}c_a^\dagger+v_{ab}c_a\right)
\f}
and
\f{equation}{
\beta_b|\psi\rangle = 0.
\f}

## Basis functions {#basisFunctions}

The cylindrical basis states are defined as

\f{equation}{
|\varphi_{m,n,n_z,d,s}^{(b_{\perp},b_z, d_0)}\rangle
         \equiv
         |\phi_{m,n,n_z,d}^{(b_{\perp},b_z,d_0)}\rangle \otimes|s \rangle
\f}

with

\f{equation}{
\langle \mathbf{r}
        \equiv
        (z,\mathbf{r_{\perp}})|\phi_{m,n,n_z,d}^{(b_{\perp},b_z,d_0)}\rangle
        \equiv
        \phi_{n_z}(z-z_d) \phi_{m,n}(\mathbf{r_{\perp}}),
\f}
\f{equation}{
z_d^{(d_0)}\equiv(\frac{1}{2}-d).d_0
\f}
and (Berger's conventions)
\f{equation}{
\phi_{n_z}(z)
         \equiv
         \frac{1}{\sqrt{b_z}}
         \frac{1}{\sqrt{2^{n_z} \sqrt{\pi}n_z!}}
         e^{-\frac{1}{2}z^2/b_z^2}H_{n_z}(z/b_z)
\f}
\f{equation}{
\phi_{m,n}(\mathbf{r_{\perp}}
         \equiv
         (r_{\perp},\theta))\equiv\frac{1}{b_{\perp}\sqrt{\pi}}
         \sqrt{\frac{n!}{(n+|m|)!}}
         e^{-\frac{1}{2}r_{\perp}^2/b_{\perp}^2}
         (r_{\perp}/b_{\perp})^{|m|}
         L_n^{|m|}(r_{\perp}^2/b_{\perp}^2)
         e^{im\theta}.
\f}
The quantities
\f{equation}{
Z(z, n_z, d)
\equiv
\phi_{n_z}(z - z_d)
= \frac{1}{\sqrt{b_z}}
         \frac{1}{\sqrt{2^{n_z} \sqrt{\pi}n_z!}}
         e^{-\frac{1}{2}(z - z_d)^2/b_z^2}H_{n_z}((z-z_d)/b_z)
\f}
\f{equation}{
R(r_\perp, m, n)
         \equiv
         \frac{1}{b_{\perp}\sqrt{\pi}}
         \sqrt{\frac{n!}{(n+|m|)!}}
         e^{-\frac{1}{2}r_{\perp}^2/b_{\perp}^2}
         (r_{\perp}/b_{\perp})^{|m|}
         L_n^{|m|}(r_{\perp}^2/b_{\perp}^2)
\f}
are calculated in the Basis::zPart() and Basis::rPart() methods.
One then has
\f{equation}{
  \langle \mathbf{r}
        \equiv
        (z,\mathbf{r_{\perp}})|\phi_{m,n,n_z,d}^{(b_{\perp},b_z,d_0)}\rangle
        \equiv
    Z(z, n_z, d)
    .
    R(r_\perp, m, n)
    .
         e^{im\theta} .
\f}

The \f$b_z\f$ and \f$b_{\perp}\f$ quantities are linked to the oscillator frequencies as
\f{equation}{
b_{\perp}\equiv\sqrt{\frac{\hbar}{M\omega_{\perp}}},
b_z\equiv\sqrt{\frac{\hbar}{M\omega_z}}
\f}
where \f$M\f$ is the nucleus mass.

## Links with HFB2CT solver

The `alf` and `bet` variables in `HFB2CT` are defined as
\f{equation}{
\begin{alignedat}{30}
\alpha&\equiv& \frac{1}{b_z^2}\\
\beta&\equiv& \frac{1}{b_{\perp}^2}.
\end{alignedat}
\f}
The `q1` variable in HFB2CT is defined as
\f{equation}{
\textrm{q1}\equiv \frac{\hbar \omega_{z}}{\hbar \omega_{\perp}}
         =
         \frac{\alpha}{\beta}
         =
         \frac{b_{\perp}^2}{b_z^2}.
\f}

## Quadratures

Gauss-Hermite :
\f{equation}{
\int_{-\infty}^\infty dx P(x) e^{-x^2} = \sum_i w_i^{\mathrm{He}} P(x^{\mathrm{He}}_i)
\f}
Gauss-Laguerre :
\f{equation}{
\int_0^\infty dx P(x) e^{-x} = \sum_i w_i^{\mathrm{La}} P(x^{\mathrm{La}}_i)
\f}
Gauss-Legendre :
\f{equation}{
\int_{-1}^1 dx P(x) = \sum_i w_i^{\mathrm{Le}} P(x^{\mathrm{Le}}_i)
\f}

We introduce the normalized basis functions (calculated in Basis::rPartNorm() and Basis::zPartNorm() ) as

\f{equation}{ R_{m,n}^{(N)}(\eta) \equiv \mathcal{N}_{m,n} L_{|m|,n}(\eta) \eta^{|m|/2} e^{-\eta/2} \f}
\f{equation}{ Z_{n_z}^{(N)}(\zeta) \equiv \mathcal{N}_{n_z} H_{n_z}(\zeta) e^{-\zeta^2/2} \f}
with
\f{equation}{ \mathcal{N}_{m,n} \equiv \sqrt{\frac{n_\perp!}{(n_\perp+|m|)!}} \f}
\f{equation}{ \mathcal{N}_{n_z} \equiv \frac{1}{\sqrt{2^{n_z}\sqrt{\pi}n_z!}}. \f}
We have the important relations
\f{equation}{
R(r_\perp, m, n)
  = \frac{1}{b_\perp\sqrt{\pi}}
  R^{(N)}_{m, n}(r_\perp^2/b_\perp^2)
\f}
\f{equation}{
Z(z, n_z, d)
  = \frac{1}{\sqrt{b_z}}
  Z^{(N)}_{n_z}((z - z_d)/b_z),
\f}
and we can write
\f{equation}{
  \langle \mathbf{r} \equiv (z,\mathbf{r_{\perp}})|\phi_{m,n,n_z,d}^{(b_{\perp},b_z,d_0)}\rangle
  = \frac{1}{b_\perp\sqrt{b_z}}
  Z_{n_z}^{(N)}((z-z_d)/b_z)
  R_{m, n}^{(N)}(r_\perp^2/b_\perp^2)e^{im\theta}.
\f}
To use the different quadratures, we also introduce the normalized reduced quantities
\f{equation}{ R_{m,n}^{(NR)}(\eta) \equiv \mathcal{N}_{m,n} L_{m,n}(\eta) \eta^{m/2} \f}
\f{equation}{ Z_{n_z}^{(NR)}(\zeta) \equiv \mathcal{N}_{n_z} H_{n_z}(\zeta). \f}

# Representations

## Orthonormal representation

Let us define
\f{equation}{
\begin{alignedat}{30}
  t_{n_{z1},n_{z2}}^{(b_1,b_2,d_1,d_2)}\equiv& \int \phi_{n_{z1}}^{(b_1)}(z-d_1)\phi_{n_{z2}}^{(b_2)}(z-d_2)dz\\
  T_{a, b}^{+-}\equiv&  t_{a, b}^{(+\frac{d_0}{2}, -\frac{d_0}{2}, b_z, b_z)}\\
  T_{a, b}^{-+}\equiv&  t_{a, b}^{(-\frac{d_0}{2}, +\frac{d_0}{2}, b_z, b_z)} = T_{b,a}^{+-}.
\end{alignedat}
\f}

It is calculated using the relations

\f{equation}{
\begin{alignedat}{30}
   \alpha_1\equiv& \frac{1}{b_1^2}\\
   \alpha_2\equiv& \frac{1}{b_2^2}\\
   (\alpha_1+\alpha_2)\sqrt{n_{z1}}t_{n_{z1},n_{z2}}^{(b_1,b_2,d_1,d_2)}=&
   (d_2-d_1)\alpha_2\sqrt{2\alpha_1}t_{n_{z1}-1,n_{z2}}^{(b_1,b_2,d_1,d_2)}+\\
   &2\sqrt{\alpha_1\alpha_2n_{z2}}t_{n_{z1}-1,n_{z2}}^{(b_1,b_2,d_1,d_2)}+
   (\alpha1-\alpha2)\sqrt{n_{z1}-1}t_{n_{z1}-2,n_{z2}}^{(b_1,b_2,d_1,d_2)}\\
   (\alpha_1+\alpha_2)\sqrt{n_z}t_{0,n_{z}}^{(b_1,b_2,d_1,d_2)}=&
   (d_1-d_2)\alpha_1\sqrt{2\alpha_2}t_{0,n_{z}-1}^{(b_1,b_2,d_1,d_2)}+
   (\alpha2-\alpha1)\sqrt{n_z-1}t_{0,n_{z}-2}^{(b_1,b_2,d_1,d_2)}\\
   t_{0,0}^{(b_1,b_2,d_1,d_2)}=&
   \sqrt{\frac{2\sqrt{\alpha_1\alpha_2}}{\alpha_1+\alpha_2}}
   e^{-\frac{(d_1-d_2)^2}{2}\frac{\alpha_1\alpha_2}{\alpha_1+\alpha_2}},
\end{alignedat}
\f}

c.f. subroutine trz12() in HFB2CT.

The basis overlap is defined as

\f{equation}{
\begin{alignedat}{30}
  T\equiv&[\langle a|b\rangle]^{a,b}\\
  =&I^{m_a,m_b}\times I^{n_a,n_b}\times \left(\begin{array}{cc}I^{n_{za},n_{zb}}&T_{10}^{n_{za},n_{zb}}\\T_{01}^{n_{za},n_{zb}}&I^{n_{za},n_{zb}}\end{array}\right)^{d_a,d_b}\times I^{s_a,s_b}.
\end{alignedat}
\f}
with
\f{equation}{
\begin{alignedat}{30}
[M_{a,b,c,d}]^{a,c}\equiv& \left(\begin{array}{ccc}M_{0,b,0,d}&M_{0,b,1,d}&\ldots\\M_{1,b,0,d}&M_{1,b,1,d}&\ldots\\\vdots&\vdots&\ddots\\\end{array}\right)\\
I^{a,b}\equiv& [\delta_{a,b}]^{a,b}\\
\left(\begin{array}{ccc}A&B&\ldots\\C&D&\ldots\\\vdots&\vdots&\ddots\end{array}\right)\times M\equiv& \left(\begin{array}{ccc}AM&BM&\ldots\\CM&DM&\ldots\\\vdots&\vdots&\ddots\end{array}\right)\\
T_{01}^{n_{za},n_{zb}}\equiv&[T^{-+}_{n_{za},n_{zb}}]^{n_{za},n_{zb}}\\
T_{10}^{n_{za},n_{zb}}\equiv&[T^{+-}_{n_{za},n_{zb}}]^{n_{za},n_{zb}}.
\end{alignedat}
\f}





@section basisTruncation Basis truncation

Truncating a basis means selecting which basis states belong to a Basis object. The basis truncation used follows the following algorithm:

We define

\f{equation}{
\nu(i) \equiv (N+2) Q^{2/3} + \frac{1}{2} - i.Q
\f}
and
\f{equation}{
m^\textrm{max} \equiv \textrm{sup} \{ i:\nu(i) \ge 1 \}
\f}
The set of included basis states is then defined by
\f{equation}{
\left\{
\begin{array}{rcl}
0 &\le& m \lt \textrm{mMax} \equiv m^\textrm{max}
\end{array}
\right.
\f}
















@class SolverHFBBroyden

TODO TODO TODO

Example of some inline math : \f$a=b^5\f$.

\f[
a=b^2\tag{Plop}
\f]

\f{eqnarray*}{
    g &=& \frac{Gm_2}{r^2} \\
      &=& \frac{(6.673 \times 10^{-11}\,\mbox{m}^3\,\mbox{kg}^{-1}\,
          \mbox{s}^{-2})(5.9736 \times 10^{24}\,\mbox{kg})}{(6371.01\,\mbox{km})^2} \\
      &=& 9.82066032\,\mbox{m/s}^2
\f}




// @mainpage
//
// @section Formalism
//
// Basis formalism
//
// Example of some inline math : \f$a=b^5\f$.
//
// \f[
// a=b^2\tag{Plop}
// \f]
//
// \f{eqnarray*}{
//     g &=& \frac{Gm_2}{r^2} \\
//       &=& \frac{(6.673 \times 10^{-11}\,\mbox{m}^3\,\mbox{kg}^{-1}\,
//           \mbox{s}^{-2})(5.9736 \times 10^{24}\,\mbox{kg})}{(6371.01\,\mbox{km})^2} \\
//       &=& 9.82066032\,\mbox{m/s}^2
// \f}
//
//
//
//
