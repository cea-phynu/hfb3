#!/usr/bin/env python3

# ==============================================================================
# ==============================================================================
# ==============================================================================

import click
import sys
from sympy import symbols, Lambda, N, Float
from sympy.printing.str import Printer
from sympy import sqrt, exp, pi, I, factorial, Integral, Abs, diff
from sympy import assoc_laguerre, hermite, oo, summation, Min, Max, Ynm
from sympy import cos, sin

# ==============================================================================
# ==============================================================================
# ==============================================================================

Printer.set_global_settings(min=-2, max=4)
keys = {}

# ==============================================================================
# ==============================================================================
# ==============================================================================


def WVAL(label, args, val, fp, prec, bar):
    global keys
    if label not in keys:
        keys[label] = 0

    for ia, arg in enumerate(args):
        fp.write(f"#define {label}_{keys[label]:02d}_ARG{ia:01d}  {str(arg)[1:-1]}\n")

    fp.write(f"#define {label}_{keys[label]:02d}_TARGET {N(val, prec)}\n")

    keys[label] += 1

    if bar:
        bar.update(1)

# ==============================================================================
# ==============================================================================
# ==============================================================================


def WBAS(label, ref, fp):
    if (ref[0] > 0.0):
        fp.write(f'#define {label}_STATE "examples/42Ca_deformed_2x9.msg.gz"\n')
    else:
        fp.write(f'#define {label}_STATE "examples/42Ca_deformed_1x11.msg.gz"\n')

# ==============================================================================
# ==============================================================================
# ==============================================================================


def WACT(label, fp):
    fp.write(f'#define {label}\n')

# ==============================================================================
# ==============================================================================
# ==============================================================================


def main():

    # ==========================================================================

    basis_zpart_1ct                              = True
    basis_zpart_2ct                              = True
    basis_rpart_1ct                              = True
    basis_rpart_2ct                              = True
    basis_z0recouv_1ct                           = True
    basis_z0recouv_2ct                           = True
    basis_irecouv_1ct                            = True
    basis_irecouv_2ct                            = True
    basis_izrecouv_1ct                           = True
    basis_izrecouv_2ct                           = True
    basis_hcyl_1ct                               = True
    basis_hcyl_2ct                               = True
    basis_overlapz_12                            = True
    basis_overlapz_22                            = True
    basis_overlapr_12                            = True
    basis_overlapr_2X                            = True
    basis_talmanz_1ct                            = True
    basis_talmanz_2ct                            = False
    fieldcentraldirect_calcgaussianintegralz_1ct = True
    fieldcentraldirect_calcgaussianintegralz_2ct = True
    multipole_operators_znrecouv_1ct             = True
    multipole_operators_znrecouv_2ct             = True
    multipole_operators_rnrecouv                 = True
    multipole_operators_rnodd                    = True
    multipole_operators_qlmho_1ct                = True
    multipole_operators_qlmho_2ct                = True

    # ==========================================================================

    # --- test_1ct.rhl.dat.gz ---
    # vd0 = 0
    # vbr = Float('2.4144410198519912', prec)
    # vbz = Float('3.1170299534289376', prec)
    # basis1ct = (vd0, vbr, vbz)
    basis1ct = (0, 2.4144410198519912, 3.1170299534289376)

    # ==========================================================================

    # --- test_2ct.rhl.dat.gz ---
    # vd0 = Float('5.0', prec)
    # vbr = Float('2.3878704919130063', prec)
    # vbz = Float('3.1087425890539820', prec)
    # basis2ct = (vd0, vbr, vbz)
    basis2ct = (5, 2.3878704919130063, 3.1087425890539820)

    # ==========================================================================

    m, n, nz, d, s      = symbols("m n n_z d s")
    ma, na, nza, da, sa = symbols("m_a n_a n^z_a d_a s_a")
    mb, nb, nzb, db, sb = symbols("m_b n_b n^z_b d_b s_b")
    mc, nc, nzc, dc, sc = symbols("m_c n_c n^z_c d_c s_c")
    md, nd, nzd, dd, sd = symbols("m_d n_d n^z_d d_d s_d")
    br, bz, d0, p       = symbols("b_\\perp b_z d_0 p")
    rp, z, zpa, zpb     = symbols("r_\\perp z Z_a Z_b")
    rpa, rpb            = symbols("R_a R_b")
    theta, phi, x, y    = symbols("theta phi x y")
    u, r, k, z          = symbols("u r k z")

    # ==========================================================================
    # ==========================================================================
    # ==========================================================================

    with click.progressbar(length=1, label=f"{'Def   - zpart':20}") as bar:

        zpart = Lambda((z, nz, d, d0, br, bz),
                       exp(-(z - (1 / 2 - d) * d0)**2 / (2 * bz**2)) *
                       hermite(nz, (z - (1 / 2 - d) * d0) / bz) /
                       sqrt(bz * sqrt(pi) * 2**nz * factorial(nz)))

        bar.update(1)

    # ==========================================================================
    # ==========================================================================
    # ==========================================================================

    with click.progressbar(length=1, label=f"{'Def   - rpart':20}") as bar:

        rpart = Lambda((rp, m, n, d0, br, bz),
                       exp(-rp**2 / (2 * br**2)) *
                       assoc_laguerre(n, m, (rp**2 / br**2)) *
                       (rp / br)**m /
                       (br * sqrt(pi)) *
                       sqrt(factorial(n) / factorial(n + m)))

        bar.update(1)

    # ==========================================================================
    # ==========================================================================
    # ==========================================================================

    with click.progressbar(length=1, label=f"{'Def   - z0recouv':20}") as bar:

        z0recouv = Lambda((nza, nzb, da, db, d0, br, bz),
                         Integral(zpart(z, nza, da, d0, br, bz) *
                                  zpart(z, nzb, db, d0, br, bz),
                                  (z, -oo, oo)))

        bar.update(1)

    # ==========================================================================
    # ==========================================================================
    # ==========================================================================

    with click.progressbar(length=1, label=f"{'Def   - irecouv':20}") as bar:

        irecouv = Lambda((nza, nzb, da, db, d0, br, bz),
                         Integral(zpart(z, nza, da, d0, br, bz) *
                                  zpart(z, nzb, db, d0, br, bz),
                                  (z, 0, oo)))

        bar.update(1)

    # ==========================================================================
    # ==========================================================================
    # ==========================================================================

    with click.progressbar(length=1, label=f"{'Def   - izrecouv':20}") as bar:

        izrecouv = Lambda((nza, nzb, da, db, d0, br, bz),
                         Integral(zpart(z, nza, da, d0, br, bz) *
                                  zpart(z, nzb, db, d0, br, bz) * z,
                                  (z, 0, oo)))

        bar.update(1)

    # ==========================================================================
    # ==========================================================================
    # ==========================================================================

    with click.progressbar(length=1, label=f"{'Def   - hcyl':20}") as bar:

        hcyl = Lambda((m, n, nza, nzb, da, db, d0, br, bz),
                      20.735 *
                      ( -Integral(          zpart(z, nza, da, d0, br, bz) *
                                  diff(zpart(z, nzb, db, d0, br, bz), z, z, evaluate=False),
                                 (z, -oo, oo)) +
                         Integral(zpart(z, nza, da, d0, br, bz) *
                                  zpart(z, nzb, db, d0, br, bz) *
                                  (Abs(z) - d0 / 2)**2 / (bz**4),
                                 (z, -oo, oo)) +
                        2 * (2 * n + m + 1) / br**2 * z0recouv(nza, nzb, da, db, d0, br, bz)
                      )
                     )

        bar.update(1)

    # ==========================================================================
    # ==========================================================================
    # ==========================================================================

    with click.progressbar(length=1, label=f"{'Def   - overlapz':20}") as bar:

        overlapz = Lambda((zpa, zpb), Integral(zpa * zpb, (z, -oo, oo)))

        bar.update(1)

    # ==========================================================================
    # ==========================================================================
    # ==========================================================================

    with click.progressbar(length=1, label=f"{'Def   - overlapr':20}") as bar:

        overlapr = Lambda((rpa, rpb), Integral(rpa * rpb * rp * 2 * pi, (rp, 0, oo)))

        bar.update(1)

    # ==========================================================================
    # ==========================================================================
    # ==========================================================================

    with click.progressbar(length=2, label=f"{'Def   - talmanz':20}") as bar:

        mnz = Lambda((nz), sqrt(2**nz / factorial(nz)))

        bar.update(1)

        talmanz = Lambda((nza, nzb, da, db, nzc, d0, br, bz),
                         1 / (mnz(nza) * mnz(nzb) * mnz(nzc)) *
                         summation(
                                   (mnz(nzd) * mnz(nzc - nzd))**2 *
                                   mnz(nza - nzd) *
                                   mnz(nzb - nzc + nzd) *
                                   z0recouv(nza - nzd, nzb - nzc + nzd, da, db, d0, br, bz),
                                   (nzd, Max(nzc - nzb, 0) , Min(nzc, nza))
                                  )
                        )

        bar.update(1)

    # ==========================================================================
    # ==========================================================================
    # ==========================================================================

    with click.progressbar(length=2, label=f"{'Def   - gIntz':20}") as bar:

        apart = Lambda((z, nz, bz),
                       exp(-z**2 / (2 * bz**2)) *
                       hermite(nz, z / bz) /
                       sqrt(bz * sqrt(pi) * 2**nz * factorial(nz)))

        bar.update(1)

        gaussianIntegralz = Lambda((da, db, dc, dd, nz, p, d0, br, bz),
                         Integral(apart(z - (
                                             (1 / 2 - da) +
                                             (1 / 2 - dc) -
                                             (1 / 2 - db) -
                                             (1 / 2 - dd)
                                            ) * d0 / (2 * sqrt(2)), 0, bz) *
                                  apart(z - (
                                             (1 / 2 - da) +
                                             (1 / 2 - dc) -
                                             (1 / 2 - db) -
                                             (1 / 2 - dd)
                                            ) * d0 / (2 * sqrt(2)), nz, bz) *
                                  exp(-2 * z**2 / p**2),
                                  (z, -oo, oo)))

        bar.update(1)

    # ==========================================================================
    # ==========================================================================
    # ==========================================================================

    with click.progressbar(length=1, label=f"{'Def   - znrecouv':20}") as bar:

        znrecouv = Lambda((n, nza, nzb, da, db, d0, br, bz),
                         Integral(zpart(z, nza, da, d0, br, bz) *
                                  zpart(z, nzb, db, d0, br, bz) * z**n,
                                  (z, -oo, oo)))

        bar.update(1)

    # ==========================================================================
    # ==========================================================================
    # ==========================================================================

    with click.progressbar(length=1, label=f"{'Def   - rnrecouv':20}") as bar:

        rnrecouv = Lambda((n, ma, mb, na, nb, d0, br, bz),
                         Integral(rpart(rp, ma, na, d0, br, bz) *
                                  rpart(rp, mb, nb, d0, br, bz) * rp**n * rp * 2 * pi,
                                  (rp, 0, oo)))

        bar.update(1)

    # ==========================================================================
    # ==========================================================================
    # ==========================================================================

    with click.progressbar(length=1, label=f"{'Def   - rnodd':20}") as bar:

        rnodd = Lambda((k, ma, na, mb, nb, d0, br, bz),
                         Integral(rpart(rp, ma, na, d0, br, bz) *
                                  rpart(rp, mb, nb, d0, br, bz) * rp ** k,
                                  (rp, 0, oo)))

        bar.update(1)

    # ==========================================================================
    # ==========================================================================
    # ==========================================================================

    valList = [(lbd, 0) for lbd in range(7)] + [(2, -2), (2, -1), (2, 1), (2, 2)]

    with click.progressbar(length=len(valList), label=f"{'Def   - qlmHO':20}") as bar:

        Ylm, Qlm, qlmHO = {}, {}, {}
        for lbd, mu in valList:
            Ylm[(lbd, mu)] = Ynm(lbd, mu, theta, phi).expand(func=True).factor()

        for lbd, mu in valList:

            q = (r**lbd * Ylm[(lbd, mu)]).expand(func=True)
            q = q.subs(exp( I * phi) * sin(theta), u).factor().subs(u, (x + I * y) / r)
            q = q.subs(exp(-I * phi) * sin(theta), u).factor().subs(u, (x - I * y) / r)
            q = q.subs(cos(theta), z / r)
            q = q.subs(x + I * y, exp( I * phi) * rp)
            q = q.subs(x - I * y, exp(-I * phi) * rp)
            q = q.subs(r, sqrt(rp**2 + z**2))
            Qlm[(lbd, mu)] = q.expand().factor()

            qlmHO[(lbd, mu)] = Lambda((ma, na, nza, da, sa, mb, nb, nzb, db, sb, d0, br, bz),
                                      Integral(
                                        Integral(
                                          Integral(zpart(z, nza, da, d0, br, bz) *
                                                zpart(z, nzb, db, d0, br, bz) *
                                                rpart(rp, ma, na, d0, br, bz) *
                                                rpart(rp, mb, nb, d0, br, bz) *
                                                exp(-I * ma * phi) *
                                                exp( I * mb * phi) *
                                                rp *
                                                Qlm[(lbd, mu)],
                                                (z, -oo, oo)),
                                          (rp, 0, oo)),
                                        (phi, 0, 2 * pi))
                                      )

            bar.update(1)

    # ==========================================================================

    fp = open("../../misc/tests/maxima_targets2.h", "w+")
    prec = 20

    # ==========================================================================

    if basis_zpart_1ct:
        with click.progressbar(length=21, label=f"{'Calc. - zpart 1ct':20}") as bar:

            WACT("BASIS_ZPART_1CT", fp)

            c = basis1ct
            WBAS("BASIS_ZPART_1CT", c, fp)
            vz  = Float('-2.238', prec)

            a = (vz,  0, 0); WVAL("BASIS_ZPART_1CT", (a,), zpart(*a, *c), fp, prec, bar)
            a = (vz,  1, 1); WVAL("BASIS_ZPART_1CT", (a,), zpart(*a, *c), fp, prec, bar)
            a = (vz,  2, 0); WVAL("BASIS_ZPART_1CT", (a,), zpart(*a, *c), fp, prec, bar)
            a = (vz,  3, 1); WVAL("BASIS_ZPART_1CT", (a,), zpart(*a, *c), fp, prec, bar)
            a = (vz,  4, 0); WVAL("BASIS_ZPART_1CT", (a,), zpart(*a, *c), fp, prec, bar)
            a = (vz,  5, 1); WVAL("BASIS_ZPART_1CT", (a,), zpart(*a, *c), fp, prec, bar)
            a = (vz,  6, 0); WVAL("BASIS_ZPART_1CT", (a,), zpart(*a, *c), fp, prec, bar)
            a = (vz,  7, 1); WVAL("BASIS_ZPART_1CT", (a,), zpart(*a, *c), fp, prec, bar)
            a = (vz,  8, 0); WVAL("BASIS_ZPART_1CT", (a,), zpart(*a, *c), fp, prec, bar)
            a = (vz,  9, 1); WVAL("BASIS_ZPART_1CT", (a,), zpart(*a, *c), fp, prec, bar)
            a = (vz, 10, 0); WVAL("BASIS_ZPART_1CT", (a,), zpart(*a, *c), fp, prec, bar)
            a = (vz, 11, 1); WVAL("BASIS_ZPART_1CT", (a,), zpart(*a, *c), fp, prec, bar)
            a = (vz, 12, 0); WVAL("BASIS_ZPART_1CT", (a,), zpart(*a, *c), fp, prec, bar)
            a = (vz, 13, 1); WVAL("BASIS_ZPART_1CT", (a,), zpart(*a, *c), fp, prec, bar)
            a = (vz, 14, 0); WVAL("BASIS_ZPART_1CT", (a,), zpart(*a, *c), fp, prec, bar)
            a = (vz, 15, 1); WVAL("BASIS_ZPART_1CT", (a,), zpart(*a, *c), fp, prec, bar)
            a = (vz, 16, 0); WVAL("BASIS_ZPART_1CT", (a,), zpart(*a, *c), fp, prec, bar)
            a = (vz, 17, 1); WVAL("BASIS_ZPART_1CT", (a,), zpart(*a, *c), fp, prec, bar)
            a = (vz, 18, 0); WVAL("BASIS_ZPART_1CT", (a,), zpart(*a, *c), fp, prec, bar)
            a = (vz, 19, 1); WVAL("BASIS_ZPART_1CT", (a,), zpart(*a, *c), fp, prec, bar)
            a = (vz, 20, 0); WVAL("BASIS_ZPART_1CT", (a,), zpart(*a, *c), fp, prec, bar)

    # ==========================================================================

    if basis_zpart_2ct:
        with click.progressbar(length=21, label=f"{'Calc. - zpart 2ct':20}") as bar:

            WACT("BASIS_ZPART_2CT", fp)

            c = basis2ct
            WBAS("BASIS_ZPART_2CT", c, fp)
            vz  = Float('-2.238', prec)

            a = (vz,  0, 0); WVAL("BASIS_ZPART_2CT", (a,), zpart(*a, *c), fp, prec, bar)
            a = (vz,  1, 1); WVAL("BASIS_ZPART_2CT", (a,), zpart(*a, *c), fp, prec, bar)
            a = (vz,  2, 0); WVAL("BASIS_ZPART_2CT", (a,), zpart(*a, *c), fp, prec, bar)
            a = (vz,  3, 1); WVAL("BASIS_ZPART_2CT", (a,), zpart(*a, *c), fp, prec, bar)
            a = (vz,  4, 0); WVAL("BASIS_ZPART_2CT", (a,), zpart(*a, *c), fp, prec, bar)
            a = (vz,  5, 1); WVAL("BASIS_ZPART_2CT", (a,), zpart(*a, *c), fp, prec, bar)
            a = (vz,  6, 0); WVAL("BASIS_ZPART_2CT", (a,), zpart(*a, *c), fp, prec, bar)
            a = (vz,  7, 1); WVAL("BASIS_ZPART_2CT", (a,), zpart(*a, *c), fp, prec, bar)
            a = (vz,  8, 0); WVAL("BASIS_ZPART_2CT", (a,), zpart(*a, *c), fp, prec, bar)
            a = (vz,  9, 1); WVAL("BASIS_ZPART_2CT", (a,), zpart(*a, *c), fp, prec, bar)
            a = (vz, 10, 0); WVAL("BASIS_ZPART_2CT", (a,), zpart(*a, *c), fp, prec, bar)
            a = (vz, 11, 1); WVAL("BASIS_ZPART_2CT", (a,), zpart(*a, *c), fp, prec, bar)
            a = (vz, 12, 0); WVAL("BASIS_ZPART_2CT", (a,), zpart(*a, *c), fp, prec, bar)
            a = (vz, 13, 1); WVAL("BASIS_ZPART_2CT", (a,), zpart(*a, *c), fp, prec, bar)
            a = (vz, 14, 0); WVAL("BASIS_ZPART_2CT", (a,), zpart(*a, *c), fp, prec, bar)
            a = (vz, 15, 1); WVAL("BASIS_ZPART_2CT", (a,), zpart(*a, *c), fp, prec, bar)
            a = (vz, 16, 0); WVAL("BASIS_ZPART_2CT", (a,), zpart(*a, *c), fp, prec, bar)
            a = (vz, 17, 1); WVAL("BASIS_ZPART_2CT", (a,), zpart(*a, *c), fp, prec, bar)
            a = (vz, 18, 0); WVAL("BASIS_ZPART_2CT", (a,), zpart(*a, *c), fp, prec, bar)
            a = (vz, 19, 1); WVAL("BASIS_ZPART_2CT", (a,), zpart(*a, *c), fp, prec, bar)
            a = (vz, 20, 0); WVAL("BASIS_ZPART_2CT", (a,), zpart(*a, *c), fp, prec, bar)

    # ==========================================================================

    if basis_rpart_1ct:
        with click.progressbar(length=21, label=f"{'Calc. - rpart 1ct':20}") as bar:

            WACT("BASIS_RPART_1CT", fp)

            c = basis1ct
            WBAS("BASIS_RPART_1CT", c, fp)
            vr  = Float('2.238', prec)

            a = (vr,  0,  2); WVAL("BASIS_RPART_1CT", (a,), rpart(*a, *c), fp, prec, bar)
            a = (vr,  1,  1); WVAL("BASIS_RPART_1CT", (a,), rpart(*a, *c), fp, prec, bar)
            a = (vr,  2,  3); WVAL("BASIS_RPART_1CT", (a,), rpart(*a, *c), fp, prec, bar)
            a = (vr,  3,  2); WVAL("BASIS_RPART_1CT", (a,), rpart(*a, *c), fp, prec, bar)
            a = (vr,  4,  1); WVAL("BASIS_RPART_1CT", (a,), rpart(*a, *c), fp, prec, bar)
            a = (vr,  5,  4); WVAL("BASIS_RPART_1CT", (a,), rpart(*a, *c), fp, prec, bar)
            a = (vr,  6,  3); WVAL("BASIS_RPART_1CT", (a,), rpart(*a, *c), fp, prec, bar)
            a = (vr,  7,  4); WVAL("BASIS_RPART_1CT", (a,), rpart(*a, *c), fp, prec, bar)
            a = (vr,  8,  2); WVAL("BASIS_RPART_1CT", (a,), rpart(*a, *c), fp, prec, bar)
            a = (vr,  9,  5); WVAL("BASIS_RPART_1CT", (a,), rpart(*a, *c), fp, prec, bar)
            a = (vr, 10,  2); WVAL("BASIS_RPART_1CT", (a,), rpart(*a, *c), fp, prec, bar)
            a = (vr, 11, 11); WVAL("BASIS_RPART_1CT", (a,), rpart(*a, *c), fp, prec, bar)
            a = (vr, 12,  3); WVAL("BASIS_RPART_1CT", (a,), rpart(*a, *c), fp, prec, bar)
            a = (vr, 13, 12); WVAL("BASIS_RPART_1CT", (a,), rpart(*a, *c), fp, prec, bar)
            a = (vr, 14,  1); WVAL("BASIS_RPART_1CT", (a,), rpart(*a, *c), fp, prec, bar)
            a = (vr, 15, 14); WVAL("BASIS_RPART_1CT", (a,), rpart(*a, *c), fp, prec, bar)
            a = (vr, 16,  3); WVAL("BASIS_RPART_1CT", (a,), rpart(*a, *c), fp, prec, bar)
            a = (vr, 17, 14); WVAL("BASIS_RPART_1CT", (a,), rpart(*a, *c), fp, prec, bar)
            a = (vr, 18,  2); WVAL("BASIS_RPART_1CT", (a,), rpart(*a, *c), fp, prec, bar)
            a = (vr, 19, 15); WVAL("BASIS_RPART_1CT", (a,), rpart(*a, *c), fp, prec, bar)
            a = (vr, 20,  5); WVAL("BASIS_RPART_1CT", (a,), rpart(*a, *c), fp, prec, bar)

    # ==========================================================================

    if basis_rpart_2ct:
        with click.progressbar(length=21, label=f"{'Calc. - rpart 2ct':20}") as bar:

            WACT("BASIS_RPART_2CT", fp)

            c = basis2ct
            WBAS("BASIS_RPART_2CT", c, fp)
            vr  = Float('-2.238', prec)

            a = (vr,  0,  2); WVAL("BASIS_RPART_2CT", (a,), rpart(*a, *c), fp, prec, bar)
            a = (vr,  1,  1); WVAL("BASIS_RPART_2CT", (a,), rpart(*a, *c), fp, prec, bar)
            a = (vr,  2,  3); WVAL("BASIS_RPART_2CT", (a,), rpart(*a, *c), fp, prec, bar)
            a = (vr,  3,  2); WVAL("BASIS_RPART_2CT", (a,), rpart(*a, *c), fp, prec, bar)
            a = (vr,  4,  1); WVAL("BASIS_RPART_2CT", (a,), rpart(*a, *c), fp, prec, bar)
            a = (vr,  5,  4); WVAL("BASIS_RPART_2CT", (a,), rpart(*a, *c), fp, prec, bar)
            a = (vr,  6,  3); WVAL("BASIS_RPART_2CT", (a,), rpart(*a, *c), fp, prec, bar)
            a = (vr,  7,  4); WVAL("BASIS_RPART_2CT", (a,), rpart(*a, *c), fp, prec, bar)
            a = (vr,  8,  2); WVAL("BASIS_RPART_2CT", (a,), rpart(*a, *c), fp, prec, bar)
            a = (vr,  9,  5); WVAL("BASIS_RPART_2CT", (a,), rpart(*a, *c), fp, prec, bar)
            a = (vr, 10,  2); WVAL("BASIS_RPART_2CT", (a,), rpart(*a, *c), fp, prec, bar)
            a = (vr, 11, 11); WVAL("BASIS_RPART_2CT", (a,), rpart(*a, *c), fp, prec, bar)
            a = (vr, 12,  3); WVAL("BASIS_RPART_2CT", (a,), rpart(*a, *c), fp, prec, bar)
            a = (vr, 13, 12); WVAL("BASIS_RPART_2CT", (a,), rpart(*a, *c), fp, prec, bar)
            a = (vr, 14,  1); WVAL("BASIS_RPART_2CT", (a,), rpart(*a, *c), fp, prec, bar)
            a = (vr, 15, 14); WVAL("BASIS_RPART_2CT", (a,), rpart(*a, *c), fp, prec, bar)
            a = (vr, 16,  3); WVAL("BASIS_RPART_2CT", (a,), rpart(*a, *c), fp, prec, bar)
            a = (vr, 17, 14); WVAL("BASIS_RPART_2CT", (a,), rpart(*a, *c), fp, prec, bar)
            a = (vr, 18,  2); WVAL("BASIS_RPART_2CT", (a,), rpart(*a, *c), fp, prec, bar)
            a = (vr, 19, 15); WVAL("BASIS_RPART_2CT", (a,), rpart(*a, *c), fp, prec, bar)
            a = (vr, 20,  5); WVAL("BASIS_RPART_2CT", (a,), rpart(*a, *c), fp, prec, bar)

    # ==========================================================================

    if basis_z0recouv_1ct:
        with click.progressbar(length=6, label=f"{'Calc. - z0recouv 1ct':20}") as bar:

            WACT("BASIS_Z0RECOUV_1CT", fp)

            c = basis1ct
            WBAS("BASIS_Z0RECOUV_1CT", c, fp)

            a = ( 5,  3, 0, 0); WVAL("BASIS_Z0RECOUV_1CT", (a,), z0recouv(*a, *c), fp, prec, bar)
            a = ( 3,  5, 0, 0); WVAL("BASIS_Z0RECOUV_1CT", (a,), z0recouv(*a, *c), fp, prec, bar)
            a = ( 2,  7, 0, 0); WVAL("BASIS_Z0RECOUV_1CT", (a,), z0recouv(*a, *c), fp, prec, bar)
            a = ( 7,  2, 0, 0); WVAL("BASIS_Z0RECOUV_1CT", (a,), z0recouv(*a, *c), fp, prec, bar)
            a = (10, 10, 0, 0); WVAL("BASIS_Z0RECOUV_1CT", (a,), z0recouv(*a, *c), fp, prec, bar)
            a = ( 9, 10, 0, 0); WVAL("BASIS_Z0RECOUV_1CT", (a,), z0recouv(*a, *c), fp, prec, bar)

    # ==========================================================================

    if basis_z0recouv_2ct:
        with click.progressbar(length=6, label=f"{'Calc. - z0recouv 2ct':20}") as bar:

            WACT("BASIS_Z0RECOUV_2CT", fp)

            c = basis2ct
            WBAS("BASIS_Z0RECOUV_2CT", c, fp)

            a = ( 5,  3, 0, 1); WVAL("BASIS_Z0RECOUV_2CT", (a,), z0recouv(*a, *c), fp, prec, bar)
            a = ( 3,  5, 0, 1); WVAL("BASIS_Z0RECOUV_2CT", (a,), z0recouv(*a, *c), fp, prec, bar)
            a = ( 2,  7, 0, 1); WVAL("BASIS_Z0RECOUV_2CT", (a,), z0recouv(*a, *c), fp, prec, bar)
            a = ( 7,  2, 0, 1); WVAL("BASIS_Z0RECOUV_2CT", (a,), z0recouv(*a, *c), fp, prec, bar)
            a = (10, 10, 0, 1); WVAL("BASIS_Z0RECOUV_2CT", (a,), z0recouv(*a, *c), fp, prec, bar)
            a = ( 9, 10, 0, 1); WVAL("BASIS_Z0RECOUV_2CT", (a,), z0recouv(*a, *c), fp, prec, bar)

    # ==========================================================================

    if basis_irecouv_1ct:
        with click.progressbar(length=4, label=f"{'Calc. - irecouv 1ct':20}") as bar:

            WACT("BASIS_IRECOUV_1CT", fp)

            c = basis1ct
            WBAS("BASIS_IRECOUV_1CT", c, fp)

            a = (9, 10, 0, 0); WVAL("BASIS_IRECOUV_1CT", (a,), irecouv(*a, *c), fp, prec, bar)
            a = (9, 10, 0, 1); WVAL("BASIS_IRECOUV_1CT", (a,), irecouv(*a, *c), fp, prec, bar)
            a = (9, 10, 1, 0); WVAL("BASIS_IRECOUV_1CT", (a,), irecouv(*a, *c), fp, prec, bar)
            a = (9, 10, 1, 1); WVAL("BASIS_IRECOUV_1CT", (a,), irecouv(*a, *c), fp, prec, bar)

    # ==========================================================================

    if basis_irecouv_2ct:
        with click.progressbar(length=4, label=f"{'Calc. - irecouv 2ct':20}") as bar:

            WACT("BASIS_IRECOUV_2CT", fp)

            c = basis2ct
            WBAS("BASIS_IRECOUV_2CT", c, fp)

            a = (9, 10, 0, 0); WVAL("BASIS_IRECOUV_2CT", (a,), irecouv(*a, *c), fp, prec, bar)
            a = (9, 10, 0, 1); WVAL("BASIS_IRECOUV_2CT", (a,), irecouv(*a, *c), fp, prec, bar)
            a = (9, 10, 1, 0); WVAL("BASIS_IRECOUV_2CT", (a,), irecouv(*a, *c), fp, prec, bar)
            a = (9, 10, 1, 1); WVAL("BASIS_IRECOUV_2CT", (a,), irecouv(*a, *c), fp, prec, bar)

    # ==========================================================================

    if basis_izrecouv_1ct:
        with click.progressbar(length=4, label=f"{'Calc. - izrecouv 1ct':20}") as bar:

            WACT("BASIS_IZRECOUV_1CT", fp)

            c = basis1ct
            WBAS("BASIS_IZRECOUV_1CT", c, fp)

            a = (11, 10, 0, 0); WVAL("BASIS_IZRECOUV_1CT", (a,), izrecouv(*a, *c), fp, prec, bar)
            a = (11, 10, 0, 1); WVAL("BASIS_IZRECOUV_1CT", (a,), izrecouv(*a, *c), fp, prec, bar)
            a = (11, 10, 1, 0); WVAL("BASIS_IZRECOUV_1CT", (a,), izrecouv(*a, *c), fp, prec, bar)
            a = (11, 10, 1, 1); WVAL("BASIS_IZRECOUV_1CT", (a,), izrecouv(*a, *c), fp, prec, bar)

    # ==========================================================================

    if basis_izrecouv_2ct:
        with click.progressbar(length=4, label=f"{'Calc. - izrecouv 2ct':20}") as bar:

            WACT("BASIS_IZRECOUV_2CT", fp)

            c = basis2ct
            WBAS("BASIS_IZRECOUV_2CT", c, fp)

            a = (11, 10, 0, 0); WVAL("BASIS_IZRECOUV_2CT", (a,), izrecouv(*a, *c), fp, prec, bar)
            a = (11, 10, 0, 1); WVAL("BASIS_IZRECOUV_2CT", (a,), izrecouv(*a, *c), fp, prec, bar)
            a = (11, 10, 1, 0); WVAL("BASIS_IZRECOUV_2CT", (a,), izrecouv(*a, *c), fp, prec, bar)
            a = (11, 10, 1, 1); WVAL("BASIS_IZRECOUV_2CT", (a,), izrecouv(*a, *c), fp, prec, bar)

    # ==========================================================================

    if basis_hcyl_1ct:
        with click.progressbar(length=4, label=f"{'Calc. - hcyl 1ct':20}") as bar:

            WACT("BASIS_HCYL_1CT", fp)

            c = basis1ct
            WBAS("BASIS_HCYL_1CT", c, fp)

            a = (2, 0, 11, 10, 0, 0); WVAL("BASIS_HCYL_1CT", (a,), hcyl(*a, *c).doit(), fp, prec, bar)
            a = (2, 0, 11, 10, 0, 1); WVAL("BASIS_HCYL_1CT", (a,), hcyl(*a, *c).doit(), fp, prec, bar)
            a = (2, 0, 11, 10, 1, 0); WVAL("BASIS_HCYL_1CT", (a,), hcyl(*a, *c).doit(), fp, prec, bar)
            a = (2, 0, 11, 10, 1, 1); WVAL("BASIS_HCYL_1CT", (a,), hcyl(*a, *c).doit(), fp, prec, bar)

    # ==========================================================================

    if basis_hcyl_2ct:
        with click.progressbar(length=8, label=f"{'Calc. - hcyl 2ct':20}") as bar:

            WACT("BASIS_HCYL_2CT", fp)

            c = basis2ct
            WBAS("BASIS_HCYL_2CT", c, fp)

            a = (2, 0, 5, 3, 0, 0); WVAL("BASIS_HCYL_2CT", (a,), hcyl(*a, *c).doit(), fp, prec, bar)
            a = (2, 0, 5, 3, 0, 1); WVAL("BASIS_HCYL_2CT", (a,), hcyl(*a, *c).doit(), fp, prec, bar)
            a = (2, 0, 5, 3, 1, 0); WVAL("BASIS_HCYL_2CT", (a,), hcyl(*a, *c).doit(), fp, prec, bar)
            a = (2, 0, 5, 3, 1, 1); WVAL("BASIS_HCYL_2CT", (a,), hcyl(*a, *c).doit(), fp, prec, bar)
            a = (3, 0, 3, 1, 0, 0); WVAL("BASIS_HCYL_2CT", (a,), hcyl(*a, *c).doit(), fp, prec, bar)
            a = (3, 0, 3, 1, 0, 1); WVAL("BASIS_HCYL_2CT", (a,), hcyl(*a, *c).doit(), fp, prec, bar)
            a = (3, 0, 3, 1, 1, 0); WVAL("BASIS_HCYL_2CT", (a,), hcyl(*a, *c).doit(), fp, prec, bar)
            a = (3, 0, 3, 1, 1, 1); WVAL("BASIS_HCYL_2CT", (a,), hcyl(*a, *c).doit(), fp, prec, bar)

    # ==========================================================================

    if basis_overlapz_12:
        with click.progressbar(length=1, label=f"{'Calc. - overlapz 12':20}") as bar:

            WACT("BASIS_OVERLAPZ_12", fp)

            ca = basis1ct
            WBAS("BASIS_OVERLAPZ_12_A", ca, fp)

            cb = basis2ct
            WBAS("BASIS_OVERLAPZ_12_B", cb, fp)

            a = (3, 0, 0, 0)

            zparta = zpart(z, a[0], a[2], *ca)
            zpartb = zpart(z, a[1], a[3], *cb)
            WVAL("BASIS_OVERLAPZ_12", (a,), overlapz(zparta, zpartb), fp, prec, bar)

    # ==========================================================================

    if basis_overlapz_22:
        with click.progressbar(length=1, label=f"{'Calc. - overlapz 22':20}") as bar:

            WACT("BASIS_OVERLAPZ_22", fp)

            ca = basis2ct
            WBAS("BASIS_OVERLAPZ_22_A", ca, fp)

            cb = basis2ct
            WBAS("BASIS_OVERLAPZ_22_B", cb, fp)

            a = (3, 5, 0, 1)
            zparta = zpart(z, a[0], a[2], *ca)
            zpartb = zpart(z, a[1], a[3], *cb)
            WVAL("BASIS_OVERLAPZ_22", (a,), overlapz(zparta, zpartb), fp, prec, bar)

    # ==========================================================================

    if basis_overlapr_12:
        with click.progressbar(length=1, label=f"{'Calc. - overlapr 12':20}") as bar:

            WACT("BASIS_OVERLAPR_12", fp)

            ca = basis1ct
            WBAS("BASIS_OVERLAPR_12_A", ca, fp)

            cb = basis2ct
            WBAS("BASIS_OVERLAPR_12_B", cb, fp)

            a = (2, 2, 2, 1)
            rparta = rpart(rp, a[0], a[2], *ca)
            rpartb = rpart(rp, a[1], a[3], *cb)
            WVAL("BASIS_OVERLAPR_12", (a,), overlapr(rparta, rpartb), fp, prec, bar)

    # ==========================================================================

    if basis_overlapr_2X:
        with click.progressbar(length=1, label=f"{'Calc. - overlapr 2X':20}") as bar:

            WACT("BASIS_OVERLAPR_2X", fp)

            ca = basis2ct
            WBAS("BASIS_OVERLAPR_2X_A", ca, fp)

            cb = (Float("0.0", prec), Float("2.4144410198519912", prec), Float("1.0", prec))
            WBAS("BASIS_OVERLAPR_2X_B", cb, fp)

            a = (2, 2, 2, 1)
            rparta = rpart(rp, a[0], a[2], *ca)
            rpartb = rpart(rp, a[1], a[3], *cb)
            WVAL("BASIS_OVERLAPR_2X", (a,), overlapr(rparta, rpartb), fp, prec, bar)

    # ==========================================================================

    if basis_talmanz_1ct:
        with click.progressbar(length=24, label=f"{'Calc. - talmanz 1ct':20}") as bar:

            WACT("BASIS_TALMANZ_1CT", fp)

            c = basis1ct
            WBAS("BASIS_TALMANZ_1CT", c, fp)

            a = (11, 11, 0, 0,  0); WVAL("BASIS_TALMANZ_1CT", (a,), talmanz(*a, *c).doit(), fp, prec, bar)
            a = (11, 11, 0, 0,  1); WVAL("BASIS_TALMANZ_1CT", (a,), talmanz(*a, *c).doit(), fp, prec, bar)
            a = (11, 11, 0, 0,  2); WVAL("BASIS_TALMANZ_1CT", (a,), talmanz(*a, *c).doit(), fp, prec, bar)
            a = (11, 11, 0, 0,  3); WVAL("BASIS_TALMANZ_1CT", (a,), talmanz(*a, *c).doit(), fp, prec, bar)
            a = (11, 11, 0, 0,  4); WVAL("BASIS_TALMANZ_1CT", (a,), talmanz(*a, *c).doit(), fp, prec, bar)
            a = (11, 11, 0, 0,  5); WVAL("BASIS_TALMANZ_1CT", (a,), talmanz(*a, *c).doit(), fp, prec, bar)
            a = (11, 11, 0, 0,  6); WVAL("BASIS_TALMANZ_1CT", (a,), talmanz(*a, *c).doit(), fp, prec, bar)
            a = (11, 11, 0, 0,  7); WVAL("BASIS_TALMANZ_1CT", (a,), talmanz(*a, *c).doit(), fp, prec, bar)
            a = (11, 11, 0, 0,  8); WVAL("BASIS_TALMANZ_1CT", (a,), talmanz(*a, *c).doit(), fp, prec, bar)
            a = (11, 11, 0, 0,  9); WVAL("BASIS_TALMANZ_1CT", (a,), talmanz(*a, *c).doit(), fp, prec, bar)
            a = (11, 11, 0, 0, 10); WVAL("BASIS_TALMANZ_1CT", (a,), talmanz(*a, *c).doit(), fp, prec, bar)
            a = (11, 11, 0, 0, 11); WVAL("BASIS_TALMANZ_1CT", (a,), talmanz(*a, *c).doit(), fp, prec, bar)
            a = (11, 11, 0, 0, 12); WVAL("BASIS_TALMANZ_1CT", (a,), talmanz(*a, *c).doit(), fp, prec, bar)
            a = (11, 11, 0, 0, 13); WVAL("BASIS_TALMANZ_1CT", (a,), talmanz(*a, *c).doit(), fp, prec, bar)
            a = (11, 11, 0, 0, 14); WVAL("BASIS_TALMANZ_1CT", (a,), talmanz(*a, *c).doit(), fp, prec, bar)
            a = (11, 11, 0, 0, 15); WVAL("BASIS_TALMANZ_1CT", (a,), talmanz(*a, *c).doit(), fp, prec, bar)
            a = (11, 11, 0, 0, 16); WVAL("BASIS_TALMANZ_1CT", (a,), talmanz(*a, *c).doit(), fp, prec, bar)
            a = (11, 11, 0, 0, 17); WVAL("BASIS_TALMANZ_1CT", (a,), talmanz(*a, *c).doit(), fp, prec, bar)
            a = (11, 11, 0, 0, 18); WVAL("BASIS_TALMANZ_1CT", (a,), talmanz(*a, *c).doit(), fp, prec, bar)
            a = (11, 11, 0, 0, 19); WVAL("BASIS_TALMANZ_1CT", (a,), talmanz(*a, *c).doit(), fp, prec, bar)
            a = (11, 11, 0, 0, 20); WVAL("BASIS_TALMANZ_1CT", (a,), talmanz(*a, *c).doit(), fp, prec, bar)
            a = (11, 11, 0, 0, 21); WVAL("BASIS_TALMANZ_1CT", (a,), talmanz(*a, *c).doit(), fp, prec, bar)
            a = (11, 11, 0, 0, 22); WVAL("BASIS_TALMANZ_1CT", (a,), talmanz(*a, *c).doit(), fp, prec, bar)
            a = (11, 11, 0, 0, 23); WVAL("BASIS_TALMANZ_1CT", (a,), talmanz(*a, *c).doit(), fp, prec, bar)

    # ==========================================================================

    if basis_talmanz_2ct:
        with click.progressbar(length=16, label=f"{'Calc. - talmanz 2ct':20}") as bar:

            WACT("BASIS_TALMANZ_2CT", fp)

            c = basis2ct
            WBAS("BASIS_TALMANZ_2CT", c, fp)

            a = (9, 8, 0, 0, 0); WVAL("BASIS_TALMANZ_2CT", (a,), talmanz(*a, *c).doit(), fp, prec, bar)
            a = (9, 8, 0, 0, 1); WVAL("BASIS_TALMANZ_2CT", (a,), talmanz(*a, *c).doit(), fp, prec, bar)
            a = (9, 8, 0, 0, 2); WVAL("BASIS_TALMANZ_2CT", (a,), talmanz(*a, *c).doit(), fp, prec, bar)
            a = (9, 8, 0, 0, 3); WVAL("BASIS_TALMANZ_2CT", (a,), talmanz(*a, *c).doit(), fp, prec, bar)
            a = (9, 8, 0, 0, 4); WVAL("BASIS_TALMANZ_2CT", (a,), talmanz(*a, *c).doit(), fp, prec, bar)
            a = (9, 8, 0, 0, 5); WVAL("BASIS_TALMANZ_2CT", (a,), talmanz(*a, *c).doit(), fp, prec, bar)
            a = (9, 8, 0, 0, 6); WVAL("BASIS_TALMANZ_2CT", (a,), talmanz(*a, *c).doit(), fp, prec, bar)
            a = (9, 8, 0, 0, 7); WVAL("BASIS_TALMANZ_2CT", (a,), talmanz(*a, *c).doit(), fp, prec, bar)
            a = (9, 8, 0, 0, 3); WVAL("BASIS_TALMANZ_2CT", (a,), talmanz(*a, *c).doit(), fp, prec, bar)
            a = (9, 8, 0, 1, 3); WVAL("BASIS_TALMANZ_2CT", (a,), talmanz(*a, *c).doit(), fp, prec, bar)
            a = (9, 8, 0, 0, 3); WVAL("BASIS_TALMANZ_2CT", (a,), talmanz(*a, *c).doit(), fp, prec, bar)
            a = (9, 8, 0, 1, 3); WVAL("BASIS_TALMANZ_2CT", (a,), talmanz(*a, *c).doit(), fp, prec, bar)
            a = (9, 8, 1, 0, 3); WVAL("BASIS_TALMANZ_2CT", (a,), talmanz(*a, *c).doit(), fp, prec, bar)
            a = (9, 8, 1, 1, 3); WVAL("BASIS_TALMANZ_2CT", (a,), talmanz(*a, *c).doit(), fp, prec, bar)
            a = (9, 8, 1, 0, 3); WVAL("BASIS_TALMANZ_2CT", (a,), talmanz(*a, *c).doit(), fp, prec, bar)
            a = (9, 8, 1, 1, 3); WVAL("BASIS_TALMANZ_2CT", (a,), talmanz(*a, *c).doit(), fp, prec, bar)

    # ==========================================================================

    if fieldcentraldirect_calcgaussianintegralz_1ct:
        with click.progressbar(length=41, label=f"{'Calc. - gIntz 1ct':20}") as bar:

            WACT("FIELDCENTRALDIRECT_CALCGAUSSIANINTEGRALZ_1CT", fp)

            c = basis1ct
            WBAS("FIELDCENTRALDIRECT_CALCGAUSSIANINTEGRALZ_1CT", c, fp)

            a = (0, 0, 0, 0,  0); b = (1.0,); WVAL("FIELDCENTRALDIRECT_CALCGAUSSIANINTEGRALZ_1CT", (a, b), gaussianIntegralz(*a, *b, *c).doit(), fp, prec, bar)
            a = (0, 0, 0, 0,  1); b = (1.0,); WVAL("FIELDCENTRALDIRECT_CALCGAUSSIANINTEGRALZ_1CT", (a, b), gaussianIntegralz(*a, *b, *c).doit(), fp, prec, bar)
            a = (0, 0, 0, 0,  2); b = (1.0,); WVAL("FIELDCENTRALDIRECT_CALCGAUSSIANINTEGRALZ_1CT", (a, b), gaussianIntegralz(*a, *b, *c).doit(), fp, prec, bar)
            a = (0, 0, 0, 0,  3); b = (1.0,); WVAL("FIELDCENTRALDIRECT_CALCGAUSSIANINTEGRALZ_1CT", (a, b), gaussianIntegralz(*a, *b, *c).doit(), fp, prec, bar)
            a = (0, 0, 0, 0,  4); b = (1.0,); WVAL("FIELDCENTRALDIRECT_CALCGAUSSIANINTEGRALZ_1CT", (a, b), gaussianIntegralz(*a, *b, *c).doit(), fp, prec, bar)
            a = (0, 0, 0, 0,  5); b = (1.0,); WVAL("FIELDCENTRALDIRECT_CALCGAUSSIANINTEGRALZ_1CT", (a, b), gaussianIntegralz(*a, *b, *c).doit(), fp, prec, bar)
            a = (0, 0, 0, 0,  6); b = (1.0,); WVAL("FIELDCENTRALDIRECT_CALCGAUSSIANINTEGRALZ_1CT", (a, b), gaussianIntegralz(*a, *b, *c).doit(), fp, prec, bar)
            a = (0, 0, 0, 0,  7); b = (1.0,); WVAL("FIELDCENTRALDIRECT_CALCGAUSSIANINTEGRALZ_1CT", (a, b), gaussianIntegralz(*a, *b, *c).doit(), fp, prec, bar)
            a = (0, 0, 0, 0,  8); b = (1.0,); WVAL("FIELDCENTRALDIRECT_CALCGAUSSIANINTEGRALZ_1CT", (a, b), gaussianIntegralz(*a, *b, *c).doit(), fp, prec, bar)
            a = (0, 0, 0, 0,  9); b = (1.0,); WVAL("FIELDCENTRALDIRECT_CALCGAUSSIANINTEGRALZ_1CT", (a, b), gaussianIntegralz(*a, *b, *c).doit(), fp, prec, bar)
            a = (0, 0, 0, 0, 10); b = (1.0,); WVAL("FIELDCENTRALDIRECT_CALCGAUSSIANINTEGRALZ_1CT", (a, b), gaussianIntegralz(*a, *b, *c).doit(), fp, prec, bar)
            a = (0, 0, 0, 0, 11); b = (1.0,); WVAL("FIELDCENTRALDIRECT_CALCGAUSSIANINTEGRALZ_1CT", (a, b), gaussianIntegralz(*a, *b, *c).doit(), fp, prec, bar)
            a = (0, 0, 0, 0, 12); b = (1.0,); WVAL("FIELDCENTRALDIRECT_CALCGAUSSIANINTEGRALZ_1CT", (a, b), gaussianIntegralz(*a, *b, *c).doit(), fp, prec, bar)
            a = (0, 0, 0, 0, 13); b = (1.0,); WVAL("FIELDCENTRALDIRECT_CALCGAUSSIANINTEGRALZ_1CT", (a, b), gaussianIntegralz(*a, *b, *c).doit(), fp, prec, bar)
            a = (0, 0, 0, 0, 14); b = (1.0,); WVAL("FIELDCENTRALDIRECT_CALCGAUSSIANINTEGRALZ_1CT", (a, b), gaussianIntegralz(*a, *b, *c).doit(), fp, prec, bar)
            a = (0, 0, 0, 0, 15); b = (1.0,); WVAL("FIELDCENTRALDIRECT_CALCGAUSSIANINTEGRALZ_1CT", (a, b), gaussianIntegralz(*a, *b, *c).doit(), fp, prec, bar)
            a = (0, 0, 0, 0, 16); b = (1.0,); WVAL("FIELDCENTRALDIRECT_CALCGAUSSIANINTEGRALZ_1CT", (a, b), gaussianIntegralz(*a, *b, *c).doit(), fp, prec, bar)
            a = (0, 0, 0, 0, 17); b = (1.0,); WVAL("FIELDCENTRALDIRECT_CALCGAUSSIANINTEGRALZ_1CT", (a, b), gaussianIntegralz(*a, *b, *c).doit(), fp, prec, bar)
            a = (0, 0, 0, 0, 18); b = (1.0,); WVAL("FIELDCENTRALDIRECT_CALCGAUSSIANINTEGRALZ_1CT", (a, b), gaussianIntegralz(*a, *b, *c).doit(), fp, prec, bar)
            a = (0, 0, 0, 0, 19); b = (1.0,); WVAL("FIELDCENTRALDIRECT_CALCGAUSSIANINTEGRALZ_1CT", (a, b), gaussianIntegralz(*a, *b, *c).doit(), fp, prec, bar)
            a = (0, 0, 0, 0, 20); b = (1.0,); WVAL("FIELDCENTRALDIRECT_CALCGAUSSIANINTEGRALZ_1CT", (a, b), gaussianIntegralz(*a, *b, *c).doit(), fp, prec, bar)
            a = (0, 0, 0, 0, 21); b = (1.0,); WVAL("FIELDCENTRALDIRECT_CALCGAUSSIANINTEGRALZ_1CT", (a, b), gaussianIntegralz(*a, *b, *c).doit(), fp, prec, bar)
            a = (0, 0, 0, 0, 22); b = (1.0,); WVAL("FIELDCENTRALDIRECT_CALCGAUSSIANINTEGRALZ_1CT", (a, b), gaussianIntegralz(*a, *b, *c).doit(), fp, prec, bar)
            a = (0, 0, 0, 0, 23); b = (1.0,); WVAL("FIELDCENTRALDIRECT_CALCGAUSSIANINTEGRALZ_1CT", (a, b), gaussianIntegralz(*a, *b, *c).doit(), fp, prec, bar)
            a = (0, 0, 0, 0, 24); b = (1.0,); WVAL("FIELDCENTRALDIRECT_CALCGAUSSIANINTEGRALZ_1CT", (a, b), gaussianIntegralz(*a, *b, *c).doit(), fp, prec, bar)
            a = (0, 0, 0, 0, 25); b = (1.0,); WVAL("FIELDCENTRALDIRECT_CALCGAUSSIANINTEGRALZ_1CT", (a, b), gaussianIntegralz(*a, *b, *c).doit(), fp, prec, bar)
            a = (0, 0, 0, 0, 26); b = (1.0,); WVAL("FIELDCENTRALDIRECT_CALCGAUSSIANINTEGRALZ_1CT", (a, b), gaussianIntegralz(*a, *b, *c).doit(), fp, prec, bar)
            a = (0, 0, 0, 0, 27); b = (1.0,); WVAL("FIELDCENTRALDIRECT_CALCGAUSSIANINTEGRALZ_1CT", (a, b), gaussianIntegralz(*a, *b, *c).doit(), fp, prec, bar)
            a = (0, 0, 0, 0, 28); b = (1.0,); WVAL("FIELDCENTRALDIRECT_CALCGAUSSIANINTEGRALZ_1CT", (a, b), gaussianIntegralz(*a, *b, *c).doit(), fp, prec, bar)
            a = (0, 0, 0, 0, 29); b = (1.0,); WVAL("FIELDCENTRALDIRECT_CALCGAUSSIANINTEGRALZ_1CT", (a, b), gaussianIntegralz(*a, *b, *c).doit(), fp, prec, bar)
            a = (0, 0, 0, 0, 30); b = (1.0,); WVAL("FIELDCENTRALDIRECT_CALCGAUSSIANINTEGRALZ_1CT", (a, b), gaussianIntegralz(*a, *b, *c).doit(), fp, prec, bar)
            a = (0, 0, 0, 0, 31); b = (1.0,); WVAL("FIELDCENTRALDIRECT_CALCGAUSSIANINTEGRALZ_1CT", (a, b), gaussianIntegralz(*a, *b, *c).doit(), fp, prec, bar)
            a = (0, 0, 0, 0, 32); b = (1.0,); WVAL("FIELDCENTRALDIRECT_CALCGAUSSIANINTEGRALZ_1CT", (a, b), gaussianIntegralz(*a, *b, *c).doit(), fp, prec, bar)
            a = (0, 0, 0, 0, 33); b = (1.0,); WVAL("FIELDCENTRALDIRECT_CALCGAUSSIANINTEGRALZ_1CT", (a, b), gaussianIntegralz(*a, *b, *c).doit(), fp, prec, bar)
            a = (0, 0, 0, 0, 34); b = (1.0,); WVAL("FIELDCENTRALDIRECT_CALCGAUSSIANINTEGRALZ_1CT", (a, b), gaussianIntegralz(*a, *b, *c).doit(), fp, prec, bar)
            a = (0, 0, 0, 0, 35); b = (1.0,); WVAL("FIELDCENTRALDIRECT_CALCGAUSSIANINTEGRALZ_1CT", (a, b), gaussianIntegralz(*a, *b, *c).doit(), fp, prec, bar)
            a = (0, 0, 0, 0, 36); b = (1.0,); WVAL("FIELDCENTRALDIRECT_CALCGAUSSIANINTEGRALZ_1CT", (a, b), gaussianIntegralz(*a, *b, *c).doit(), fp, prec, bar)
            a = (0, 0, 0, 0, 37); b = (1.0,); WVAL("FIELDCENTRALDIRECT_CALCGAUSSIANINTEGRALZ_1CT", (a, b), gaussianIntegralz(*a, *b, *c).doit(), fp, prec, bar)
            a = (0, 0, 0, 0, 38); b = (1.0,); WVAL("FIELDCENTRALDIRECT_CALCGAUSSIANINTEGRALZ_1CT", (a, b), gaussianIntegralz(*a, *b, *c).doit(), fp, prec, bar)
            a = (0, 0, 0, 0, 39); b = (1.0,); WVAL("FIELDCENTRALDIRECT_CALCGAUSSIANINTEGRALZ_1CT", (a, b), gaussianIntegralz(*a, *b, *c).doit(), fp, prec, bar)
            a = (0, 0, 0, 0, 40); b = (1.0,); WVAL("FIELDCENTRALDIRECT_CALCGAUSSIANINTEGRALZ_1CT", (a, b), gaussianIntegralz(*a, *b, *c).doit(), fp, prec, bar)

    # ==========================================================================

    if fieldcentraldirect_calcgaussianintegralz_2ct:
        with click.progressbar(length=41, label=f"{'Calc. - gIntz 2ct':20}") as bar:

            WACT("FIELDCENTRALDIRECT_CALCGAUSSIANINTEGRALZ_2CT", fp)

            c = basis2ct
            WBAS("FIELDCENTRALDIRECT_CALCGAUSSIANINTEGRALZ_2CT", c, fp)

            a = (0, 0, 0, 0,  0); b = (1.0,); WVAL("FIELDCENTRALDIRECT_CALCGAUSSIANINTEGRALZ_2CT", (a, b), gaussianIntegralz(*a, *b, *c).doit(), fp, prec, bar)
            a = (0, 0, 0, 1,  1); b = (1.0,); WVAL("FIELDCENTRALDIRECT_CALCGAUSSIANINTEGRALZ_2CT", (a, b), gaussianIntegralz(*a, *b, *c).doit(), fp, prec, bar)
            a = (0, 0, 1, 0,  2); b = (1.0,); WVAL("FIELDCENTRALDIRECT_CALCGAUSSIANINTEGRALZ_2CT", (a, b), gaussianIntegralz(*a, *b, *c).doit(), fp, prec, bar)
            a = (0, 0, 1, 1,  3); b = (1.0,); WVAL("FIELDCENTRALDIRECT_CALCGAUSSIANINTEGRALZ_2CT", (a, b), gaussianIntegralz(*a, *b, *c).doit(), fp, prec, bar)
            a = (0, 1, 0, 0,  4); b = (1.0,); WVAL("FIELDCENTRALDIRECT_CALCGAUSSIANINTEGRALZ_2CT", (a, b), gaussianIntegralz(*a, *b, *c).doit(), fp, prec, bar)
            a = (0, 1, 0, 1,  5); b = (1.0,); WVAL("FIELDCENTRALDIRECT_CALCGAUSSIANINTEGRALZ_2CT", (a, b), gaussianIntegralz(*a, *b, *c).doit(), fp, prec, bar)
            a = (0, 1, 1, 0,  6); b = (1.0,); WVAL("FIELDCENTRALDIRECT_CALCGAUSSIANINTEGRALZ_2CT", (a, b), gaussianIntegralz(*a, *b, *c).doit(), fp, prec, bar)
            a = (0, 1, 1, 1,  7); b = (1.0,); WVAL("FIELDCENTRALDIRECT_CALCGAUSSIANINTEGRALZ_2CT", (a, b), gaussianIntegralz(*a, *b, *c).doit(), fp, prec, bar)
            a = (1, 0, 0, 0,  8); b = (1.0,); WVAL("FIELDCENTRALDIRECT_CALCGAUSSIANINTEGRALZ_2CT", (a, b), gaussianIntegralz(*a, *b, *c).doit(), fp, prec, bar)
            a = (1, 0, 0, 1,  9); b = (1.0,); WVAL("FIELDCENTRALDIRECT_CALCGAUSSIANINTEGRALZ_2CT", (a, b), gaussianIntegralz(*a, *b, *c).doit(), fp, prec, bar)
            a = (1, 0, 1, 0, 10); b = (1.0,); WVAL("FIELDCENTRALDIRECT_CALCGAUSSIANINTEGRALZ_2CT", (a, b), gaussianIntegralz(*a, *b, *c).doit(), fp, prec, bar)
            a = (1, 0, 1, 1, 11); b = (1.0,); WVAL("FIELDCENTRALDIRECT_CALCGAUSSIANINTEGRALZ_2CT", (a, b), gaussianIntegralz(*a, *b, *c).doit(), fp, prec, bar)
            a = (1, 1, 0, 0, 12); b = (1.0,); WVAL("FIELDCENTRALDIRECT_CALCGAUSSIANINTEGRALZ_2CT", (a, b), gaussianIntegralz(*a, *b, *c).doit(), fp, prec, bar)
            a = (1, 1, 0, 1, 13); b = (1.0,); WVAL("FIELDCENTRALDIRECT_CALCGAUSSIANINTEGRALZ_2CT", (a, b), gaussianIntegralz(*a, *b, *c).doit(), fp, prec, bar)
            a = (1, 1, 1, 0, 14); b = (1.0,); WVAL("FIELDCENTRALDIRECT_CALCGAUSSIANINTEGRALZ_2CT", (a, b), gaussianIntegralz(*a, *b, *c).doit(), fp, prec, bar)
            a = (1, 1, 1, 1, 15); b = (1.0,); WVAL("FIELDCENTRALDIRECT_CALCGAUSSIANINTEGRALZ_2CT", (a, b), gaussianIntegralz(*a, *b, *c).doit(), fp, prec, bar)
            a = (0, 0, 0, 0, 16); b = (1.0,); WVAL("FIELDCENTRALDIRECT_CALCGAUSSIANINTEGRALZ_2CT", (a, b), gaussianIntegralz(*a, *b, *c).doit(), fp, prec, bar)
            a = (0, 0, 0, 1, 17); b = (1.0,); WVAL("FIELDCENTRALDIRECT_CALCGAUSSIANINTEGRALZ_2CT", (a, b), gaussianIntegralz(*a, *b, *c).doit(), fp, prec, bar)
            a = (0, 0, 1, 0, 18); b = (1.0,); WVAL("FIELDCENTRALDIRECT_CALCGAUSSIANINTEGRALZ_2CT", (a, b), gaussianIntegralz(*a, *b, *c).doit(), fp, prec, bar)
            a = (0, 0, 1, 1, 19); b = (1.0,); WVAL("FIELDCENTRALDIRECT_CALCGAUSSIANINTEGRALZ_2CT", (a, b), gaussianIntegralz(*a, *b, *c).doit(), fp, prec, bar)
            a = (0, 1, 0, 0, 20); b = (1.0,); WVAL("FIELDCENTRALDIRECT_CALCGAUSSIANINTEGRALZ_2CT", (a, b), gaussianIntegralz(*a, *b, *c).doit(), fp, prec, bar)
            a = (0, 1, 0, 1, 21); b = (1.0,); WVAL("FIELDCENTRALDIRECT_CALCGAUSSIANINTEGRALZ_2CT", (a, b), gaussianIntegralz(*a, *b, *c).doit(), fp, prec, bar)
            a = (0, 1, 1, 0, 22); b = (1.0,); WVAL("FIELDCENTRALDIRECT_CALCGAUSSIANINTEGRALZ_2CT", (a, b), gaussianIntegralz(*a, *b, *c).doit(), fp, prec, bar)
            a = (0, 1, 1, 1, 23); b = (1.0,); WVAL("FIELDCENTRALDIRECT_CALCGAUSSIANINTEGRALZ_2CT", (a, b), gaussianIntegralz(*a, *b, *c).doit(), fp, prec, bar)
            a = (1, 0, 0, 0, 24); b = (1.0,); WVAL("FIELDCENTRALDIRECT_CALCGAUSSIANINTEGRALZ_2CT", (a, b), gaussianIntegralz(*a, *b, *c).doit(), fp, prec, bar)
            a = (1, 0, 0, 1, 25); b = (1.0,); WVAL("FIELDCENTRALDIRECT_CALCGAUSSIANINTEGRALZ_2CT", (a, b), gaussianIntegralz(*a, *b, *c).doit(), fp, prec, bar)
            a = (1, 0, 1, 0, 26); b = (1.0,); WVAL("FIELDCENTRALDIRECT_CALCGAUSSIANINTEGRALZ_2CT", (a, b), gaussianIntegralz(*a, *b, *c).doit(), fp, prec, bar)
            a = (1, 0, 1, 1, 27); b = (1.0,); WVAL("FIELDCENTRALDIRECT_CALCGAUSSIANINTEGRALZ_2CT", (a, b), gaussianIntegralz(*a, *b, *c).doit(), fp, prec, bar)
            a = (1, 1, 0, 0, 28); b = (1.0,); WVAL("FIELDCENTRALDIRECT_CALCGAUSSIANINTEGRALZ_2CT", (a, b), gaussianIntegralz(*a, *b, *c).doit(), fp, prec, bar)
            a = (1, 1, 0, 1, 29); b = (1.0,); WVAL("FIELDCENTRALDIRECT_CALCGAUSSIANINTEGRALZ_2CT", (a, b), gaussianIntegralz(*a, *b, *c).doit(), fp, prec, bar)
            a = (1, 1, 1, 0, 30); b = (1.0,); WVAL("FIELDCENTRALDIRECT_CALCGAUSSIANINTEGRALZ_2CT", (a, b), gaussianIntegralz(*a, *b, *c).doit(), fp, prec, bar)
            a = (1, 1, 1, 1, 31); b = (1.0,); WVAL("FIELDCENTRALDIRECT_CALCGAUSSIANINTEGRALZ_2CT", (a, b), gaussianIntegralz(*a, *b, *c).doit(), fp, prec, bar)
            a = (0, 0, 0, 0, 32); b = (1.0,); WVAL("FIELDCENTRALDIRECT_CALCGAUSSIANINTEGRALZ_2CT", (a, b), gaussianIntegralz(*a, *b, *c).doit(), fp, prec, bar)
            a = (0, 0, 0, 1, 33); b = (1.0,); WVAL("FIELDCENTRALDIRECT_CALCGAUSSIANINTEGRALZ_2CT", (a, b), gaussianIntegralz(*a, *b, *c).doit(), fp, prec, bar)
            a = (0, 0, 1, 0, 34); b = (1.0,); WVAL("FIELDCENTRALDIRECT_CALCGAUSSIANINTEGRALZ_2CT", (a, b), gaussianIntegralz(*a, *b, *c).doit(), fp, prec, bar)
            a = (0, 0, 1, 1, 35); b = (1.0,); WVAL("FIELDCENTRALDIRECT_CALCGAUSSIANINTEGRALZ_2CT", (a, b), gaussianIntegralz(*a, *b, *c).doit(), fp, prec, bar)
            a = (0, 1, 0, 0, 36); b = (1.0,); WVAL("FIELDCENTRALDIRECT_CALCGAUSSIANINTEGRALZ_2CT", (a, b), gaussianIntegralz(*a, *b, *c).doit(), fp, prec, bar)
            a = (0, 1, 0, 1, 37); b = (1.0,); WVAL("FIELDCENTRALDIRECT_CALCGAUSSIANINTEGRALZ_2CT", (a, b), gaussianIntegralz(*a, *b, *c).doit(), fp, prec, bar)
            a = (0, 1, 1, 0, 38); b = (1.0,); WVAL("FIELDCENTRALDIRECT_CALCGAUSSIANINTEGRALZ_2CT", (a, b), gaussianIntegralz(*a, *b, *c).doit(), fp, prec, bar)
            a = (0, 1, 1, 1, 39); b = (1.0,); WVAL("FIELDCENTRALDIRECT_CALCGAUSSIANINTEGRALZ_2CT", (a, b), gaussianIntegralz(*a, *b, *c).doit(), fp, prec, bar)
            a = (1, 0, 0, 0, 40); b = (1.0,); WVAL("FIELDCENTRALDIRECT_CALCGAUSSIANINTEGRALZ_2CT", (a, b), gaussianIntegralz(*a, *b, *c).doit(), fp, prec, bar)

    # ==========================================================================

    if multipole_operators_znrecouv_1ct:
        with click.progressbar(length=49, label=f"{'Calc. - znrecouv 1ct':20}") as bar:

            WACT("MULTIPOLE_OPERATORS_ZNRECOUV_1CT", fp)

            c = basis1ct
            WBAS("MULTIPOLE_OPERATORS_ZNRECOUV_1CT", c, fp)

            for n in range(7):
                a = (n, 0, 0, 0, 0); WVAL("MULTIPOLE_OPERATORS_ZNRECOUV_1CT", (a,), znrecouv(*a, *c), fp, prec, bar)
                a = (n, 1, 1, 0, 0); WVAL("MULTIPOLE_OPERATORS_ZNRECOUV_1CT", (a,), znrecouv(*a, *c), fp, prec, bar)
                a = (n, 1, 3, 0, 0); WVAL("MULTIPOLE_OPERATORS_ZNRECOUV_1CT", (a,), znrecouv(*a, *c), fp, prec, bar)
                a = (n, 2, 4, 0, 0); WVAL("MULTIPOLE_OPERATORS_ZNRECOUV_1CT", (a,), znrecouv(*a, *c), fp, prec, bar)
                a = (n, 3, 4, 0, 0); WVAL("MULTIPOLE_OPERATORS_ZNRECOUV_1CT", (a,), znrecouv(*a, *c), fp, prec, bar)
                a = (n, 4, 3, 0, 0); WVAL("MULTIPOLE_OPERATORS_ZNRECOUV_1CT", (a,), znrecouv(*a, *c), fp, prec, bar)
                a = (n, 6, 7, 0, 0); WVAL("MULTIPOLE_OPERATORS_ZNRECOUV_1CT", (a,), znrecouv(*a, *c), fp, prec, bar)

    # ==========================================================================

    if multipole_operators_znrecouv_2ct:
        with click.progressbar(length=49, label=f"{'Calc. - znrecouv 2ct':20}") as bar:

            WACT("MULTIPOLE_OPERATORS_ZNRECOUV_2CT", fp)

            c = basis2ct
            WBAS("MULTIPOLE_OPERATORS_ZNRECOUV_2CT", c, fp)

            for n in range(7):
                a = (n, 0, 0, 0, 0); WVAL("MULTIPOLE_OPERATORS_ZNRECOUV_2CT", (a,), znrecouv(*a, *c), fp, prec, bar)
                a = (n, 1, 1, 0, 1); WVAL("MULTIPOLE_OPERATORS_ZNRECOUV_2CT", (a,), znrecouv(*a, *c), fp, prec, bar)
                a = (n, 1, 3, 1, 0); WVAL("MULTIPOLE_OPERATORS_ZNRECOUV_2CT", (a,), znrecouv(*a, *c), fp, prec, bar)
                a = (n, 2, 4, 1, 1); WVAL("MULTIPOLE_OPERATORS_ZNRECOUV_2CT", (a,), znrecouv(*a, *c), fp, prec, bar)
                a = (n, 3, 4, 0, 1); WVAL("MULTIPOLE_OPERATORS_ZNRECOUV_2CT", (a,), znrecouv(*a, *c), fp, prec, bar)
                a = (n, 4, 3, 0, 0); WVAL("MULTIPOLE_OPERATORS_ZNRECOUV_2CT", (a,), znrecouv(*a, *c), fp, prec, bar)
                a = (n, 6, 7, 1, 0); WVAL("MULTIPOLE_OPERATORS_ZNRECOUV_2CT", (a,), znrecouv(*a, *c), fp, prec, bar)

    # ==========================================================================

    if multipole_operators_rnrecouv:
        with click.progressbar(length=36, label=f"{'Calc. - rnrecouv':20}") as bar:

            WACT("MULTIPOLE_OPERATORS_RNRECOUV", fp)

            c = basis1ct
            WBAS("MULTIPOLE_OPERATORS_RNRECOUV", c, fp)

            for n in range(1, 7):
                a = (n, 0, 0, 1, 0); WVAL("MULTIPOLE_OPERATORS_RNRECOUV", (a,), rnrecouv(*a, *c), fp, prec, bar)
                a = (n, 1, 1, 3, 0); WVAL("MULTIPOLE_OPERATORS_RNRECOUV", (a,), rnrecouv(*a, *c), fp, prec, bar)
                a = (n, 1, 1, 0, 1); WVAL("MULTIPOLE_OPERATORS_RNRECOUV", (a,), rnrecouv(*a, *c), fp, prec, bar)
                a = (n, 2, 2, 1, 3); WVAL("MULTIPOLE_OPERATORS_RNRECOUV", (a,), rnrecouv(*a, *c), fp, prec, bar)
                a = (n, 3, 3, 2, 2); WVAL("MULTIPOLE_OPERATORS_RNRECOUV", (a,), rnrecouv(*a, *c), fp, prec, bar)
                a = (n, 4, 4, 0, 1); WVAL("MULTIPOLE_OPERATORS_RNRECOUV", (a,), rnrecouv(*a, *c), fp, prec, bar)
                a = (n, 6, 6, 2, 2); WVAL("MULTIPOLE_OPERATORS_RNRECOUV", (a,), rnrecouv(*a, *c), fp, prec, bar)

    # ==========================================================================

    if multipole_operators_rnodd:
        with click.progressbar(length=16, label=f"{'Calc. - rnodd':20}") as bar:

            WACT("MULTIPOLE_OPERATORS_RNODD", fp)

            c = basis1ct
            WBAS("MULTIPOLE_OPERATORS_RNODD", c, fp)

            for n in range(2, 4):
                a = (n, 0, 0, 0, 0); WVAL("MULTIPOLE_OPERATORS_RNODD", (a,), rnodd(*a, *c), fp, prec, bar)
                a = (n, 0, 1, 0, 0); WVAL("MULTIPOLE_OPERATORS_RNODD", (a,), rnodd(*a, *c), fp, prec, bar)
                a = (n, 0, 2, 0, 0); WVAL("MULTIPOLE_OPERATORS_RNODD", (a,), rnodd(*a, *c), fp, prec, bar)
                a = (n, 0, 3, 1, 1); WVAL("MULTIPOLE_OPERATORS_RNODD", (a,), rnodd(*a, *c), fp, prec, bar)
                a = (n, 0, 4, 1, 1); WVAL("MULTIPOLE_OPERATORS_RNODD", (a,), rnodd(*a, *c), fp, prec, bar)
                a = (n, 0, 5, 1, 2); WVAL("MULTIPOLE_OPERATORS_RNODD", (a,), rnodd(*a, *c), fp, prec, bar)
                a = (n, 2, 3, 2, 2); WVAL("MULTIPOLE_OPERATORS_RNODD", (a,), rnodd(*a, *c), fp, prec, bar)
                a = (n, 2, 2, 3, 3); WVAL("MULTIPOLE_OPERATORS_RNODD", (a,), rnodd(*a, *c), fp, prec, bar)

    # ==========================================================================

    if multipole_operators_qlmho_1ct:
        with click.progressbar(length=8*len(valList), label=f"{'Calc. - qlmHO 1ct':20}") as bar:

            WACT("MULTIPOLE_OPERATORS_QLMHO_1CT", fp)

            c = basis1ct
            WBAS("MULTIPOLE_OPERATORS_QLMHO_1CT", c, fp)

            for l, m in valList:
                ma = 2
                mb = ma - m
                a = (ma, 1, 1, 0, 0); b = (mb, 1, 0, 0, 0); WVAL("MULTIPOLE_OPERATORS_QLMHO_1CT", ((l, m), a, b), qlmHO[(l, m)](*a, *b, *c).doit(), fp, prec, bar)
                a = (ma, 1, 2, 0, 0); b = (mb, 0, 1, 0, 0); WVAL("MULTIPOLE_OPERATORS_QLMHO_1CT", ((l, m), a, b), qlmHO[(l, m)](*a, *b, *c).doit(), fp, prec, bar)
                a = (ma, 2, 0, 0, 0); b = (mb, 2, 0, 0, 0); WVAL("MULTIPOLE_OPERATORS_QLMHO_1CT", ((l, m), a, b), qlmHO[(l, m)](*a, *b, *c).doit(), fp, prec, bar)
                a = (ma, 2, 0, 0, 0); b = (mb, 1, 2, 0, 0); WVAL("MULTIPOLE_OPERATORS_QLMHO_1CT", ((l, m), a, b), qlmHO[(l, m)](*a, *b, *c).doit(), fp, prec, bar)
                a = (ma, 1, 2, 0, 0); b = (mb, 1, 1, 0, 0); WVAL("MULTIPOLE_OPERATORS_QLMHO_1CT", ((l, m), a, b), qlmHO[(l, m)](*a, *b, *c).doit(), fp, prec, bar)
                a = (ma, 2, 2, 0, 0); b = (mb, 2, 0, 0, 0); WVAL("MULTIPOLE_OPERATORS_QLMHO_1CT", ((l, m), a, b), qlmHO[(l, m)](*a, *b, *c).doit(), fp, prec, bar)
                a = (ma, 2, 0, 0, 0); b = (mb, 2, 0, 0, 0); WVAL("MULTIPOLE_OPERATORS_QLMHO_1CT", ((l, m), a, b), qlmHO[(l, m)](*a, *b, *c).doit(), fp, prec, bar)
                a = (ma, 2, 1, 0, 0); b = (mb, 0, 2, 0, 0); WVAL("MULTIPOLE_OPERATORS_QLMHO_1CT", ((l, m), a, b), qlmHO[(l, m)](*a, *b, *c).doit(), fp, prec, bar)

    # ==========================================================================

    if multipole_operators_qlmho_2ct:
        with click.progressbar(length=8*len(valList), label=f"{'Calc. - qlmHO 2ct':20}") as bar:

            WACT("MULTIPOLE_OPERATORS_QLMHO_2CT", fp)

            c = basis2ct
            WBAS("MULTIPOLE_OPERATORS_QLMHO_2CT", c, fp)

            for l, m in valList:
                ma = 2
                mb = ma - m
                a = (ma, 1, 1, 0, 0); b = (mb, 1, 0, 0, 0); WVAL("MULTIPOLE_OPERATORS_QLMHO_2CT", ((l, m), a, b), qlmHO[(l, m)](*a, *b, *c).doit(), fp, prec, bar)
                a = (ma, 1, 2, 0, 0); b = (mb, 0, 1, 0, 0); WVAL("MULTIPOLE_OPERATORS_QLMHO_2CT", ((l, m), a, b), qlmHO[(l, m)](*a, *b, *c).doit(), fp, prec, bar)
                a = (ma, 2, 0, 0, 0); b = (mb, 2, 0, 0, 0); WVAL("MULTIPOLE_OPERATORS_QLMHO_2CT", ((l, m), a, b), qlmHO[(l, m)](*a, *b, *c).doit(), fp, prec, bar)
                a = (ma, 2, 0, 0, 0); b = (mb, 1, 2, 0, 0); WVAL("MULTIPOLE_OPERATORS_QLMHO_2CT", ((l, m), a, b), qlmHO[(l, m)](*a, *b, *c).doit(), fp, prec, bar)
                a = (ma, 1, 2, 0, 0); b = (mb, 1, 1, 0, 0); WVAL("MULTIPOLE_OPERATORS_QLMHO_2CT", ((l, m), a, b), qlmHO[(l, m)](*a, *b, *c).doit(), fp, prec, bar)
                a = (ma, 2, 2, 0, 0); b = (mb, 2, 0, 0, 0); WVAL("MULTIPOLE_OPERATORS_QLMHO_2CT", ((l, m), a, b), qlmHO[(l, m)](*a, *b, *c).doit(), fp, prec, bar)
                a = (ma, 2, 0, 0, 0); b = (mb, 2, 0, 0, 0); WVAL("MULTIPOLE_OPERATORS_QLMHO_2CT", ((l, m), a, b), qlmHO[(l, m)](*a, *b, *c).doit(), fp, prec, bar)
                a = (ma, 2, 1, 0, 0); b = (mb, 0, 2, 0, 0); WVAL("MULTIPOLE_OPERATORS_QLMHO_2CT", ((l, m), a, b), qlmHO[(l, m)](*a, *b, *c).doit(), fp, prec, bar)

# =============================================================1================
# ==============================================================================
# ==============================================================================


if __name__ == "__main__":

    main()
