//==============================================================================
// HFB3
// Copyright CEA, DAM F-91297 Arpajon, France
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//==============================================================================

#include "io_berger.h"
#include "solver_hfb_broyden.h"
#include "tools.h"
#include <fstream>
#include <stdlib.h>

/** \file
 *  \brief Methods of the IOberger class.
 */

/*                            \
 *                             t   .
 *                              ?  (
 *                              l  `,  .
 *                              `t  t  >
 *                               ?  )  >
 *                            _,,:cc,,,;.
 *                      ,,c$$$$$$;$$$$$$$$$$$c,
 *                  ,c$$$$$$$$$$$$$$$$$$$$$$$$$$h,
 *               ,c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$c,
 *            ,d$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$,
 *          ,$$$$$",d$$$$$$$$$$$$$$$$$$$$$$$$$hc`?$$$$$c
 *         ,$$$$$u$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$h,?$$$$L
 *         J$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$b$$$$$>
 *         $$$$$$$P$$?$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$>
 *         $$$$$$$$`$b`$$$$$$$$$$$$$$$$$$$$$"$"d$$$$$$$$$
 *         `$$$$$$c ",nn,?$$$$$$$$$$$$$$$",,`r",d$$$$$$$'
 *          `$$$$$$$ MMMMM`$$$$$$$$$$$$F,MMMb`$$$$$$$$$'
 *           `$$$$$${MMTTM.?$$$$$$$$$$$,MMMMM $$$$$$$$'
 *            `$$$$$,M";;; `$$$$$$$$$$'M",,`",$$$$$$$'
 *              ?$$$$,{( `) $$$$$$$$$$ ,('`) J$$$$$P'
 *               ?$$$$,(  ) $$$$$$$$$$ ('  ),$$$$$F
 *                `$$$$.`-' $$$$$$$$$$,`--',$$$$$'
 *                  $$$$$hh$$$$$????$$$hc$$$$$$$'  I tawt I taw an ugly F77 code !
 *                 d$$$$$$$$$ `======' $$$$$$$$$;
 *                 $$$$$$$$$$$$c,,,,c$$$$$$$$$$$'
 *                  "?$$$$P"" "$$$$$$?????$$??"
 *                             $$$$$
 *                            4$$$$$c
 *                          ,$$$$$$$"c
 *                         z${$$$$$$$`$,
 *                        z${$$$$$$$$$`$c
 *                       {$>$$$$$$$$$$;?$
 *                       `${$$$$$$$$$$$:$
 *                        ?L$$$$$$$$$$$:$
 *                         ?$$$$$$$$$$$d'
 *                         `$$$$$$$$$$F
 *                          `?$c`??3$F
 *                           CCC   CCC
 *     ,;{CC>C>;,,.          `CC   CCC              .,,,,,.
 *    CCC{CCCCCCCCCCCCCCCCCCCCCC   CC>CCCCCCCCCCCCCCCCCCCCCC>,
 *   {>,CCCCCCCCCCCCCCCCCCCCCCCC, `CCCCCCCCCCCCCCCCCCCCCC(`{CC
 *     CCCCCCCCCCCCCCCCCCCCCCC>'   `'CCCCCCCCCCCCCCCCCCCCCCC`C
 *      `''{{CCCCCCCC>>''               `'{{CCCCCCCCCCCCCCC>
 */

//==============================================================================
//==============================================================================
//==============================================================================

/** The constructor.
 *
 *  Open a rhl.dat file, build the corresponding DataTree.
 */

IOberger::IOberger(void)
{
  DBG_ENTER;

  DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** The destructor.
 */

IOberger::~IOberger(void)
{
  DBG_ENTER;

  if (buffer != NULL)
  {
    free(buffer);
  }

  DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Read the basis parameters and create a local Basis instance.
 */

void IOberger::readBasis(void)
{
  DBG_ENTER;

  alpha         = -999.0;
  beta          = -999.0;
  d_0           = -999.0;
  nOscil        = -999;
  n_zMaxImposed = -999;
  gq            = -999.0;
  jdx           = -999;

  noeRead("alf", "D", 1, &alpha);
  noeRead("bet", "D", 1, &beta);
  noeRead("d1", "D", 1, &d_0);
  noeRead("nn", "I", 1, &nOscil);
  noeRead("maxnz", "I", 1, &n_zMaxImposed);
  noeRead("gq", "D", 1, &gq);
  noeRead("jdx", "I", 1, &jdx);
  char cdummy = 0;
  v2ct = false;
  noeRead("v2ct", "L", 1, &cdummy);

  if (cdummy != 0) v2ct = true;

  if (  (alpha         < 0.0)
        || (beta          < 0.0)
        || (d_0           < 0.0)
        || (nOscil        < 0  )
        || (n_zMaxImposed < 0  )
        || (gq            < 0.0)
        || (jdx           < 0  ))
  {
    ERROR("cannot construct basis");
  }

  double b_r = 1 / sqrt(beta);
  double b_z = 1 / sqrt(alpha);

  if (v2ct)
  {
    basis = Basis(d_0, b_r, b_z, nOscil, n_zMaxImposed, gq);
  }
  else
  {
    basis = Basis(0, b_r, b_z, nOscil, n_zMaxImposed, gq);
  }

  mxMax = basis.mMax;

  DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Read misc. quantities.
 */

void IOberger::readMisc(void)
{
  DBG_ENTER;

  double toto[27];
  double erearr[9];
  double eneb[3 * 4 * 8];
  double ecdm2[3];
  double ecoulx;
  char text[512];
  char type[512];
  double v1, v2;
  INT i1;
  double xlmb[4];


  try
  {
    noeRead("xlmb", "D", 4, &xlmb);
  }
  catch (...)
  {
  }

  chemPot = arma::vec(2);
  chemPot(0) = xlmb[2];
  chemPot(1) = xlmb[3];

  try
  {
    noeRead("xbz", "D", 1, &zNumber);
    noeRead("xbn", "D", 1, &nNumber);
    noeRead("wtab(iwtab,09)", "D", 1, &convergence);
  }
  catch (...)
  {
  }

  niterx = -1;
  nitert = -1;

  try
  {
    noeRead("iterx", "I", 1, &niterx);
    noeRead("itert", "I", 1, &nitert);
  }
  catch (...)
  {
  }

  try
  {
    noeRead("finalcvg", "D", 1, &convergence);
  }
  catch (...)
  {
  }

  noeRead("qsciss", "D", 9, &toto);
  qsciss = toto[8];

  if (noeRead("nomjob", "S", 1, text)) jobname = text;

  if (noeRead("nomro1", "S", 1, text)) parent = text;

  nbcon = 0;
  noeRead("nbcon", "I", 1, &nbcon);

  if (nbcon <= 0)
  {
    std::stringstream ss;
    ss << "dimension 'nbcon' not found in IOberger::noeRead()";
    Tools::warning(ss.str());
  }
  else
  {
    constraint_val = arma::zeros(nbcon);
    constraint_pond = arma::zeros(nbcon);
    constraint_iter1c = IVEC(nbcon, arma::fill::zeros);
    constraint_lambda = arma::zeros(nbcon);
    noeRead("xlamci", "D", nbcon, &toto);

    for (INT i = 0; i < nbcon; i++)
    {
      sprintf(text, "constraint %02d", i + 1);
      noeRead(text, "SDDI", 1, type, &v1, &v2, &i1);
      std::string stype = type;

      // conversion from "cdm(imp)" to "Q10"
      if (stype == "cdm(imp)") stype = "Q10";

      constraint_type.push_back(stype);
      constraint_val(i) = v1;
      constraint_pond(i) = v2;
      constraint_iter1c(i) = i1;
      constraint_lambda(i) = toto[i];
    }
  }

  noeRead("ehbtot", "D", 1, &eneTot);
  /*
    the berger2ct's values come from the eij[12 * idene + 3 * idct + iso] array where:
      - idene is the energy id :
        . 0 -> cin+pot -> without pairing but with cdm2
        . 1 -> cinet.
        . 2 -> poten. -> without pairing and without cdm2
        . 3 -> coul.
        . 4 -> cent.dir
        . 5 -> cent.ech
        . 6 -> spin-orb
        . 7 -> densite
      - idct is the id of center :
        . 0 -> 1
        . 1 -> 12
        . 2 -> 2
        . 3 -> t
      - iso is the isospin :
        . 0 -> neutron
        . 1 -> proton
        . 2 -> total
    this array is saved in zfinit.f from eij(iso, idct, idene).
  */
  noeRead("eij", "D", 3 * 4 * 8, eneb);
  noeRead("ecoulx", "D", 1, &ecoulx);
  noeRead("ecdm2", "D", 3, ecdm2);
  noeRead("xebcs", "D", 9, toto);
  noeRead("e3c", "D", 9, erearr);
  energiesNeut["Kinetic"          ] = eneb[12 * 1 + 9 ];
  energiesNeut["Coulomb_Direct"   ] = eneb[12 * 3 + 9 ];
  energiesNeut["Central_Direct"   ] = eneb[12 * 4 + 9 ];
  energiesNeut["Central_Exchange" ] = eneb[12 * 5 + 9 ];
  energiesNeut["Spin-orbit"       ] = eneb[12 * 6 + 9 ];
  energiesNeut["Density"          ] = eneb[12 * 7 + 9 ];
  energiesNeut["2-Body_COM_Cor."  ] = ecdm2[0];
  energiesNeut["Coulomb_Exchange" ] = 0.;
  energiesNeut["Pairing_Central"  ] = toto[6];
  energiesNeut["Rearrangement"    ] = erearr[6] * -1.0;
  energiesProt["Kinetic"          ] = eneb[12 * 1 + 10];
  energiesProt["Coulomb_Direct"   ] = eneb[12 * 3 + 10];
  energiesProt["Central_Direct"   ] = eneb[12 * 4 + 10];
  energiesProt["Central_Exchange" ] = eneb[12 * 5 + 10];
  energiesProt["Spin-orbit"       ] = eneb[12 * 6 + 10];
  energiesProt["Density"          ] = eneb[12 * 7 + 10];
  energiesProt["2-Body_COM_Cor."  ] = ecdm2[1];
  energiesProt["Coulomb_Exchange" ] = ecoulx;
  energiesProt["Pairing_Central"  ] = toto[7];
  energiesProt["Rearrangement"    ] = erearr[7] * -1.0;
  // These values can be used by the D12S interaction (limit D1S with D2 form)
  energiesNeut["Density_D2"       ] = energiesNeut["Density"];
  energiesProt["Density_D2"       ] = energiesProt["Density"];
  energiesNeut["Rearrangement_D2" ] = energiesNeut["Rearrangement"];
  energiesProt["Rearrangement_D2" ] = energiesProt["Rearrangement"];
  //  energiesNeut["Pairing_D2"   ] = energiesNeut["Pairing_Central"];
  //  energiesProt["Pairing_D2"   ] = energiesProt["Pairing_Central"];
  corrections = arma::zeros(6);
  noeRead("ztroa", "D", 1, &v1);
  corrections(0) = v1;
  noeRead("ztrog", "D", 1, &v1);
  corrections(1) = v1;
  noeRead("ztroqa", "D", 1, &v1);
  corrections(2) = v1;
  noeRead("ztroqg", "D", 1, &v1);
  corrections(3) = v1;
  noeRead("ztroja", "D", 1, &v1);
  corrections(4) = v1;
  noeRead("ztrojg", "D", 1, &v1);
  corrections(5) = v1;
  noeRead("masses0", "D", 15, toto);
  inertia0 = arma::zeros(15);

  for (INT i = 0; i < 15; i++)
    inertia0(i) = toto[i];

  noeRead("masses23", "D", 12, toto);
  inertia23 = arma::zeros(12);

  for (INT i = 0; i < 12; i++)
    inertia23(i) = toto[i];

  noeRead("masses24", "D", 12, toto);
  inertia24 = arma::zeros(12);

  for (INT i = 0; i < 12; i++)
    inertia24(i) = toto[i];

  noeRead("masses34", "D", 12, toto);
  inertia34 = arma::zeros(12);

  for (INT i = 0; i < 12; i++)
    inertia34(i) = toto[i];

  noeRead("masses234", "D", 27, toto);
  inertia234 = arma::zeros(27);

  for (INT i = 0; i < 27; i++)
    inertia234(i) = toto[i];





  INT ke = -1;
  noeRead("ke", "I", 1, &ke);

  if (ke <= 0)
  {
    std::stringstream ss;
    ss << "dimension 'ke' not found in IOberger::noeRead()";
    Tools::warning(ss.str());
  }
  else
  {
    Qnumbers QNid("id", ke);
    matDn = arma::zeros(ke, basis.HOqn.nb);
    matDp = arma::zeros(ke, basis.HOqn.nb);
    INT nord[ke];
    INT nordw[ke];
    double xoc[ke];
    double xocw[ke];
    double xenr[ke];
    double xenrw[ke];
    INT kenr[basis.mMax + 1];
    noeRead("nord", "I", ke, nord);
    noeRead("nordw", "I", ke, nordw);
    noeRead("xoc", "D", ke, xoc);
    noeRead("xocw", "D", ke, xocw);
    noeRead("xenr", "D", ke, xenr);
    noeRead("xenrw", "D", ke, xenrw);
    noeRead("kenr", "I", basis.mMax + 1, kenr);
    vecVn = arma::zeros(ke);
    vecVp = arma::zeros(ke);
    eneQPn = arma::zeros(ke);
    eneQPp = arma::zeros(ke);
    vecOmegan = IVEC(ke, arma::fill::zeros);
    vecOmegap = IVEC(ke, arma::fill::zeros);
    double toto2[jdx * jdx];
    INT istate = 0;

    if (mxMax == -1) ERROR("mxMax == -1");

    for (INT mo = 0; mo < mxMax; mo++)
    {
      arma::mat xav  = arma::zeros(jdx, jdx);
      arma::mat xavw = arma::zeros(jdx, jdx);
      char label[128];

      sprintf(label, "xav(%2.2d)", mo + 1);
      noeRead(label, "D", jdx * jdx, toto2);
      INT id = 0;

      for (INT i = 0; i < jdx; i++)
      {
        for (INT j = 0; j < jdx; j++)
        {
          xav(i, j) = toto2[id];
          id++;
        }
      }

      sprintf(label, "xavw(%2.2d)", mo + 1);
      noeRead(label, "D", jdx * jdx, toto2);
      id = 0;

      for (INT i = 0; i < jdx; i++)
      {
        for (INT j = 0; j < jdx; j++)
        {
          xavw(i, j) = toto2[id];
          id++;
        }
      }

      INT i = 0;

      for (INT a = kenr[mo]; a < kenr[mo + 1]; a += 1)
      {
        INT j = 0;
        INT mbstop = mo + 2;

        if (mo >= basis.mMax - 1) mbstop = basis.mMax;

        INT rankn = -1;

        for (INT rank = 0; rank < kenr[basis.mMax]; rank++)
        {
          if (nord[rank] == istate + 1)
          {
            rankn = rank;
            break;
          }
        }

        INT rankp = -1;

        for (INT rank = 0; rank < kenr[basis.mMax]; rank++)
        {
          if (nordw[rank] == istate + 1)
          {
            rankp = rank;
            break;
          }
        }

        if ((rankn != -1) && (rankp != -1))
        {
          vecVn(rankn) = xoc[istate];
          vecOmegan(rankn) = mo;
          eneQPn(istate) = xenr[istate];
          vecVp(rankp) = xocw[istate];
          eneQPp(istate) = xenrw[istate];
          vecOmegap(rankp) = mo;

          for (INT mb = mo; mb < mbstop; mb++)
          {
            for (INT nb = 0; nb < basis.nMax(mb); nb++)
            {
              for (INT db = 0; db < basis.dMax; db++)
              {
                for (INT n_zb = 0; n_zb < basis.n_zMax(mb, nb); n_zb++)
                {
                  INT sb = mb - mo;
                  INT b = basis.HOqn.find({mb, nb, n_zb, db, sb});
                  matDn(rankn, b) = xav(i, j);
                  matDp(rankp, b) = xavw(i, j);
                  j++;

                  ASSERT(j <= jdx, "j > jdx");
                }
              }
            }
          }
        }

        i++;
        istate++;
      }
    }
  }

  DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Read rho and kappa(NEUTRON) matrices.
 */

void IOberger::readRhoKappa(void)
{
  DBG_ENTER;
  INT nro0 = 0;
  INT nro1 = 0;
  INT nro2 = 0;
  noeRead("nro0", "I", 1, &nro0);
  noeRead("nro1", "I", 1, &nro1);
  noeRead("nro2", "I", 1, &nro2);
  arma::vec ro0;
  arma::vec ro1;
  arma::vec ro2;
  arma::vec ro0w;
  arma::vec ro1w;
  arma::vec ro2w;
  arma::vec xcp0;
  arma::vec xcp1;
  arma::vec xcp2;
  arma::vec xcp0w;
  arma::vec xcp1w;
  arma::vec xcp2w;

  if (nro0 != 0)
  {
    double toto1[nro0];

    noeRead("ro0", "H", nro0, toto1);
    ro0 = arma::zeros(nro0);

    for (INT i = 0; i < nro0; i++) ro0(i) = toto1[i];

    noeRead("ro0w", "H", nro0, toto1);
    ro0w = arma::zeros(nro0);

    for (INT i = 0; i < nro0; i++) ro0w(i) = toto1[i];

    noeRead("xcp0", "H", nro0, toto1);
    xcp0 = arma::zeros(nro0);

    for (INT i = 0; i < nro0; i++) xcp0(i) = toto1[i];

    noeRead("xcp0w", "H", nro0, toto1);
    xcp0w = arma::zeros(nro0);

    for (INT i = 0; i < nro0; i++) xcp0w(i) = toto1[i];

  }

  if (nro1 != 0)
  {
    double toto2[nro1];

    noeRead("ro1", "H", nro1, toto2);
    ro1 = arma::zeros(nro1);

    for (INT i = 0; i < nro1; i++) ro1(i) = toto2[i];

    noeRead("ro1w", "H", nro1, toto2);
    ro1w = arma::zeros(nro1);

    for (INT i = 0; i < nro1; i++) ro1w(i) = toto2[i];

    noeRead("xcp1", "H", nro1, toto2);
    xcp1 = arma::zeros(nro1);

    for (INT i = 0; i < nro1; i++) xcp1(i) = toto2[i];

    noeRead("xcp1w", "H", nro1, toto2);
    xcp1w = arma::zeros(nro1);

    for (INT i = 0; i < nro1; i++) xcp1w(i) = toto2[i];

  }

  if (nro2 != 0)
  {
    double toto3[nro2];

    noeRead("ro2", "H", nro2, toto3);
    ro2 = arma::zeros(nro2);

    for (INT i = 0; i < nro2; i++) ro2(i) = toto3[i];

    noeRead("ro2w", "H", nro2, toto3);
    ro2w = arma::zeros(nro2);

    for (INT i = 0; i < nro2; i++) ro2w(i) = toto3[i];

    noeRead("xcp2", "H", nro2, toto3);
    xcp2 = arma::zeros(nro2);

    for (INT i = 0; i < nro2; i++) xcp2(i) = toto3[i];

    noeRead("xcp2w", "H", nro2, toto3);
    xcp2w = arma::zeros(nro2);

    for (INT i = 0; i < nro2; i++) xcp2w(i) = toto3[i];

  }

  rho(NEUTRON) = arma::zeros(basis.HOqn.nb, basis.HOqn.nb);
  rho(PROTON ) = arma::zeros(basis.HOqn.nb, basis.HOqn.nb);
  kappa(NEUTRON) = arma::zeros(basis.HOqn.nb, basis.HOqn.nb);
  kappa(PROTON ) = arma::zeros(basis.HOqn.nb, basis.HOqn.nb);
  IVEC nzMax = basis.n_zMaxBerger;

  INT md = 0;
  INT nd = 0;
  INT nzd = 0;
  INT dd = 0;
  INT mb = 0;
  INT nb = 0;
  INT nzb = 0;
  INT db = 0;
  INT m = 0;
  INT mo = 0;
  INT iswap = 1;
  INT i0 = 0;
  INT i2 = 0;
  INT i3 = 0;

  IMAT mati0 = arma::ones<IMAT >(basis.HOqn.nb, basis.HOqn.nb) * 0;
  IMAT mati1 = arma::ones<IMAT >(basis.HOqn.nb, basis.HOqn.nb) * 0;
  IMAT mati2 = arma::ones<IMAT >(basis.HOqn.nb, basis.HOqn.nb) * 0;


  for (mo = 1; mo <= mxMax; mo++)
  {
    INT nx = (mxMax - mo) / 2 + 1;
    INT nnx = (mxMax - mo - 1) / 2 + 1;
    m = mo;
    md = mo;
    mb = mo;
reboucle:

    for (nb = 1; nb <= nx; nb++)
    {
      for (nd = 1; nd <= nb; nd++)
      {
        INT nzbx1 = nzMax(nb * 2 + m - 3);
        INT nzbx2 = nzMax(nb * 2 + m - 3);
        INT nzdx1 = nzMax(nd * 2 + m - 3);
        INT nzdx2 = nzMax(nd * 2 + m - 3);

        if (!v2ct) nzbx1 = 0;

        if (!v2ct) nzdx1 = 0;

        if (nzbx1 != 0)
        {
          for (nzb = 1; nzb <= nzbx1; nzb++)
          {
            INT nzdx = 0;

            if (nb == nd) nzdx = nzb;
            else nzdx = nzdx1;

            for (nzd = 1; nzd <= nzdx; nzd++)
            {
              if (nzd <= nzdx1)
              {
                i3++;

                if (v2ct)
                {
                  dd = 0;
                  db = 0;
                }

                double vrhon = ro1[i3 - 1];
                double vrhop = ro1w[i3 - 1];
                double vkappan = xcp1[i3 - 1];
                double vkappap = xcp1w[i3 - 1];
                INT sd = md - mo;
                INT sb = mb - mo;
                INT idd = basis.HOqn.find({md - 1, nd - 1, nzd - 1, dd, sd});
                INT idb = basis.HOqn.find({mb - 1, nb - 1, nzb - 1, db, sb});

                mati1(idd, idb) = i3;
                mati1(idb, idd) = i3;

                //                INFO("I1: %04d (%2d %2d %2d %2d %2d) (%2d %2d %2d %2d %2d)", i3, md - 1, nd - 1, nzd - 1, dd, sd, mb - 1, nb - 1, nzb - 1, db, sb);

                rho(NEUTRON)(idd, idb) = vrhon;
                rho(PROTON )(idd, idb) = vrhop;
                rho(NEUTRON)(idb, idd) = vrhon;
                rho(PROTON )(idb, idd) = vrhop;
                kappa(NEUTRON)(idd, idb) = vkappan;
                kappa(PROTON )(idd, idb) = vkappap;
                kappa(NEUTRON)(idb, idd) = vkappan;
                kappa(PROTON )(idb, idd) = vkappap;
              }

              if (nb != nd)
              {
                i2++;

                if (v2ct)
                {
                  dd = 1;
                  db = 0;
                }

                double vrhon = ro2[i2 - 1];
                double vrhop = ro2w[i2 - 1];
                double vkappan = xcp2[i2 - 1];
                double vkappap = xcp2w[i2 - 1];
                INT sd = md - mo;
                INT sb = mb - mo;
                INT idd = basis.HOqn.find({md - 1, nd - 1, nzd - 1, dd, sd});
                INT idb = basis.HOqn.find({mb - 1, nb - 1, nzb - 1, db, sb});

                mati2(idd, idb) = i2;
                mati2(idb, idd) = i2;

                //                INFO("I2: %04d (%2d %2d %2d %2d %2d) (%2d %2d %2d %2d %2d)", i2, md - 1, nd - 1, nzd - 1, dd, sd, mb - 1, nb - 1, nzb - 1, db, sb);

                rho(NEUTRON)(idd, idb) = vrhon;
                rho(PROTON )(idd, idb) = vrhop;
                rho(NEUTRON)(idb, idd) = vrhon;
                rho(PROTON )(idb, idd) = vrhop;
                kappa(NEUTRON)(idd, idb) = vkappan;
                kappa(PROTON )(idd, idb) = vkappap;
                kappa(NEUTRON)(idb, idd) = vkappan;
                kappa(PROTON )(idb, idd) = vkappap;
              }
            }
          }
        }

        for (nzb = 1; nzb <= nzbx2; nzb++)
        {
          for (nzd = 1; nzd <= nzdx2; nzd++)
          {
            if (nzd <= nzdx1)
            {
              i2++;

              if (v2ct)
              {
                dd = 0;
                db = 1;
              }

              double vrhon = ro2[i2 - 1];
              double vrhop = ro2w[i2 - 1];
              double vkappan = xcp2[i2 - 1];
              double vkappap = xcp2w[i2 - 1];
              INT sd = md - mo;
              INT sb = mb - mo;
              INT idd = basis.HOqn.find({md - 1, nd - 1, nzd - 1, dd, sd});
              INT idb = basis.HOqn.find({mb - 1, nb - 1, nzb - 1, db, sb});

              mati2(idd, idb) = i2;
              mati2(idb, idd) = i2;

              //              INFO("I2: %04d (%2d %2d %2d %2d %2d) (%2d %2d %2d %2d %2d)", i2, md - 1, nd - 1, nzd - 1, dd, sd, mb - 1, nb - 1, nzb - 1, db, sb);

              rho(NEUTRON)(idd, idb) = vrhon;
              rho(PROTON )(idd, idb) = vrhop;
              rho(NEUTRON)(idb, idd) = vrhon;
              rho(PROTON )(idb, idd) = vrhop;
              kappa(NEUTRON)(idd, idb) = vkappan;
              kappa(PROTON )(idd, idb) = vkappap;
              kappa(NEUTRON)(idb, idd) = vkappan;
              kappa(PROTON )(idb, idd) = vkappap;
            }

            if ((nb != nd) || (nzd <= nzb))
            {
              i0++;

              if (v2ct)
              {
                dd = 1;
                db = 1;
              }

              double vrhon = ro0[i0 - 1];
              double vrhop = ro0w[i0 - 1];
              double vkappan = xcp0[i0 - 1];
              double vkappap = xcp0w[i0 - 1];
              INT sd = md - mo;
              INT sb = mb - mo;
              INT idd = basis.HOqn.find({md - 1, nd - 1, nzd - 1, dd, sd});
              INT idb = basis.HOqn.find({mb - 1, nb - 1, nzb - 1, db, sb});

              mati0(idd, idb) = i0;
              mati0(idb, idd) = i0;

              //              INFO("I0: %04d (%2d %2d %2d %2d %2d) (%2d %2d %2d %2d %2d)", i0, md - 1, nd - 1, nzd - 1, dd, sd, mb - 1, nb - 1, nzb - 1, db, sb);

              rho(NEUTRON)(idd, idb) = vrhon;
              rho(PROTON )(idd, idb) = vrhop;
              rho(NEUTRON)(idb, idd) = vrhon;
              rho(PROTON )(idb, idd) = vrhop;
              kappa(NEUTRON)(idd, idb) = vkappan;
              kappa(PROTON )(idd, idb) = vkappap;
              kappa(NEUTRON)(idb, idd) = vkappan;
              kappa(PROTON )(idb, idd) = vkappap;
            }
          }
        }
      }
    }

    //==============================================================================
    iswap = 1 - iswap;

    if (iswap == 1) continue;

    md = mo;
    mb = mo + 1;

    if (mo == mxMax) continue;

    for (nb = 1; nb <= nnx; nb++)
    {
      for (nd = 1; nd <= nx; nd++)
      {
        INT nzbx1 = nzMax(nb * 2 + m + 1 - 3);
        INT nzbx2 = nzMax(nb * 2 + m + 1 - 3);
        INT nzdx1 = nzMax(nd * 2 + m - 3);
        INT nzdx2 = nzMax(nd * 2 + m - 3);

        if (!v2ct) nzbx1 = 0;

        if (!v2ct) nzdx1 = 0;

        if (nzbx1 != 0)
        {
          for (nzb = 1; nzb <= nzbx1; nzb++)
          {
            for (nzd = 1; nzd <= nzdx2; nzd++)
            {
              if (nzd <= nzdx1)
              {
                i3++;

                if (v2ct)
                {
                  dd = 0;
                  db = 0;
                }

                double vrhon = ro1[i3 - 1];
                double vrhop = ro1w[i3 - 1];
                double vkappan = xcp1[i3 - 1];
                double vkappap = xcp1w[i3 - 1];
                INT sd = md - mo;
                INT sb = mb - mo;
                INT idd = basis.HOqn.find({md - 1, nd - 1, nzd - 1, dd, sd});
                INT idb = basis.HOqn.find({mb - 1, nb - 1, nzb - 1, db, sb});

                mati1(idd, idb) = i3;
                mati1(idb, idd) = i3;

                //                INFO("I1: %04d (%2d %2d %2d %2d %2d) (%2d %2d %2d %2d %2d)", i3, md - 1, nd - 1, nzd - 1, dd, sd, mb - 1, nb - 1, nzb - 1, db, sb);

                rho(NEUTRON)(idd, idb) = vrhon;
                rho(PROTON )(idd, idb) = vrhop;
                rho(NEUTRON)(idb, idd) = vrhon;
                rho(PROTON )(idb, idd) = vrhop;
                kappa(NEUTRON)(idd, idb) = vkappan;
                kappa(PROTON )(idd, idb) = vkappap;
                kappa(NEUTRON)(idb, idd) = vkappan;
                kappa(PROTON )(idb, idd) = vkappap;
              }

              i2++;

              if (v2ct)
              {
                dd = 1;
                db = 0;
              }

              double vrhon = ro2[i2 - 1];
              double vrhop = ro2w[i2 - 1];
              double vkappan = xcp2[i2 - 1];
              double vkappap = xcp2w[i2 - 1];
              INT sd = md - mo;
              INT sb = mb - mo;
              INT idd = basis.HOqn.find({md - 1, nd - 1, nzd - 1, dd, sd});
              INT idb = basis.HOqn.find({mb - 1, nb - 1, nzb - 1, db, sb});

              mati2(idd, idb) = i2;
              mati2(idb, idd) = i2;

              //              INFO("I2: %04d (%2d %2d %2d %2d %2d) (%2d %2d %2d %2d %2d)", i2, md - 1, nd - 1, nzd - 1, dd, sd, mb - 1, nb - 1, nzb - 1, db, sb);

              rho(NEUTRON)(idd, idb) = vrhon;
              rho(PROTON )(idd, idb) = vrhop;
              rho(NEUTRON)(idb, idd) = vrhon;
              rho(PROTON )(idb, idd) = vrhop;
              kappa(NEUTRON)(idd, idb) = vkappan;
              kappa(PROTON )(idd, idb) = vkappap;
              kappa(NEUTRON)(idb, idd) = vkappan;
              kappa(PROTON )(idb, idd) = vkappap;
            }
          }
        }

        for (nzb = 1; nzb <= nzbx2; nzb++)
        {
          for (nzd = 1; nzd <= nzdx2; nzd++)
          {
            if (nzd <= nzdx1)
            {
              i2++;

              if (v2ct)
              {
                dd = 0;
                db = 1;
              }

              double vrhon = ro2[i2 - 1];
              double vrhop = ro2w[i2 - 1];
              double vkappan = xcp2[i2 - 1];
              double vkappap = xcp2w[i2 - 1];
              INT sd = md - mo;
              INT sb = mb - mo;
              INT idd = basis.HOqn.find({md - 1, nd - 1, nzd - 1, dd, sd});
              INT idb = basis.HOqn.find({mb - 1, nb - 1, nzb - 1, db, sb});

              mati2(idd, idb) = i2;
              mati2(idb, idd) = i2;

              //              INFO("I2: %04d (%2d %2d %2d %2d %2d) (%2d %2d %2d %2d %2d)", i2, md - 1, nd - 1, nzd - 1, dd, sd, mb - 1, nb - 1, nzb - 1, db, sb);

              rho(NEUTRON)(idd, idb) = vrhon;
              rho(PROTON )(idd, idb) = vrhop;
              rho(NEUTRON)(idb, idd) = vrhon;
              rho(PROTON )(idb, idd) = vrhop;
              kappa(NEUTRON)(idd, idb) = vkappan;
              kappa(PROTON )(idd, idb) = vkappap;
              kappa(NEUTRON)(idb, idd) = vkappan;
              kappa(PROTON )(idb, idd) = vkappap;
            }

            i0++;

            if (v2ct)
            {
              dd = 1;
              db = 1;
            }

            double vrhon = ro0[i0 - 1];
            double vrhop = ro0w[i0 - 1];
            double vkappan = xcp0[i0 - 1];
            double vkappap = xcp0w[i0 - 1];
            INT sd = md - mo;
            INT sb = mb - mo;
            INT idd = basis.HOqn.find({md - 1, nd - 1, nzd - 1, dd, sd});
            INT idb = basis.HOqn.find({mb - 1, nb - 1, nzb - 1, db, sb});

            mati0(idd, idb) = i0;
            mati0(idb, idd) = i0;

            //            INFO("I0: %04d (%2d %2d %2d %2d %2d) (%2d %2d %2d %2d %2d)", i0, md - 1, nd - 1, nzd - 1, dd, sd, mb - 1, nb - 1, nzb - 1, db, sb);

            rho(NEUTRON)(idd, idb) = vrhon;
            rho(PROTON )(idd, idb) = vrhop;
            rho(NEUTRON)(idb, idd) = vrhon;
            rho(PROTON )(idb, idd) = vrhop;
            kappa(NEUTRON)(idd, idb) = vkappan;
            kappa(PROTON )(idd, idb) = vkappap;
            kappa(NEUTRON)(idb, idd) = vkappan;
            kappa(PROTON )(idb, idd) = vkappap;
          }
        }
      }
    }

    md = mo + 1;
    mb = mo + 1;
    m = mo + 1;
    nx = nnx;
    goto reboucle;
  }

  // check if Omega-symmetry is present in rho and kappa(NEUTRON) matrices.
  basis.HOqn.checkOmegaSym(rho(NEUTRON), "rho(NEUTRON)");
  basis.HOqn.checkOmegaSym(rho(PROTON ), "rho(PROTON )");
  basis.HOqn.checkOmegaSym(kappa(NEUTRON), "kappa(NEUTRON)");
  basis.HOqn.checkOmegaSym(kappa(PROTON ), "kappa(PROTON )");

  DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Hexadecimal to double format converter.
 *
 *  \param doubleVal The double value.
 *  \param size The size of the hexadecimal value.
 *  \param hexaVal The hexadecimal value.
 */

void IOberger::hexaToDouble(void *doubleVal, INT size, const char *hexaVal)
{
  INT i;
  unsigned char *ptr;
  char toto[3];
  INT c;
  ptr = static_cast<unsigned char *>(doubleVal);

  for (i = 0; i < size; i++)
  {
    toto[0] = hexaVal[i * 2];
    toto[1] = hexaVal[i * 2 + 1];
    toto[2] = 0;
    sscanf(toto, "%X", &c);
    ptr[size - 1 - i] = (unsigned char)c;
  }
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Build an index of the rhl.dat file.
 *
 *  \return The number of labels.
 */

INT IOberger::noeIndex(void)
{
  DBG_ENTER;

  std::size_t found = 0;

  while (true)
  {
    found = content.find('\n', found + 1);

    if (found == std::string::npos)
    {
      break;
    }

    if (found + 1 < content.size())
    {
      if (content[found + 1] != '#')
      {
        std::size_t found2 = content.find(':', found + 1);

        if (found2 != std::string::npos)
        {
          labelPos[content.substr(found + 1, found2 - found - 2)] = int(found2 + 1);
          nbLabels++;
        }
      }
    }
  }

  DBG_RETURN(nbLabels);
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Search the index of the rhl.dat file for a given label.
 *
 * \param label The label.
 * \return The index value.
 */

INT IOberger::noeSearch(const std::string label) const
{
  if (labelPos.count(label) == 0)
  {
    return -1;
  }

  return labelPos.at(label);
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Interpret the content of an rhl.dat.gz file.
 *
 * \return The number of labels (-1 in case of error).
 */

INT IOberger::noeOpen(void)
{
  DBG_ENTER;

  std::stringstream stream;
  stream << content;
  char word[11];
  INT i;

  for (i = 0; i < 10; i++)
  {
    word[i] = 0;

    if (!stream.eof()) word[i] = (char)(stream.get());
  }

  word[10] = 0;

  if (strcmp(word, "#noe000001") != 0) ERROR("wrong format");

  DBG_RETURN(noeIndex());
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Read a label from a rhl.dat file.
 */

bool IOberger::noeRead(const std::string label, const std::string type, INT n, ...)
{
  DBG_ENTER;

  INT ii;
  INT i, j;
  void **v;
  char ctoto;
  char toto[512];
  UINT nsub = type.length();
  va_list Numbers;
  va_start(Numbers, n);
  v = new void *[nsub];

  ASSERT(v != NULL, "can not allocate memory");

  for (INT i = 0; i < nsub; i++)
  {
    v[i] = va_arg(Numbers, void *);
  }

  ii = noeSearch(label);

  if (ii == -1)
  {
    std::stringstream ss;
    ss << "label '" << label << "' not found in IOberger::noeRead()";
    Tools::warning(ss.str());
    DBG_RETURN(false);
  }

  std::size_t istart = ii + 1;

  for (i = 0; i < n; i++)
    for (j = 0; j < nsub; j++)
    {
      if (type[j] == 'S')
      {
        std::size_t found1 = content.find("'", istart + 1);
        std::string tutu2 = content.substr(istart + 1, found1 - istart - 1);
        istart = found1 + 2;
        strcpy((static_cast<char *>(v[j])), tutu2.c_str());
      }
      else
      {
        std::size_t found1 = content.find(' ', istart + 1);
        std::string tutu2 = content.substr(istart, found1 - istart);
        istart = found1 + 1;
        strcpy(toto, tutu2.c_str());

        if (type[j] == 'I')
        {
          sscanf(toto, "%d ", &(static_cast<INT *>((v[j])))[i]);
        }

#ifdef ALWAYS_HEXA

        if ((type[j] == 'H') || (type[j] == 'D'))
        {
          hexaToDouble(&(static_cast<double *>((v[j])))[i], sizeof(double), toto);
        }

#else

        if (type[j] == 'H')
        {
          h2f(&((double *)(v[j]))[i], sizeof(double), toto);
        }

        if (type[j] == 'D') sscanf(toto, "%lf ", &((double *)(v[j]))[i]);

#endif

        if (type[j] == 'L')
        {
          sscanf(toto, "%c ", &ctoto);

          if (ctoto == 'T') ((static_cast<char *>(v[j]))[i]) = -1;
          else            ((static_cast<char *>(v[j]))[i]) = 0;
        }

        if (type[j] == 'C') sscanf(toto, "'%c' ", &(static_cast<char *>(v[j]))[i]);
      }

      /*
            if (type[j] == 'S')
            {
              INT k = 0;
              stream.seekg(idebut);
              ctoto = stream.get();
              ctoto = stream.get();

              while (ctoto != '\'')
              {
                (static_cast<char *>(v[j]))[k++] = ctoto;
                ctoto = stream.get();
              }

              ctoto = stream.get();
              (static_cast<char *>(v[j]))[k] = 0;
            }
      */
    }

  delete[] v;

  DBG_RETURN(true);
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Return a DataTree instance from BERGER2CT formatted content.
 *
 * \param _content The content of a rhl.dat file.
 */

DataTree IOberger::fromContent(const std::string &_content)
{
  DBG_ENTER;
  DataTree result;

  if (!checkFileType(_content))
  {
    DBG_RETURN(result);
  }

  content = _content;
  nbLabels = 0;
  noeOpen();
  noeIndex();
  readBasis();
  readMisc();
  readRhoKappa();

  result.merge(basis.getDataTree("initial/")); // keep a copy of the initial basis parameters in case of a basis conversion
  result.merge(basis.getDataTree(""));
  result.set("HFB/rho", rho);
  result.set("HFB/kappa", kappa);
  result.set("HFB/N", (INT)(nNumber));
  result.set("HFB/Z", (INT)(zNumber));
  result.set("HFB/HOtoHFNeut", matDn);
  result.set("HFB/HOtoHFProt", matDp);
  result.set("HFB/indivOccupNeut", vecVn);
  result.set("HFB/indivOccupProt", vecVp);
  result.set("HFB/indivEnerNeut", eneQPn);
  result.set("HFB/indivEnerProt", eneQPp);
  result.set("HFB/chemPot", chemPot);

  result.set("HFB/eneTot", eneTot);

  IMAT oaiHFn = arma::zeros<IMAT >(vecOmegan.n_rows, 2);
  IMAT oaiHFp = arma::zeros<IMAT >(vecOmegap.n_rows, 2);
  oaiHFn.col(0) = vecOmegan;
  oaiHFp.col(0) = vecOmegap;
  result.set("HFB/indivOmegaNeut", oaiHFn);
  result.set("HFB/indivOmegaProt", oaiHFp);
  result.set("berger/cvg", convergence);
  result.set("berger/iterx", niterx);
  result.set("berger/itert", nitert);
  result.set("berger/qsciss", qsciss);
  result.set("berger/parent", parent);
  result.set("misc/jobname", jobname);

  for (INT i = 0; i < nbcon; i++)
  {
    if      (constraint_type[i] == "Q10") result.set("constraints/q10t", constraint_val[i]);
    else if (constraint_type[i] == "q20") result.set("constraints/q20t", constraint_val[i]);
    else if (constraint_type[i] == "q30") result.set("constraints/q30t", constraint_val[i]);
    else if (constraint_type[i] == "q40") result.set("constraints/q40t", constraint_val[i]);
    else if (constraint_type[i] == "q50") result.set("constraints/q50t", constraint_val[i]);
    else if (constraint_type[i] == "q60") result.set("constraints/q60t", constraint_val[i]);
    else    ERROR("Unknown constraint type: '" + constraint_type[i] + "'");

    if      (constraint_type[i] == "Q10") result.set("constraints_lambda/q10t", constraint_lambda[i]);
    else if (constraint_type[i] == "q20") result.set("constraints_lambda/q20t", constraint_lambda[i]);
    else if (constraint_type[i] == "q30") result.set("constraints_lambda/q30t", constraint_lambda[i]);
    else if (constraint_type[i] == "q40") result.set("constraints_lambda/q40t", constraint_lambda[i]);
    else if (constraint_type[i] == "q50") result.set("constraints_lambda/q50t", constraint_lambda[i]);
    else if (constraint_type[i] == "q60") result.set("constraints_lambda/q60t", constraint_lambda[i]);
    else    ERROR("Unknown constraint type: '" + constraint_type[i] + "'");
  }

  for (auto &ene: energiesNeut)
  {
    result.set("berger/enerNeut_" + ene.first, ene.second);
  }
  for (auto &ene: energiesProt)
  {
    result.set("berger/enerProt_" + ene.first, ene.second);
  }

  result.set("berger/enerCorrections", corrections);

  if (inertia0.size() != 0)
  {
    result.set("berger/inertia0", inertia0);
    result.set("berger/inertia23", inertia23);
    result.set("berger/inertia24", inertia24);
    result.set("berger/inertia34", inertia34);
    result.set("berger/inertia234", inertia234);
  }

  // TODO: implement the exact transformation instead

  /* Tools::mesg("IOberg", "will now perform one iter. of SolverHFBBroyden"); */
  /* SolverHFBBroyden solverHFBBroyden(result); */
  /* solverHFBBroyden.nextIter(); */
  /* result.merge(solverHFBBroyden.state.getDataTree()); */

  DBG_RETURN(result);
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Check file type
 */

bool IOberger::checkFileType(const std::string &_content)
{
  DBG_ENTER;

  if (_content.substr(0, 10) != "#noe000001")
  {
    DBG_RETURN(false);
  }

  DBG_RETURN(true);
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Save a DataTree to a .rhl file.
 */

void IOberger::saveState(const State &state, const std::string &filename)
{
  DBG_ENTER;

  //============================================================================
  //= Calculate sizes ==========================================================
  //============================================================================

  INT id0 = 0;
  INT id1 = 0;

  for (INT om = 0; om < state.basis.mMax; om++)
  {
    for (INT sd = 0; sd < state.basis.sMax; sd++)
    {
      INT md = om + sd;

      if (md >= state.basis.mMax) continue;

      for (INT sb = 0; sb < state.basis.sMax; sb++)
      {
        INT mb = om + sb;

        if (mb >= state.basis.mMax) continue;

        if (mb < md) continue;

        for (INT nb = 0; nb < state.basis.nMax(mb); nb++)
        {
          INT ndMax = (md == mb) ? nb + 1 : state.basis.nMax(md);

          for (INT nd = 0; nd < ndMax; nd++)
          {
            for (INT n_zb = 0; n_zb < state.basis.n_zMax(mb, nb); n_zb++)
            {
              INT n_zdMax = ((md == mb) && (nd == nb)) ? n_zb + 1 : state.basis.n_zMax(md, nd);

              for (INT n_zd = 0; n_zd < n_zdMax; n_zd++)
              {

                if (state.basis.dMax == 2)
                {
                  {
                    //                    INT dd = 1;
                    //                    INT db = 1;
                    id0++;
                    //                    INFO("ID0: %04d (%2d %2d %2d %2d %2d) (%2d %2d %2d %2d %2d)", id0, md, nd, n_zd, dd, sd, mb, nb, n_zb, db, sb);
                  }
                  {
                    //                    INT dd = 0;
                    //                    INT db = 0;
                    id1++;
                    //                    INFO("ID1: %04d (%2d %2d %2d %2d %2d) (%2d %2d %2d %2d %2d)", id1, md, nd, n_zd, dd, sd, mb, nb, n_zb, db, sb);
                  }
                }
                else
                {
                  {
                    //                    INT dd = 0;
                    //                    INT db = 0;
                    id0++;
                    //                    INFO("ID0: %04d (%2d %2d %2d %2d %2d) (%2d %2d %2d %2d %2d)", id0, md, nd, n_zd, dd, sd, mb, nb, n_zb, db, sb);
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  INT id2 = 0;

  if (state.basis.dMax == 2)
  {
    for (INT om = 0; om < state.basis.mMax; om++)
    {
      for (INT sd = 0; sd < state.basis.sMax; sd++)
      {
        INT md = om + sd;

        if (md >= state.basis.mMax) continue;

        for (INT sb = 0; sb < state.basis.sMax; sb++)
        {
          INT mb = om + sb;

          if (mb >= state.basis.mMax) continue;

          if (mb < md) continue;

          for (INT nb = 0; nb < state.basis.nMax(mb); nb++)
          {

            INT ndMax = (md == mb) ? nb + 1 : state.basis.nMax(md);

            for (INT nd = 0; nd < ndMax; nd++)
            {
              for (INT db = 0; db < state.basis.dMax; db++)
              {
                //                INT dd = 1 - db;
                if ((md == mb) && (nd == nb) && (db == 0)) continue;

                for (INT n_zb = 0; n_zb < state.basis.n_zMax(mb, nb); n_zb++)
                {

                  INT n_zdMax = state.basis.n_zMax(md, nd);

                  for (INT n_zd = 0; n_zd < n_zdMax; n_zd++)
                  {
                    id2++;

                    //                  INFO("ID2: %04d (%2d %2d %2d %2d %2d) (%2d %2d %2d %2d %2d)", id2, md, nd, n_zd, dd, sd, mb, nb, n_zb, db, sb);
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  //============================================================================
  //= Allocate vectors =========================================================
  //============================================================================

  Tools::debug(PF("IOberger: Berger2ct sizes: %d %d %d", id0, id1, id2));

  arma::vec ro0n = arma::zeros(id0);
  arma::vec ro1n = arma::zeros(id1);
  arma::vec ro2n = arma::zeros(id2);
  arma::vec ro0p = arma::zeros(id0);
  arma::vec ro1p = arma::zeros(id1);
  arma::vec ro2p = arma::zeros(id2);
  arma::vec xcp0n = arma::zeros(id0);
  arma::vec xcp1n = arma::zeros(id1);
  arma::vec xcp2n = arma::zeros(id2);
  arma::vec xcp0p = arma::zeros(id0);
  arma::vec xcp1p = arma::zeros(id1);
  arma::vec xcp2p = arma::zeros(id2);

  //============================================================================
  //= Fill vectors =============================================================
  //============================================================================

  INT i0 = 0;
  INT i1 = 0;

  for (INT om = 0; om < state.basis.mMax; om++)
  {
    for (INT sd = 0; sd < state.basis.sMax; sd++)
    {
      INT md = om + sd;

      if (md >= state.basis.mMax) continue;

      for (INT sb = 0; sb < state.basis.sMax; sb++)
      {
        INT mb = om + sb;

        if (mb >= state.basis.mMax) continue;

        if (mb < md) continue;

        for (INT nb = 0; nb < state.basis.nMax(mb); nb++)
        {
          INT ndMax = (md == mb) ? nb + 1 : state.basis.nMax(md);

          for (INT nd = 0; nd < ndMax; nd++)
          {
            for (INT n_zb = 0; n_zb < state.basis.n_zMax(mb, nb); n_zb++)
            {
              INT n_zdMax = ((md == mb) && (nd == nb)) ? n_zb + 1 : state.basis.n_zMax(md, nd);

              for (INT n_zd = 0; n_zd < n_zdMax; n_zd++)
              {

                if (state.basis.dMax == 2)
                {
                  {
                    INT dd = 1;
                    INT db = 1;

                    INT b = state.basis.HOqn.find({mb, nb, n_zb, db, sb});
                    INT d = state.basis.HOqn.find({md, nd, n_zd, dd, sd});
                    ro0n(i0) = state.rho(NEUTRON)(b, d);
                    ro0p(i0) = state.rho(PROTON )(b, d);
                    xcp0n(i0) = state.kappa(NEUTRON)(b, d);
                    xcp0p(i0) = state.kappa(PROTON )(b, d);
                    i0++;
                  }
                  {
                    INT dd = 0;
                    INT db = 0;

                    INT b = state.basis.HOqn.find({mb, nb, n_zb, db, sb});
                    INT d = state.basis.HOqn.find({md, nd, n_zd, dd, sd});
                    ro1n(i1) = state.rho(NEUTRON)(b, d);
                    ro1p(i1) = state.rho(PROTON )(b, d);
                    xcp1n(i1) = state.kappa(NEUTRON)(b, d);
                    xcp1p(i1) = state.kappa(PROTON )(b, d);
                    i1++;
                  }
                }
                else
                {
                  {
                    INT dd = 0;
                    INT db = 0;

                    INT b = state.basis.HOqn.find({mb, nb, n_zb, db, sb});
                    INT d = state.basis.HOqn.find({md, nd, n_zd, dd, sd});
                    ro0n(i0) = state.rho(NEUTRON)(b, d);
                    ro0p(i0) = state.rho(PROTON )(b, d);
                    xcp0n(i0) = state.kappa(NEUTRON)(b, d);
                    xcp0p(i0) = state.kappa(PROTON )(b, d);
                    //                    INFO("i0: %d b: %d d: %d rh0n: %e", i0, b, d, ro0n(i0));
                    i0++;
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  INT i2 = 0;

  if (state.basis.dMax == 2)
  {
    for (INT om = 0; om < state.basis.mMax; om++)
    {
      for (INT sd = 0; sd < state.basis.sMax; sd++)
      {
        INT md = om + sd;

        if (md >= state.basis.mMax) continue;

        for (INT sb = 0; sb < state.basis.sMax; sb++)
        {
          INT mb = om + sb;

          if (mb >= state.basis.mMax) continue;

          if (mb < md) continue;

          for (INT nb = 0; nb < state.basis.nMax(mb); nb++)
          {

            INT ndMax = (md == mb) ? nb + 1 : state.basis.nMax(md);

            for (INT nd = 0; nd < ndMax; nd++)
            {
              for (INT db = 0; db < state.basis.dMax; db++)
              {
                INT dd = 1 - db;

                if ((md == mb) && (nd == nb) && (db == 0)) continue;

                for (INT n_zb = 0; n_zb < state.basis.n_zMax(mb, nb); n_zb++)
                {

                  INT n_zdMax = state.basis.n_zMax(md, nd);
                  INT b = state.basis.HOqn.find({mb, nb, n_zb, db, sb});

                  for (INT n_zd = 0; n_zd < n_zdMax; n_zd++)
                  {
                    INT d = state.basis.HOqn.find({md, nd, n_zd, dd, sd});
                    ro2n(i2) = state.rho(NEUTRON)(b, d);
                    ro2p(i2) = state.rho(PROTON )(b, d);
                    xcp2n(i2) = state.kappa(NEUTRON)(b, d);
                    xcp2p(i2) = state.kappa(PROTON )(b, d);
                    i2++;
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  //============================================================================
  //= Write file ===============================================================
  //============================================================================

  std::ofstream fp(filename.c_str(), std::ios::binary);

  bufferAppend(230507); // magic number
  bufferAppend(state.basis.mMax); // mox

  bufferAppend(                    id1); // nro1
  bufferAppend(                    id2); // nro2
  bufferAppend(                    id0); // nro0
  bufferAppend((state.basis.dMax == 2) ? state.basis.mMax : 0); // mx1
  bufferAppend(    state.basis.mMax); // mx2
  bufferAppend(    state.basis.mMax + 1); // mx21
  bufferAppend((INT)state.basis.omegaIndexHO(0).n_rows); // jdxi
  bufferAppend(                      0); // vbidon
  bufferAppend(2 - state.basis.dMax); // v1ct
  bufferAppend(state.basis.dMax - 1); // v2ct
  bufferAppend(                      0); // vsyimp
  bufferAppend(  state.basis.nOscil); // nn
  bufferAppend(     state.basis.g_q); // gq
  bufferAppend(     state.basis.d_0); // d1
  bufferAppend(state.basis.hwBerger); // om1
  bufferAppend( state.basis.qBerger); // q1

  for (INT i = 0; i < state.basis.mMax; i++)
  {
    bufferAppend((INT)state.basis.n_zMaxBerger(i)); // nztu(i)
  }

  bufferAppend(                    state.basis.alfBerger); // alf
  bufferAppend(                    state.basis.betBerger); // bet
  bufferAppend(                      0); // nnl
  bufferAppend(                      0); // nnh
  bufferAppend(                      0); // nwrit
  bufferAppend(                      0); // nwrit1
  bufferAppend(std::string("01234567")); // corps
  bufferAppend((double)(state.sys.nProt + state.sys.nNeut)); // xba
  bufferAppend((double)(state.sys.nProt)); // xbz
  bufferAppend((double)(state.sys.nNeut)); // xbn
  bufferAppend(                    0.0); // xba1
  bufferAppend(                    0.0); // xbz1
  bufferAppend(                    0.0); // xbn1
  bufferAppend(                    0.0); // xbn2
  bufferAppend(                    0.0); // xbz2
  bufferAppend(                      0); // in12
  bufferAppend(                      1); // jbcs
  bufferAppend(                      0); // vzpe
  bufferAppend(                    0.0); // tempr
  bufferAppend(                      0); // inter
  bufferAppend(                    0.0); // qport
  bufferAppend(state.constraints.size()); // nbcoit

  for (auto &c : state.constraints)
  {
    bufferAppend(c.second.lambda); // xlamci
  }

  for (UINT i = state.constraints.size(); i < 9; i++)
  {
    bufferAppend(0.0); // xlamci
  }

  if (state.chemPot.n_rows > 1)
  {
    bufferAppend(state.chemPot(NEUTRON)); // xlmb(1,2) ??
    bufferAppend(state.chemPot( PROTON)); // xlmb(2,2) ??
  }
  else
  {
    bufferAppend(          -5.0); // xlmb(1,2) ??
    bufferAppend(          -5.0); // xlmb(1,2) ??
  }

  bufferAppend(                    9.999); // xlmb(1,1) ??
  bufferAppend(                    9.999); // xlmb(2,1) ??

  binWrite(fp, " (tableau 1)");

  bufferAppend(-1);
  binWrite(fp, " (tableau 2)");
  bufferAppend(-1);
  binWrite(fp, " (tableau 3)");
  bufferAppend(0);
  binWrite(fp, " (tableau 4)");
  bufferAppend(-1);
  binWrite(fp, " (tableau 5)");

  if (state.basis.dMax == 1)
  {
    for (INT i = 0; i < id0; i++) bufferAppend(ro0n(i));

    for (INT i = 0; i < id0; i++) bufferAppend(ro0p(i));

    for (INT i = 0; i < id0; i++) bufferAppend(xcp0n(i));

    for (INT i = 0; i < id0; i++) bufferAppend(xcp0p(i));

    binWrite(fp, " (ro0, ro0w, xcp0, xcp0w)");
  }
  else
  {
    for (INT i = 0; i < id1; i++) bufferAppend(ro1n(i));

    for (INT i = 0; i < id1; i++) bufferAppend(ro1p(i));

    for (INT i = 0; i < id2; i++) bufferAppend(ro2n(i));

    for (INT i = 0; i < id2; i++) bufferAppend(ro2p(i));

    for (INT i = 0; i < id0; i++) bufferAppend(ro0n(i));

    for (INT i = 0; i < id0; i++) bufferAppend(ro0p(i));

    for (INT i = 0; i < id1; i++) bufferAppend(xcp1n(i));

    for (INT i = 0; i < id1; i++) bufferAppend(xcp1p(i));

    for (INT i = 0; i < id2; i++) bufferAppend(xcp2n(i));

    for (INT i = 0; i < id2; i++) bufferAppend(xcp2p(i));

    for (INT i = 0; i < id0; i++) bufferAppend(xcp0n(i));

    for (INT i = 0; i < id0; i++) bufferAppend(xcp0p(i));

    binWrite(fp, " (ro1, ro1w, ro2, ro2w, ro0, ro0w, xcp1, xcp1w, xcp2, xcp2w, xcp0, xcp0w)");
  }

  fp.close();

  Tools::mesg("IOberg", "State saved in " + PF_GREEN(filename));

  DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Write the buffer to a binary file.
 */

void IOberger::binWrite(std::ofstream &fp, const std::string &mesg)
{
  Tools::debug(PF("IOberger: writing block of %d bytes", bufferSize) + mesg);

  fp.write(reinterpret_cast<char *>(&bufferSize), sizeof(INT));
  fp.write(buffer, bufferSize);
  fp.write(reinterpret_cast<char *>(&bufferSize), sizeof(INT));
  free(buffer);
  buffer = NULL;
  bufferSize = 0;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Add some value to the binary buffer.
 */

void IOberger::bufferAppend(const std::string &s)
{
  char *new_buffer = (char *)realloc(buffer, bufferSize + s.size());

  ASSERT(new_buffer != NULL, "can not allocate memory");

  buffer = new_buffer;

  unsigned char *ptr = (unsigned char *)(s.c_str());

  for (INT i = 0; i < s.size(); i++)
  {
    buffer[bufferSize + i] = *ptr;
    ptr ++;
  }

  bufferSize += s.size();
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Add some value to the binary buffer.
 */

template<class T>
void IOberger::bufferAppend(T val)
{
  char *new_buffer = (char *)realloc(buffer, bufferSize + sizeof(T));

  ASSERT(new_buffer != NULL, "can not allocate memory");

  buffer = new_buffer;

  unsigned char *ptr = (unsigned char *)(&val);

  for (INT i = 0; i < sizeof(T); i++)
  {
    buffer[bufferSize + i] = *ptr;
    ptr ++;
  }

  bufferSize += sizeof(T);
}
