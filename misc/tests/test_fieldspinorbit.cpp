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

/** \file
 *  \brief Test suite for the FieldSpinOrbit class.
 */

#include "gtest/gtest.h"
#include "field_spin_orbit.h"
#include "state.h"
#include "interaction.h"

//==============================================================================

TEST(FieldSpinOrbit, calcField)
{
  {
    //1ct state
    double err = 1e-6;
    DataTree dataTree = DataTree("examples/42Ca_deformed_1x11.msg.gz");

    dataTree.validate();

    State state(dataTree);

    Interaction interaction(dataTree, &state);

    auto fp = interaction.getParametersFromDataTree(dataTree, "spin-orbit")[0];
    fp["wso"] = 130.0;

    FieldSpinOrbit Sp(fp, &state);

    Sp.calcField();
    Sp.calcLegacy();

    arma::mat DifN = arma::abs(Sp.fieldLegacy(NEUTRON) - Sp.field(NEUTRON, Field::DIRECT)*(1+7.210290e-06));
    arma::mat DifP = arma::abs(Sp.fieldLegacy(PROTON ) - Sp.field(PROTON , Field::DIRECT)*(1+7.210290e-06));

    double dif = MAX(DifN.max(),DifP.max());

    ASSERT_NEAR(dif, 0, err);
  }
  {
    //2ct state
    double err = 1e-6;
    DataTree dataTree("examples/42Ca_deformed_2x9.msg.gz");
    State state(dataTree);

    Interaction interaction(dataTree, &state);

    auto fp = interaction.getParametersFromDataTree(dataTree, "spin-orbit")[0];
    fp["wso"] = 130.0;

    FieldSpinOrbit Sp(fp, &state);

    Sp.calcField();
    Sp.calcLegacy();

    arma::mat DifN = arma::abs(Sp.fieldLegacy(NEUTRON) - Sp.field(NEUTRON, Field::DIRECT)*(1+7.210290e-06));
    arma::mat DifP = arma::abs(Sp.fieldLegacy(PROTON ) - Sp.field(PROTON , Field::DIRECT)*(1+7.210290e-06));

    double dif = MAX(DifN.max(),DifP.max());

    ASSERT_NEAR(dif, 0, err);
  }
}

//==============================================================================

TEST(FieldSpinOrbit, calcXz)
{
  DataTree dataTree("examples/42Ca_deformed_2x9.msg.gz");
  State state(dataTree);

  Interaction interaction(dataTree, &state);

  auto fp = interaction.getParametersFromDataTree(dataTree, "spin-orbit")[0];
  fp["wso"] = 130.0;

  {
    double err = 1e-8;
    double Error = 0.0;
    State state ("examples/42Ca_deformed_2x9.msg.gz");
    FieldSpinOrbit Sp(fp, &state);

    Sp.calcXz();
    Sp.calcYz();

    state.basis.calcIz();
    INT dMax = state.basis.dMax;
    INT n_zMax = state.basis.n_zGlobalMax;

    //Tests Xz(1).
    for(INT d_alpha = 0 ; d_alpha < dMax ; d_alpha++)
    {
      for(INT d_gamma = 0 ; d_gamma < dMax ; d_gamma++)
      {
        for(INT n_zalpha = 0 ; n_zalpha < n_zMax ; n_zalpha++)
        {
          for(INT n_zgamma = 0 ; n_zgamma < n_zMax ; n_zgamma++)
          {
            for(INT n_za = 0 ; n_za < 2*n_zMax+2 ; n_za++)
            {
              for(INT d_beta = 0 ; d_beta < dMax ; d_beta++)
              {
                for(INT d_delta = 0 ; d_delta < dMax ; d_delta++)
                {
                  Error += abs(Sp.Xz(1,n_za,d_beta+2*d_delta)(n_zalpha,n_zgamma,d_alpha+2*d_gamma) - Sp.Xz(0,n_za,d_beta+2*d_delta)(n_zgamma,n_zalpha,d_gamma+2*d_alpha));
                }//d_delta
              }//d_beta
            }//n_za
          }//n_zgamma
        }//n_zalpha
      }//d_gamma
    }//d_alpha

    ASSERT_NEAR(Error, 0, err);
  }
}

//==============================================================================

TEST(FieldSpinOrbit, calcYz)
{
  DataTree dataTree("examples/42Ca_deformed_2x9.msg.gz");
  State state(dataTree);

  Interaction interaction(dataTree, &state);

  auto fp = interaction.getParametersFromDataTree(dataTree, "spin-orbit")[0];
  fp["wso"] = 130.0;

  {
    double err = 1e-6;
    double Error_1  = 0.0;
    double Error_20 = 0.0;
    double Error_22 = 0.0;

    State state ("examples/42Ca_deformed_2x9.msg.gz");
    FieldSpinOrbit Sp(fp, &state);

    Sp.calcXz();
    Sp.calcYz();
    state.basis.calcIz();
    INT dMax = state.basis.dMax;
    INT DMax = (dMax-1)*3+1;
    INT n_zMax = state.basis.n_zGlobalMax;

    //Reshapes Xz(2).
    Multi<arma::vec> xz_2;
    for(INT d_delta = 0 ; d_delta < dMax ; d_delta++)
    {
      for(INT d_beta = 0 ; d_beta < dMax ; d_beta++)
      {
        for(INT d_alpha = 0 ; d_alpha < dMax ; d_alpha++)
        {
          for(INT d_gamma = 0 ; d_gamma < dMax ; d_gamma++)
          {
            for(INT n_zbeta = 0 ; n_zbeta < n_zMax ; n_zbeta++)
            {
              for(INT n_zdelta = 0 ; n_zdelta < n_zMax ; n_zdelta++)
              {
                arma::vec &xvec = xz_2(n_zdelta,n_zbeta,d_delta+2*d_beta,d_alpha+2*d_gamma);
                xvec = arma::zeros(2*n_zMax+2);
                for(INT n_za = 0 ; n_za < 2*n_zMax +2 ; n_za++)
                {
                  xvec(n_za) = Sp.Xz(2,n_za,d_alpha+2*d_gamma)(n_zdelta,n_zbeta,d_delta+2*d_beta);
                }//n_za
              }//n_zdelta
            }//n_zbeta
          }//d_gamma
        }//d_alpha
      }//d_beta
    }//d_delta
    //Test Yz(0)
    Multi<arma::cube> Ixz_20;
    for(INT d_alpha = 0 ; d_alpha < dMax ; d_alpha++)
    {
      for(INT d_gamma = 0 ; d_gamma < dMax ; d_gamma++)
      {
        for(INT n_zbeta = 0 ; n_zbeta < n_zMax ; n_zbeta++)
        {
          for(INT n_zdelta = 0 ; n_zdelta < n_zMax ; n_zdelta++)
          {
            Ixz_20(n_zdelta,n_zbeta,d_alpha+2*d_gamma) = arma::zeros(n_zMax,n_zMax,DMax);
          }//n_zdelta
        }//n_zbeta
      }//d_gamma
    }//d_alpha
    for(INT d_delta = 0 ; d_delta < dMax ; d_delta++)
    {
      for(INT d_beta = 0 ; d_beta < dMax ; d_beta++)
      {
        for(INT d_alpha = 0 ; d_alpha < dMax ; d_alpha++)
        {
          for(INT d_gamma = 0 ; d_gamma < dMax ; d_gamma++)
          {
            for(INT n_zbeta = 0 ; n_zbeta < n_zMax ; n_zbeta++)
            {
              for(INT n_zdelta = 0 ; n_zdelta < n_zMax ; n_zdelta++)
              {
                arma::vec &xvec = xz_2(n_zdelta,n_zbeta,d_delta+2*d_beta,d_alpha+2*d_gamma);
                for(INT n_zalpha = 0 ; n_zalpha < n_zMax ; n_zalpha++)
                {
                  for(INT n_zgamma = 0 ; n_zgamma < n_zMax ; n_zgamma++)
                  {
                    Ixz_20(n_zalpha,n_zgamma,d_alpha+2*d_gamma)(n_zdelta,n_zbeta,d_delta+2*d_beta) = arma::accu(Sp.Yz(0,n_zalpha,n_zgamma,d_alpha+2*d_gamma)%xvec);
                    Error_20 += abs(Ixz_20(n_zalpha,n_zgamma,d_alpha+2*d_gamma)(n_zdelta,n_zbeta,d_delta+2*d_beta) - state.basis.Iz(n_zalpha,n_zgamma,d_alpha+2*d_gamma)(n_zbeta,n_zdelta,d_beta+2*d_delta));
                  }//n_zgamma
                }//n_zalpha
              }//n_zdelta
            }//n_zbeta
          }//d_gamma
        }//d_alpha
      }//d_beta
    }//d_delta

    //Tests Yz(2).
    Multi<arma::cube> Ixz_22;
    for(INT d_alpha = 0 ; d_alpha < dMax ; d_alpha++)
    {
      for(INT d_gamma = 0 ; d_gamma < dMax ; d_gamma++)
      {
        for(INT n_zbeta = 0 ; n_zbeta < n_zMax ; n_zbeta++)
        {
          for(INT n_zdelta = 0 ; n_zdelta < n_zMax ; n_zdelta++)
          {
            Ixz_22(n_zdelta,n_zbeta,d_alpha+2*d_gamma) = arma::zeros(n_zMax,n_zMax,DMax);
          }//n_zdelta
        }//n_zbeta
      }//d_gamma
    }//d_alpha
    //Merges J-z and talmanz.
    for(INT d_delta = 0 ; d_delta < dMax ; d_delta++)
    {
      for(INT d_beta = 0 ; d_beta < dMax ; d_beta++)
      {
        for(INT d_alpha = 0 ; d_alpha < dMax ; d_alpha++)
        {
          for(INT d_gamma = 0 ; d_gamma < dMax ; d_gamma++)
          {
            for(INT n_zbeta = 0 ; n_zbeta < n_zMax ; n_zbeta++)
            {
              for(INT n_zdelta = 0 ; n_zdelta < n_zMax ; n_zdelta++)
              {
                arma::vec &xvec = xz_2(n_zdelta,n_zbeta,d_delta+2*d_beta,d_alpha+2*d_gamma);
                for(INT n_zalpha = 0 ; n_zalpha < n_zMax ; n_zalpha++)
                {
                  for(INT n_zgamma = 0 ; n_zgamma < n_zMax ; n_zgamma++)
                  {
                    Ixz_22(n_zdelta,n_zbeta,d_delta+2*d_beta)(n_zalpha,n_zgamma,d_alpha+2*d_gamma) = arma::accu(Sp.Yz(2,n_zalpha,n_zgamma,d_alpha+2*d_gamma)%xvec);
                    Error_22 += abs(Ixz_22(n_zdelta,n_zbeta,d_delta+2*d_beta)(n_zalpha,n_zgamma,d_alpha+2*d_gamma) - state.basis.NablaIz(n_zdelta,n_zbeta,d_delta+2*d_beta)(n_zalpha,n_zgamma,d_alpha+2*d_gamma));
                  }//n_zgamma
                }//n_zalpha
              }//n_zdelta
            }//n_zbeta
          }//d_gamma
        }//d_alpha
      }//d_beta
    }//d_delta

    //Tests Yz(1).
    for(INT d_alpha = 0 ; d_alpha < dMax ; d_alpha++)
    {
      for(INT d_gamma = 0 ; d_gamma < dMax ; d_gamma++)
      {
        for(INT n_zalpha = 0 ; n_zalpha < n_zMax ; n_zalpha++)
        {
          for(INT n_zgamma = 0 ; n_zgamma < n_zMax ; n_zgamma++)
          {
            Error_1 += arma::accu(arma::abs(Sp.Yz(1,n_zalpha,n_zgamma,d_alpha+2*d_gamma) - Sp.Yz(2,n_zgamma,n_zalpha,d_gamma+2*d_alpha)));
          }//n_zgamma
        }//n_zalpha
      }//d_gamma
    }//d_alpha

    double Error_Global = MAX(Error_22,Error_20);
    Error_Global = MAX(Error_Global,Error_1);
    ASSERT_NEAR(Error_Global, 0, err);
  }
}

//==============================================================================

TEST(FieldSpinOrbit, calcXr)
{
  DataTree dataTree("examples/42Ca_deformed_2x9.msg.gz");
  State state(dataTree);

  Interaction interaction(dataTree, &state);

  auto fp = interaction.getParametersFromDataTree(dataTree, "spin-orbit")[0];
  fp["wso"] = 130.0;

  {
    double err = 1e-10;
    double Error_0 = 0.0;
    double Error_1 = 0.0;
    double Error_2 = 0.0;
    double Error_3 = 0.0;

    State state ("examples/42Ca_deformed_2x9.msg.gz");
    FieldSpinOrbit Sp(fp, &state);

    Sp.calcXr();
    Sp.calcYr();
    state.basis.calcIr();

    INT maxn = INT(state.basis.moshinskyr(0,0).n_elem);
    INT mMax = state.basis.mMax;

    //Tests Xr(3).
    for(INT m = 0 ; m < mMax ; m++)
    {
      for(INT n_alpha = 0 ; n_alpha < state.basis.nMax(m) ; n_alpha++)
      {
       	for(INT n_gamma = 0 ; n_gamma < state.basis.nMax(m) ; n_gamma++)
        {
          for(INT n_a = 0 ; n_a < maxn ; n_a++)
          {
            Error_3 += abs(Sp.Xr(3,n_a,m)(n_alpha,n_gamma) - state.basis.talmanr(m,n_alpha,m,n_gamma)(n_a));
          }//n_a
	}//n_gamma
      }//n_alpha
    }//m
    //Tests Xr(0).
    for(INT m = 0 ; m < mMax-1 ; m++)
    {
      for(INT n_alpha = 0 ; n_alpha < state.basis.nMax(m) ; n_alpha++)
      {
       	for(INT n_gamma = 0 ; n_gamma < state.basis.nMax(m+1) ; n_gamma++)
        {
          for(INT n_a = 0 ; n_a < maxn ; n_a++)
          {
            Error_0 += abs(Sp.Xr(0,n_a,m)(n_alpha,n_gamma) - Sp.Yr(2,m+1,n_gamma,n_alpha)(n_a));
          }//n_a
	}//n_gamma
      }//n_alpha
    }//m
    //Tests Xr(1).
    for(INT m = 0 ; m < mMax-1 ; m++)
    {
      for(INT n_alpha = 0 ; n_alpha < state.basis.nMax(m) ; n_alpha++)
      {
       	for(INT n_gamma = 0 ; n_gamma < state.basis.nMax(m+1) ; n_gamma++)
        {
          for(INT n_a = 0 ; n_a < maxn ; n_a++)
          {
            Error_1 += abs(Sp.Xr(1,n_a,m)(n_alpha,n_gamma) - Sp.Yr(3,m+1,n_gamma,n_alpha)(n_a));
          }//n_a
	}//n_gamma
      }//n_alpha
    }//m
    //Tests Xr(2).
    for(INT m = 1 ; m < mMax ; m++)
    {
      for(INT n_alpha = 0 ; n_alpha < state.basis.nMax(m) ; n_alpha++)
      {
       	for(INT n_gamma = 0 ; n_gamma < state.basis.nMax(m) ; n_gamma++)
        {
          for(INT n_a = 0 ; n_a < maxn ; n_a++)
          {
            Error_2 += abs(Sp.Xr(2,n_a,m)(n_alpha,n_gamma) - Sp.Yr(1,m,n_gamma,n_alpha)(n_a));
          }//n_a
	}//n_gamma
      }//n_alpha
    }//m

    double Error_Global = MAX(Error_3,Error_2);
    Error_Global = MAX(Error_Global,Error_1);
    Error_Global = MAX(Error_Global,Error_0);
    ASSERT_NEAR(Error_Global, 0, err);

  }
}

//==============================================================================

TEST(FieldSpinOrbit, calcYr)
{
  DataTree dataTree("examples/42Ca_deformed_2x9.msg.gz");
  State state(dataTree);

  Interaction interaction(dataTree, &state);

  auto fp = interaction.getParametersFromDataTree(dataTree, "spin-orbit")[0];
  fp["wso"] = 130.000;

  {
    double err = 1e-10;
    double Error_32 = 0.0;
    double Error_33 = 0.0;
    double Error_31 = 0.0;
    double Error_20 = 0.0;

    State state ("examples/42Ca_deformed_2x9.msg.gz");
    FieldSpinOrbit Sp(fp, &state);

    Sp.calcXr();
    Sp.calcYr();
    state.basis.calcIr();

    INT maxn = INT(state.basis.moshinskyr(0,0).n_elem);
    INT mMax = state.basis.mMax;

    //Reshapes Xr(3).
    Multi<arma::vec> xr_3;
    for(INT m = 0 ; m < mMax ; m++)
    {
      for(INT n_alpha = 0 ; n_alpha < state.basis.nMax(m) ; n_alpha++)
      {
        for(INT n_gamma = 0 ; n_gamma < state.basis.nMax(m) ; n_gamma++)
        {
          xr_3(m,n_alpha,n_gamma) = arma::zeros(maxn);
          for(INT n_a = 0 ; n_a < maxn ; n_a++)
          {
            xr_3(m,n_alpha,n_gamma)(n_a) = Sp.Xr(3,n_a,m)(n_alpha,n_gamma);
          }//n_a
        }//n_gamma
      }//n_alpha
    }//m

    //Tests Yr(2).
    Multi<arma::mat> Ixr_32;
    for(INT m = 0 ; m < mMax ; m++)
    {
      for(INT mp = 1 ; mp < mMax ; mp++)
      {
        for(INT n_alpha = 0 ; n_alpha < state.basis.nMax(m) ; n_alpha++)
        {
          for(INT n_gamma = 0 ; n_gamma < state.basis.nMax(m) ; n_gamma++)
          {
            Ixr_32(m,n_alpha,n_gamma,mp) = arma::zeros(state.basis.nMax(mp),state.basis.nMax(mp-1));
            for(INT n_delta = 0 ; n_delta < state.basis.nMax(mp) ; n_delta++)
            {
              for(INT n_beta = 0 ; n_beta < state.basis.nMax(mp-1) ; n_beta++)
              {
                Ixr_32(m,n_alpha,n_gamma,mp)(n_delta,n_beta) = arma::accu(xr_3(m,n_alpha,n_gamma)%Sp.Yr(2,mp,n_delta,n_beta));
              }//n_beta
            }//n_delta
            arma::mat Dif = arma::abs(Ixr_32(m,n_alpha,n_gamma,mp) - state.basis.Ir_plus(m,n_alpha,n_gamma,mp));
            Error_32 += Dif.max();
          }//n_gamma
        }//n_alpha
      }//mp
    }//m

    //Tests Yr(3).
    Multi<arma::mat> Ixr_33;
    for(INT m = 0 ; m < mMax ; m++)
    {
      for(INT mp = 0 ; mp < mMax-1 ; mp++)
      {
        for(INT n_alpha = 0 ; n_alpha < state.basis.nMax(m) ; n_alpha++)
        {
          for(INT n_gamma = 0 ; n_gamma < state.basis.nMax(m) ; n_gamma++)
          {
            Ixr_33(m,n_alpha,n_gamma,mp) = arma::zeros(state.basis.nMax(mp),state.basis.nMax(mp+1));
            for(INT n_delta = 0 ; n_delta < state.basis.nMax(mp) ; n_delta++)
            {
              for(INT n_beta = 0 ; n_beta < state.basis.nMax(mp+1) ; n_beta++)
              {
                Ixr_33(m,n_alpha,n_gamma,mp)(n_delta,n_beta) = arma::accu(xr_3(m,n_alpha,n_gamma)%Sp.Yr(3,mp+1,n_beta,n_delta));
              }//n_beta
            }//n_delta
            arma::mat Dif = arma::abs(Ixr_33(m,n_alpha,n_gamma,mp) - state.basis.Ir_minus(m,n_alpha,n_gamma,mp));
            Error_33 += Dif.max();
          }//n_gamma
        }//n_alpha
      }//mp
    }//m

    //Tests Yr(1).
    Multi<arma::mat> Ixr_31;
    for(INT m = 0 ; m < mMax ; m++)
    {
      for(INT mp = 1 ; mp < mMax ; mp++)
      {
        for(INT n_alpha = 0 ; n_alpha < state.basis.nMax(m) ; n_alpha++)
        {
          for(INT n_gamma = 0 ; n_gamma < state.basis.nMax(m) ; n_gamma++)
          {
            Ixr_31(m,n_alpha,n_gamma,mp) = arma::zeros(state.basis.nMax(mp),state.basis.nMax(mp));
            for(INT n_delta = 0 ; n_delta < state.basis.nMax(mp) ; n_delta++)
            {
              for(INT n_beta = 0 ; n_beta < state.basis.nMax(mp) ; n_beta++)
              {
                Ixr_31(m,n_alpha,n_gamma,mp)(n_delta,n_beta) = arma::accu(xr_3(m,n_alpha,n_gamma)%Sp.Yr(1,mp,n_beta,n_delta));
              }//n_beta
            }//n_delta
            arma::mat SumI = state.basis.Ir_mm(m,n_alpha,n_gamma,mp-1) - state.basis.Ir_pp(m,n_alpha,n_gamma,mp+1);
            arma::mat Dif = arma::abs(Ixr_31(m,n_alpha,n_gamma,mp) - SumI);
            Error_31 += Dif.max();
          }//n_gamma
        }//n_alpha
      }//mp
    }//m

    //Reshapes Xr(2);
    Multi<arma::vec> xr_2;
    for(INT m = 1 ; m < mMax ; m++)
    {
      for(INT n_alpha = 0 ; n_alpha < state.basis.nMax(m) ; n_alpha++)
      {
        for(INT n_gamma = 0 ; n_gamma < state.basis.nMax(m) ; n_gamma++)
        {
          xr_2(m,n_alpha,n_gamma) = arma::zeros(maxn);
          for(INT n_a = 0 ; n_a < maxn ; n_a++)
          {
            xr_2(m,n_alpha,n_gamma)(n_a) = Sp.Xr(2,n_a,m)(n_alpha,n_gamma);
          }//n_a
        }//n_gamma
      }//n_alpha
    }//m

    //Tests Yr(0).
    Multi<arma::mat> Ixr_20;
    for(INT m = 0 ; m < mMax ; m++)
    {
      for(INT mp = 1 ; mp < mMax ; mp++)
      {
        for(INT n_alpha = 0 ; n_alpha < state.basis.nMax(m) ; n_alpha++)
        {
          for(INT n_gamma = 0 ; n_gamma < state.basis.nMax(m) ; n_gamma++)
          {
            Ixr_20(m,n_alpha,n_gamma,mp) = arma::zeros(state.basis.nMax(mp),state.basis.nMax(mp));
            for(INT n_delta = 0 ; n_delta < state.basis.nMax(mp) ; n_delta++)
            {
              for(INT n_beta = 0 ; n_beta < state.basis.nMax(mp) ; n_beta++)
              {
                Ixr_20(m,n_alpha,n_gamma,mp)(n_beta,n_delta) = arma::accu(xr_2(mp,n_beta,n_delta)%Sp.Yr(0,m,n_alpha,n_gamma));
              }//n_beta
            }//n_delta
            arma::mat SumI = state.basis.Ir_mm(m,n_alpha,n_gamma,mp-1) - state.basis.Ir_pp(m,n_alpha,n_gamma,mp+1);
            arma::mat Dif = arma::abs(Ixr_20(m,n_alpha,n_gamma,mp) - SumI);
            Error_20 += Dif.max();
          }//n_gamma
        }//n_alpha
      }//mp
    }//m

    double Error_Global = MAX(Error_31,Error_32);
    Error_Global = MAX(Error_Global,Error_33);
    Error_Global = MAX(Error_Global,Error_20);
    ASSERT_NEAR(Error_Global, 0, err);

  }
}


