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
 *  \brief Test suite for the MultipoleOperators class.
 */

#include "gtest/gtest.h"
#include "multipole_operators.h"
#include "maxima_targets.h"

//==============================================================================

TEST(MultipoleOperators, calcQl0Matrices)
{
  State state("examples/42Ca_deformed_1x11.msg.gz");
  MultipoleOperators multipoleOperators(state);

  double err = 1e-09;

  ASSERT_NEAR(multipoleOperators.matZk(0)(state.basis.HOqn.find({0, 0, 8, 0, 0}), state.basis.HOqn.find({0, 0, 9, 0, 0})), 0.000000000000000000000000000000, err);
  ASSERT_NEAR(multipoleOperators.matZk(0)(state.basis.HOqn.find({0, 0, 8, 0, 0}), state.basis.HOqn.find({0, 0,10, 0, 0})), 0.000000000000000000000000000000, err);
  ASSERT_NEAR(multipoleOperators.matZk(0)(state.basis.HOqn.find({0, 0,10, 0, 0}), state.basis.HOqn.find({0, 0,10, 0, 0})), 1.000000000000000000000000000000, err);

  ASSERT_NEAR(multipoleOperators.matZk(1)(state.basis.HOqn.find({0, 0, 8, 0, 0}), state.basis.HOqn.find({0, 0, 9, 0, 0})), MULTIPOLE_OPERATORS_CALCQL0MATRICES_TARGET00, err);
  ASSERT_NEAR(multipoleOperators.matZk(1)(state.basis.HOqn.find({0, 0, 8, 0, 0}), state.basis.HOqn.find({0, 0, 0, 0, 0})), MULTIPOLE_OPERATORS_CALCQL0MATRICES_TARGET01, err);
  ASSERT_NEAR(multipoleOperators.matZk(1)(state.basis.HOqn.find({0, 0, 1, 0, 0}), state.basis.HOqn.find({0, 0, 1, 0, 0})), MULTIPOLE_OPERATORS_CALCQL0MATRICES_TARGET02, err);
  ASSERT_NEAR(multipoleOperators.matZk(2)(state.basis.HOqn.find({0, 0, 1, 0, 0}), state.basis.HOqn.find({0, 0, 2, 0, 0})), MULTIPOLE_OPERATORS_CALCQL0MATRICES_TARGET03, err);
  ASSERT_NEAR(multipoleOperators.matZk(2)(state.basis.HOqn.find({0, 0, 1, 0, 0}), state.basis.HOqn.find({0, 0, 2, 0, 0})), MULTIPOLE_OPERATORS_CALCQL0MATRICES_TARGET04, err);
  ASSERT_NEAR(multipoleOperators.matZk(2)(state.basis.HOqn.find({0, 0, 1, 0, 0}), state.basis.HOqn.find({0, 0, 3, 0, 0})), MULTIPOLE_OPERATORS_CALCQL0MATRICES_TARGET05, err);
  ASSERT_NEAR(multipoleOperators.matZk(3)(state.basis.HOqn.find({0, 0, 2, 0, 0}), state.basis.HOqn.find({0, 0, 5, 0, 0})), MULTIPOLE_OPERATORS_CALCQL0MATRICES_TARGET06, err);
  ASSERT_NEAR(multipoleOperators.matZk(3)(state.basis.HOqn.find({0, 0, 2, 0, 0}), state.basis.HOqn.find({0, 0, 6, 0, 0})), MULTIPOLE_OPERATORS_CALCQL0MATRICES_TARGET07, err);
  ASSERT_NEAR(multipoleOperators.matZk(3)(state.basis.HOqn.find({0, 0, 3, 0, 0}), state.basis.HOqn.find({0, 0, 4, 0, 0})), MULTIPOLE_OPERATORS_CALCQL0MATRICES_TARGET08, err);
  ASSERT_NEAR(multipoleOperators.matZk(4)(state.basis.HOqn.find({0, 0, 3, 0, 0}), state.basis.HOqn.find({0, 0, 5, 0, 0})), MULTIPOLE_OPERATORS_CALCQL0MATRICES_TARGET09, err);
  ASSERT_NEAR(multipoleOperators.matZk(4)(state.basis.HOqn.find({0, 0, 7, 0, 0}), state.basis.HOqn.find({0, 0, 2, 0, 0})), MULTIPOLE_OPERATORS_CALCQL0MATRICES_TARGET10, err);
  ASSERT_NEAR(multipoleOperators.matZk(4)(state.basis.HOqn.find({0, 0, 7, 0, 0}), state.basis.HOqn.find({0, 0, 3, 0, 0})), MULTIPOLE_OPERATORS_CALCQL0MATRICES_TARGET11, err);
  ASSERT_NEAR(multipoleOperators.matZk(5)(state.basis.HOqn.find({0, 0, 1, 0, 0}), state.basis.HOqn.find({0, 0, 2, 0, 0})), MULTIPOLE_OPERATORS_CALCQL0MATRICES_TARGET12, err);
  ASSERT_NEAR(multipoleOperators.matZk(5)(state.basis.HOqn.find({0, 0, 1, 0, 0}), state.basis.HOqn.find({0, 0, 3, 0, 0})), MULTIPOLE_OPERATORS_CALCQL0MATRICES_TARGET13, err);
  ASSERT_NEAR(multipoleOperators.matZk(5)(state.basis.HOqn.find({0, 0, 0, 0, 0}), state.basis.HOqn.find({0, 0, 3, 0, 0})), MULTIPOLE_OPERATORS_CALCQL0MATRICES_TARGET14, err);
  ASSERT_NEAR(multipoleOperators.matZk(6)(state.basis.HOqn.find({0, 0, 0, 0, 0}), state.basis.HOqn.find({0, 0, 0, 0, 0})), MULTIPOLE_OPERATORS_CALCQL0MATRICES_TARGET15, err);
  ASSERT_NEAR(multipoleOperators.matZk(6)(state.basis.HOqn.find({0, 0, 0, 0, 0}), state.basis.HOqn.find({0, 0, 1, 0, 0})), MULTIPOLE_OPERATORS_CALCQL0MATRICES_TARGET16, err);

  ASSERT_NEAR(multipoleOperators.matR2k(1)(state.basis.HOqn.find({ 0, 1, 0, 0, 0}), state.basis.HOqn.find({ 0, 1, 0, 0, 0})), MULTIPOLE_OPERATORS_CALCQL0MATRICES_TARGET17, err);
  ASSERT_NEAR(multipoleOperators.matR2k(1)(state.basis.HOqn.find({ 0, 1, 0, 0, 0}), state.basis.HOqn.find({ 0, 2, 0, 0, 0})), MULTIPOLE_OPERATORS_CALCQL0MATRICES_TARGET18, err);
  ASSERT_NEAR(multipoleOperators.matR2k(1)(state.basis.HOqn.find({ 8, 0, 0, 0, 0}), state.basis.HOqn.find({ 8, 1, 0, 0, 0})), MULTIPOLE_OPERATORS_CALCQL0MATRICES_TARGET19, err);
  ASSERT_NEAR(multipoleOperators.matR2k(1)(state.basis.HOqn.find({10, 0, 0, 0, 0}), state.basis.HOqn.find({10, 0, 0, 0, 0})), MULTIPOLE_OPERATORS_CALCQL0MATRICES_TARGET20, err);
  ASSERT_NEAR(multipoleOperators.matR2k(2)(state.basis.HOqn.find({ 0, 1, 0, 0, 0}), state.basis.HOqn.find({ 0, 1, 0, 0, 0})), MULTIPOLE_OPERATORS_CALCQL0MATRICES_TARGET21, err);
  ASSERT_NEAR(multipoleOperators.matR2k(2)(state.basis.HOqn.find({ 0, 1, 0, 0, 0}), state.basis.HOqn.find({ 0, 2, 0, 0, 0})), MULTIPOLE_OPERATORS_CALCQL0MATRICES_TARGET22, err);
  ASSERT_NEAR(multipoleOperators.matR2k(2)(state.basis.HOqn.find({ 8, 0, 0, 0, 0}), state.basis.HOqn.find({ 8, 1, 0, 0, 0})), MULTIPOLE_OPERATORS_CALCQL0MATRICES_TARGET23, err);
  ASSERT_NEAR(multipoleOperators.matR2k(2)(state.basis.HOqn.find({10, 0, 0, 0, 0}), state.basis.HOqn.find({10, 0, 0, 0, 0})), MULTIPOLE_OPERATORS_CALCQL0MATRICES_TARGET24, err);
  ASSERT_NEAR(multipoleOperators.matR2k(3)(state.basis.HOqn.find({ 0, 1, 0, 0, 0}), state.basis.HOqn.find({ 0, 1, 0, 0, 0})), MULTIPOLE_OPERATORS_CALCQL0MATRICES_TARGET25, err);
  ASSERT_NEAR(multipoleOperators.matR2k(3)(state.basis.HOqn.find({ 0, 1, 0, 0, 0}), state.basis.HOqn.find({ 0, 2, 0, 0, 0})), MULTIPOLE_OPERATORS_CALCQL0MATRICES_TARGET26, err);
  ASSERT_NEAR(multipoleOperators.matR2k(3)(state.basis.HOqn.find({ 8, 0, 0, 0, 0}), state.basis.HOqn.find({ 8, 1, 0, 0, 0})), MULTIPOLE_OPERATORS_CALCQL0MATRICES_TARGET27, err);
  ASSERT_NEAR(multipoleOperators.matR2k(3)(state.basis.HOqn.find({10, 0, 0, 0, 0}), state.basis.HOqn.find({10, 0, 0, 0, 0})), MULTIPOLE_OPERATORS_CALCQL0MATRICES_TARGET28, err);

  ASSERT_NEAR(multipoleOperators.ql0(2)(state.basis.HOqn.find({0, 1, 8, 0, 0}), state.basis.HOqn.find({0, 1, 9, 0, 0})), MULTIPOLE_OPERATORS_CALCQL0MATRICES_TARGET29, err);
  ASSERT_NEAR(multipoleOperators.ql0(2)(state.basis.HOqn.find({0, 1, 8, 0, 0}), state.basis.HOqn.find({0, 2, 0, 0, 0})), MULTIPOLE_OPERATORS_CALCQL0MATRICES_TARGET30, err);
  ASSERT_NEAR(multipoleOperators.ql0(2)(state.basis.HOqn.find({1, 3, 1, 0, 0}), state.basis.HOqn.find({1, 0, 1, 0, 0})), MULTIPOLE_OPERATORS_CALCQL0MATRICES_TARGET31, err);
  ASSERT_NEAR(multipoleOperators.ql0(2)(state.basis.HOqn.find({1, 3, 1, 0, 0}), state.basis.HOqn.find({1, 0, 2, 0, 0})), MULTIPOLE_OPERATORS_CALCQL0MATRICES_TARGET32, err);
  ASSERT_NEAR(multipoleOperators.ql0(2)(state.basis.HOqn.find({1, 3, 1, 0, 0}), state.basis.HOqn.find({1, 3, 2, 0, 0})), MULTIPOLE_OPERATORS_CALCQL0MATRICES_TARGET33, err);
  ASSERT_NEAR(multipoleOperators.ql0(2)(state.basis.HOqn.find({1, 3, 1, 0, 0}), state.basis.HOqn.find({1, 3, 3, 0, 0})), MULTIPOLE_OPERATORS_CALCQL0MATRICES_TARGET34, err);
  ASSERT_NEAR(multipoleOperators.ql0(2)(state.basis.HOqn.find({2, 1, 2, 0, 0}), state.basis.HOqn.find({2, 1, 5, 0, 0})), MULTIPOLE_OPERATORS_CALCQL0MATRICES_TARGET35, err);
  ASSERT_NEAR(multipoleOperators.ql0(2)(state.basis.HOqn.find({2, 1, 2, 0, 0}), state.basis.HOqn.find({2, 1, 6, 0, 0})), MULTIPOLE_OPERATORS_CALCQL0MATRICES_TARGET36, err);

  ASSERT_NEAR(multipoleOperators.ql0(3)(state.basis.HOqn.find({0, 1, 8, 0, 0}), state.basis.HOqn.find({0, 1, 9, 0, 0})), MULTIPOLE_OPERATORS_CALCQL0MATRICES_TARGET37, err);
  ASSERT_NEAR(multipoleOperators.ql0(3)(state.basis.HOqn.find({0, 1, 8, 0, 0}), state.basis.HOqn.find({0, 2, 0, 0, 0})), MULTIPOLE_OPERATORS_CALCQL0MATRICES_TARGET38, err);
  ASSERT_NEAR(multipoleOperators.ql0(3)(state.basis.HOqn.find({1, 3, 1, 0, 0}), state.basis.HOqn.find({1, 0, 1, 0, 0})), MULTIPOLE_OPERATORS_CALCQL0MATRICES_TARGET39, err);
  ASSERT_NEAR(multipoleOperators.ql0(3)(state.basis.HOqn.find({1, 3, 1, 0, 0}), state.basis.HOqn.find({1, 0, 2, 0, 0})), MULTIPOLE_OPERATORS_CALCQL0MATRICES_TARGET40, err);
  ASSERT_NEAR(multipoleOperators.ql0(3)(state.basis.HOqn.find({1, 3, 1, 0, 0}), state.basis.HOqn.find({1, 3, 2, 0, 0})), MULTIPOLE_OPERATORS_CALCQL0MATRICES_TARGET41, err);
  ASSERT_NEAR(multipoleOperators.ql0(3)(state.basis.HOqn.find({1, 3, 1, 0, 0}), state.basis.HOqn.find({1, 3, 3, 0, 0})), MULTIPOLE_OPERATORS_CALCQL0MATRICES_TARGET42, err);
  ASSERT_NEAR(multipoleOperators.ql0(3)(state.basis.HOqn.find({2, 1, 2, 0, 0}), state.basis.HOqn.find({2, 1, 5, 0, 0})), MULTIPOLE_OPERATORS_CALCQL0MATRICES_TARGET43, err);
  ASSERT_NEAR(multipoleOperators.ql0(3)(state.basis.HOqn.find({2, 1, 2, 0, 0}), state.basis.HOqn.find({2, 1, 6, 0, 0})), MULTIPOLE_OPERATORS_CALCQL0MATRICES_TARGET44, err);

  ASSERT_NEAR(multipoleOperators.ql0(4)(state.basis.HOqn.find({0, 1, 8, 0, 0}), state.basis.HOqn.find({0, 1, 9, 0, 0})), MULTIPOLE_OPERATORS_CALCQL0MATRICES_TARGET45, err);
  ASSERT_NEAR(multipoleOperators.ql0(4)(state.basis.HOqn.find({0, 1, 8, 0, 0}), state.basis.HOqn.find({0, 2, 0, 0, 0})), MULTIPOLE_OPERATORS_CALCQL0MATRICES_TARGET46, err);
  ASSERT_NEAR(multipoleOperators.ql0(4)(state.basis.HOqn.find({1, 3, 1, 0, 0}), state.basis.HOqn.find({1, 0, 1, 0, 0})), MULTIPOLE_OPERATORS_CALCQL0MATRICES_TARGET47, err);
  ASSERT_NEAR(multipoleOperators.ql0(4)(state.basis.HOqn.find({1, 3, 1, 0, 0}), state.basis.HOqn.find({1, 0, 2, 0, 0})), MULTIPOLE_OPERATORS_CALCQL0MATRICES_TARGET48, err);
  ASSERT_NEAR(multipoleOperators.ql0(4)(state.basis.HOqn.find({1, 3, 1, 0, 0}), state.basis.HOqn.find({1, 3, 2, 0, 0})), MULTIPOLE_OPERATORS_CALCQL0MATRICES_TARGET49, err);
  ASSERT_NEAR(multipoleOperators.ql0(4)(state.basis.HOqn.find({1, 3, 1, 0, 0}), state.basis.HOqn.find({1, 3, 3, 0, 0})), MULTIPOLE_OPERATORS_CALCQL0MATRICES_TARGET50, err);
  ASSERT_NEAR(multipoleOperators.ql0(4)(state.basis.HOqn.find({2, 1, 2, 0, 0}), state.basis.HOqn.find({2, 1, 5, 0, 0})), MULTIPOLE_OPERATORS_CALCQL0MATRICES_TARGET51, err);
  ASSERT_NEAR(multipoleOperators.ql0(4)(state.basis.HOqn.find({2, 1, 2, 0, 0}), state.basis.HOqn.find({2, 1, 6, 0, 0})), MULTIPOLE_OPERATORS_CALCQL0MATRICES_TARGET52, err);

  ASSERT_NEAR(multipoleOperators.ql0(5)(state.basis.HOqn.find({0, 1, 8, 0, 0}), state.basis.HOqn.find({0, 1, 9, 0, 0})), MULTIPOLE_OPERATORS_CALCQL0MATRICES_TARGET53, err);
  ASSERT_NEAR(multipoleOperators.ql0(5)(state.basis.HOqn.find({0, 1, 8, 0, 0}), state.basis.HOqn.find({0, 2, 0, 0, 0})), MULTIPOLE_OPERATORS_CALCQL0MATRICES_TARGET54, err);
  ASSERT_NEAR(multipoleOperators.ql0(5)(state.basis.HOqn.find({1, 3, 1, 0, 0}), state.basis.HOqn.find({1, 0, 1, 0, 0})), MULTIPOLE_OPERATORS_CALCQL0MATRICES_TARGET55, err);
  ASSERT_NEAR(multipoleOperators.ql0(5)(state.basis.HOqn.find({1, 3, 1, 0, 0}), state.basis.HOqn.find({1, 0, 2, 0, 0})), MULTIPOLE_OPERATORS_CALCQL0MATRICES_TARGET56, err);
  ASSERT_NEAR(multipoleOperators.ql0(5)(state.basis.HOqn.find({1, 3, 1, 0, 0}), state.basis.HOqn.find({1, 3, 2, 0, 0})), MULTIPOLE_OPERATORS_CALCQL0MATRICES_TARGET57, err);
  ASSERT_NEAR(multipoleOperators.ql0(5)(state.basis.HOqn.find({1, 3, 1, 0, 0}), state.basis.HOqn.find({1, 3, 3, 0, 0})), MULTIPOLE_OPERATORS_CALCQL0MATRICES_TARGET58, err);
  ASSERT_NEAR(multipoleOperators.ql0(5)(state.basis.HOqn.find({2, 1, 2, 0, 0}), state.basis.HOqn.find({2, 1, 5, 0, 0})), MULTIPOLE_OPERATORS_CALCQL0MATRICES_TARGET59, err);
  ASSERT_NEAR(multipoleOperators.ql0(5)(state.basis.HOqn.find({2, 1, 2, 0, 0}), state.basis.HOqn.find({2, 1, 6, 0, 0})), MULTIPOLE_OPERATORS_CALCQL0MATRICES_TARGET60, err);

  ASSERT_NEAR(multipoleOperators.ql0(6)(state.basis.HOqn.find({0, 1, 8, 0, 0}), state.basis.HOqn.find({0, 1, 9, 0, 0})), MULTIPOLE_OPERATORS_CALCQL0MATRICES_TARGET61, err);
  ASSERT_NEAR(multipoleOperators.ql0(6)(state.basis.HOqn.find({0, 1, 8, 0, 0}), state.basis.HOqn.find({0, 2, 0, 0, 0})), MULTIPOLE_OPERATORS_CALCQL0MATRICES_TARGET62, err);
  ASSERT_NEAR(multipoleOperators.ql0(6)(state.basis.HOqn.find({1, 3, 1, 0, 0}), state.basis.HOqn.find({1, 0, 1, 0, 0})), MULTIPOLE_OPERATORS_CALCQL0MATRICES_TARGET63, err);
  ASSERT_NEAR(multipoleOperators.ql0(6)(state.basis.HOqn.find({1, 3, 1, 0, 0}), state.basis.HOqn.find({1, 0, 2, 0, 0})), MULTIPOLE_OPERATORS_CALCQL0MATRICES_TARGET64, err);
  ASSERT_NEAR(multipoleOperators.ql0(6)(state.basis.HOqn.find({1, 3, 1, 0, 0}), state.basis.HOqn.find({1, 3, 2, 0, 0})), MULTIPOLE_OPERATORS_CALCQL0MATRICES_TARGET65, err);
  ASSERT_NEAR(multipoleOperators.ql0(6)(state.basis.HOqn.find({1, 3, 1, 0, 0}), state.basis.HOqn.find({1, 3, 3, 0, 0})), MULTIPOLE_OPERATORS_CALCQL0MATRICES_TARGET66, err);
  ASSERT_NEAR(multipoleOperators.ql0(6)(state.basis.HOqn.find({2, 1, 2, 0, 0}), state.basis.HOqn.find({2, 1, 5, 0, 0})), MULTIPOLE_OPERATORS_CALCQL0MATRICES_TARGET67, err);
  ASSERT_NEAR(multipoleOperators.ql0(6)(state.basis.HOqn.find({2, 1, 2, 0, 0}), state.basis.HOqn.find({2, 1, 6, 0, 0})), MULTIPOLE_OPERATORS_CALCQL0MATRICES_TARGET68, err);
}

//==============================================================================

TEST(MultipoleOperators, MultipoleOperators)
{
  {
    // 1ct state
    State state("examples/42Ca_deformed_1x11.msg.gz");
    MultipoleOperators multipoleOperators(state);

    double err = 1e-08;

    // ALL_OUT;
    // INFO("%14.10f", multipoleOperators.qlm(0, 0));
    // INFO("%14.10f", multipoleOperators.qlm(1, 0));
    // INFO("%14.10f", multipoleOperators.qlm(2, 0));
    // INFO("%14.10f", multipoleOperators.qlm(3, 0));
    // INFO("%14.10f", multipoleOperators.qlm(4, 0));
    // INFO("%14.10f", multipoleOperators.qlm(0, 1));
    // INFO("%14.10f", multipoleOperators.qlm(1, 1));
    // INFO("%14.10f", multipoleOperators.qlm(2, 1));
    // INFO("%14.10f", multipoleOperators.qlm(3, 1));
    // INFO("%14.10f", multipoleOperators.qlm(4, 1));
    // INFO("%14.10f", multipoleOperators.qlm(0, 2));
    // INFO("%14.10f", multipoleOperators.qlm(1, 2));
    // INFO("%14.10f", multipoleOperators.qlm(2, 2));
    // INFO("%14.10f", multipoleOperators.qlm(3, 2));
    // INFO("%14.10f", multipoleOperators.qlm(4, 2));

    ASSERT_NEAR(multipoleOperators.qlm(0, 0),   6.2060854344, err);
    ASSERT_NEAR(multipoleOperators.qlm(1, 0),  -0.0000000086, err);
    ASSERT_NEAR(multipoleOperators.qlm(2, 0),   6.3796768665, err);
    ASSERT_NEAR(multipoleOperators.qlm(3, 0),   0.0000016739, err);
    ASSERT_NEAR(multipoleOperators.qlm(4, 0), 137.3177593081, err);
    ASSERT_NEAR(multipoleOperators.qlm(0, 1),   5.6418958363, err);
    ASSERT_NEAR(multipoleOperators.qlm(1, 1),   0.0000000086, err);
    ASSERT_NEAR(multipoleOperators.qlm(2, 1),   3.6203231335, err);
    ASSERT_NEAR(multipoleOperators.qlm(3, 1),   0.0000019036, err);
    ASSERT_NEAR(multipoleOperators.qlm(4, 1),  40.7764698021, err);
    ASSERT_NEAR(multipoleOperators.qlm(0, 2),  11.8479812707, err);
    ASSERT_NEAR(multipoleOperators.qlm(1, 2),  -0.0000000000, err);
    ASSERT_NEAR(multipoleOperators.qlm(2, 2),  10.0000000000, err);
    ASSERT_NEAR(multipoleOperators.qlm(3, 2),   0.0000035775, err);
    ASSERT_NEAR(multipoleOperators.qlm(4, 2), 178.0942291102, err);
  }

  {
    // 2ct state
    State state("examples/42Ca_deformed_2x9.msg.gz");
    MultipoleOperators multipoleOperators(state);

    double err = 1e-08;

    // ALL_OUT;
    // INFO("%14.10f", multipoleOperators.qlm(0, 0));
    // INFO("%14.10f", multipoleOperators.qlm(1, 0));
    // INFO("%14.10f", multipoleOperators.qlm(2, 0));
    // INFO("%14.10f", multipoleOperators.qlm(3, 0));
    // INFO("%14.10f", multipoleOperators.qlm(4, 0));
    // INFO("%14.10f", multipoleOperators.qlm(0, 1));
    // INFO("%14.10f", multipoleOperators.qlm(1, 1));
    // INFO("%14.10f", multipoleOperators.qlm(2, 1));
    // INFO("%14.10f", multipoleOperators.qlm(3, 1));
    // INFO("%14.10f", multipoleOperators.qlm(4, 1));
    // INFO("%14.10f", multipoleOperators.qlm(0, 2));
    // INFO("%14.10f", multipoleOperators.qlm(1, 2));
    // INFO("%14.10f", multipoleOperators.qlm(2, 2));
    // INFO("%14.10f", multipoleOperators.qlm(3, 2));
    // INFO("%14.10f", multipoleOperators.qlm(4, 2));

    ASSERT_NEAR(multipoleOperators.qlm(0, 0),  6.2060836799, err);
    ASSERT_NEAR(multipoleOperators.qlm(1, 0), -0.0000001621, err);
    ASSERT_NEAR(multipoleOperators.qlm(2, 0),  6.3539566216, err);
    ASSERT_NEAR(multipoleOperators.qlm(3, 0), -0.0001746998, err);
    ASSERT_NEAR(multipoleOperators.qlm(4, 0), 64.2799825149, err);
    ASSERT_NEAR(multipoleOperators.qlm(0, 1),  5.6418958352, err);
    ASSERT_NEAR(multipoleOperators.qlm(1, 1),  0.0000001621, err);
    ASSERT_NEAR(multipoleOperators.qlm(2, 1),  3.6460433784, err);
    ASSERT_NEAR(multipoleOperators.qlm(3, 1), -0.0001754108, err);
    ASSERT_NEAR(multipoleOperators.qlm(4, 1), 15.9107568066, err);
    ASSERT_NEAR(multipoleOperators.qlm(0, 2), 11.8479795151, err);
    ASSERT_NEAR(multipoleOperators.qlm(1, 2), -0.0000000000, err);
    ASSERT_NEAR(multipoleOperators.qlm(2, 2), 10.0000000000, err);
    ASSERT_NEAR(multipoleOperators.qlm(3, 2), -0.0003501106, err);
    ASSERT_NEAR(multipoleOperators.qlm(4, 2), 80.1907393215, err);
  }
}

//==============================================================================

/*
TEST(MultipoleOperators, printQlm)
{
}
*/

//==============================================================================

/*
TEST(MultipoleOperators, info)
{
}
*/



