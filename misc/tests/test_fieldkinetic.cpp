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
 *  \brief Test suite for the FieldKinetic class.
 */

#include "gtest/gtest.h"
#include "interaction.h"
#include "field_kinetic.h"
#include "state.h"

//==============================================================================

TEST(FieldKinetic, calcField)
{
  State state("examples/42Ca_deformed_2x9.msg.gz");

  FieldKinetic field(Field::Parameters(), &state);
  field.calcField();

  // Target values calculated with Mathematica (cf. misc/mathematica/test_kinetic.wls)
  ASSERT_NEAR(field.field(NEUTRON, Field::DIRECT)(196, 142), 0.0000000000000, 1e-9);
  ASSERT_NEAR(field.field(PROTON , Field::DIRECT)(182, 192), 5.8185479366527, 1e-9);
}
