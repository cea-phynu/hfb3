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
 *  \brief Meta test suite (tests not directly related to a method).
 */

#include "gtest/gtest.h"
#include "interaction.h"
#include "state.h"
#include "interaction.h"

#include "tools.h"

//==============================================================================

TEST(Meta, energy_D1S_1ct)
{
  double err = 1e-6;
  DataTree dataTree("examples/42Ca_deformed_1x11.msg.gz");
  State state(dataTree);

  Interaction interaction("D1S", &state);
  interaction.calcEnergies();

  // ALL_OUT;
  // INFO("%14.6f", interaction("kinetic"         , 0)->energy(NEUTRON, Field::DIRECT  ));
  // INFO("%14.6f", interaction("coulomb (Slater)", 0)->energy(NEUTRON, Field::DIRECT  ));
  // INFO("%14.6f", interaction("central"         , 0)->energy(NEUTRON, Field::DIRECT  ));
  // INFO("%14.6f", interaction("central"         , 1)->energy(NEUTRON, Field::DIRECT  ));
  // INFO("%14.6f", interaction("central"         , 0)->energy(NEUTRON, Field::EXCHANGE));
  // INFO("%14.6f", interaction("central"         , 1)->energy(NEUTRON, Field::EXCHANGE));
  // INFO("%14.6f", interaction("spin-orbit"      , 0)->energy(NEUTRON, Field::DIRECT  ));
  // INFO("%14.6f", interaction("density"         , 0)->energy(NEUTRON, Field::DIRECT  ));
  // INFO("%14.6f", interaction("2-body COM cor." , 0)->energy(NEUTRON, Field::DIRECT  ));
  // INFO("%14.6f", interaction("2-body COM cor." , 0)->energy(NEUTRON, Field::PAIRING ));
  // INFO("%14.6f", interaction("coulomb (Slater)", 0)->energy(NEUTRON, Field::EXCHANGE));
  // INFO("%14.6f", interaction("central"         , 0)->energy(NEUTRON, Field::PAIRING ));
  // INFO("%14.6f", interaction("central"         , 1)->energy(NEUTRON, Field::PAIRING ));
  // INFO("%14.6f", interaction("rearrangement"   , 0)->energy(NEUTRON, Field::DIRECT  ));
  // INFO("%14.6f", interaction("kinetic"         , 0)->energy(PROTON , Field::DIRECT  ));
  // INFO("%14.6f", interaction("coulomb (Slater)", 0)->energy(PROTON , Field::DIRECT  ));
  // INFO("%14.6f", interaction("central"         , 0)->energy(PROTON , Field::DIRECT  ));
  // INFO("%14.6f", interaction("central"         , 1)->energy(PROTON , Field::DIRECT  ));
  // INFO("%14.6f", interaction("central"         , 0)->energy(PROTON , Field::EXCHANGE));
  // INFO("%14.6f", interaction("central"         , 1)->energy(PROTON , Field::EXCHANGE));
  // INFO("%14.6f", interaction("spin-orbit"      , 0)->energy(PROTON , Field::DIRECT  ));
  // INFO("%14.6f", interaction("density"         , 0)->energy(PROTON , Field::DIRECT  ));
  // INFO("%14.6f", interaction("2-body COM cor." , 0)->energy(PROTON , Field::DIRECT  ));
  // INFO("%14.6f", interaction("2-body COM cor." , 0)->energy(PROTON , Field::PAIRING ));
  // INFO("%14.6f", interaction("coulomb (Slater)", 0)->energy(PROTON , Field::EXCHANGE));
  // INFO("%14.6f", interaction("central"         , 0)->energy(PROTON , Field::PAIRING ));
  // INFO("%14.6f", interaction("central"         , 1)->energy(PROTON , Field::PAIRING ));
  // INFO("%14.6f", interaction("rearrangement"   , 0)->energy(PROTON , Field::DIRECT  ));

  ASSERT_NEAR(interaction("kinetic"         , 0)->energy(NEUTRON, Field::DIRECT  ),   369.860559, err);
  ASSERT_NEAR(interaction("coulomb (Slater)", 0)->energy(NEUTRON, Field::DIRECT  ),     0.000000, err);
  ASSERT_NEAR(interaction("central"         , 0)->energy(NEUTRON, Field::DIRECT  ), -1001.729780, err);
  ASSERT_NEAR(interaction("central"         , 1)->energy(NEUTRON, Field::DIRECT  ),   -44.636063, err);
  ASSERT_NEAR(interaction("central"         , 0)->energy(NEUTRON, Field::EXCHANGE),   456.097494, err);
  ASSERT_NEAR(interaction("central"         , 1)->energy(NEUTRON, Field::EXCHANGE),  -527.734175, err);
  ASSERT_NEAR(interaction("spin-orbit"      , 0)->energy(NEUTRON, Field::DIRECT  ),    -7.100354, err);
  ASSERT_NEAR(interaction("density"         , 0)->energy(NEUTRON, Field::DIRECT  ),   552.709234, err);
  ASSERT_NEAR(interaction("2-body COM cor." , 0)->energy(NEUTRON, Field::DIRECT  ),     4.621814, err);
  ASSERT_NEAR(interaction("2-body COM cor." , 0)->energy(NEUTRON, Field::PAIRING ),     0.205243, err);
  ASSERT_NEAR(interaction("coulomb (Slater)", 0)->energy(NEUTRON, Field::EXCHANGE),     0.000000, err);
  ASSERT_NEAR(interaction("central"         , 0)->energy(NEUTRON, Field::PAIRING ),     4.492858, err);
  ASSERT_NEAR(interaction("central"         , 1)->energy(NEUTRON, Field::PAIRING ),    -8.902637, err);
  ASSERT_NEAR(interaction("rearrangement"   , 0)->energy(NEUTRON, Field::DIRECT  ),    96.139862, err);
  ASSERT_NEAR(interaction("kinetic"         , 0)->energy(PROTON , Field::DIRECT  ),   309.733474, err);
  ASSERT_NEAR(interaction("coulomb (Slater)", 0)->energy(PROTON , Field::DIRECT  ),    79.044427, err);
  ASSERT_NEAR(interaction("central"         , 0)->energy(PROTON , Field::DIRECT  ), -1009.472105, err);
  ASSERT_NEAR(interaction("central"         , 1)->energy(PROTON , Field::DIRECT  ),   -21.056042, err);
  ASSERT_NEAR(interaction("central"         , 0)->energy(PROTON , Field::EXCHANGE),   449.152348, err);
  ASSERT_NEAR(interaction("central"         , 1)->energy(PROTON , Field::EXCHANGE),  -514.596069, err);
  ASSERT_NEAR(interaction("spin-orbit"      , 0)->energy(PROTON , Field::DIRECT  ),    -1.675028, err);
  ASSERT_NEAR(interaction("density"         , 0)->energy(PROTON , Field::DIRECT  ),   552.709234, err);
  ASSERT_NEAR(interaction("2-body COM cor." , 0)->energy(PROTON , Field::DIRECT  ),     3.673729, err);
  ASSERT_NEAR(interaction("2-body COM cor." , 0)->energy(PROTON , Field::PAIRING ),     0.000000, err);
  ASSERT_NEAR(interaction("coulomb (Slater)", 0)->energy(PROTON , Field::EXCHANGE),    -7.440938, err);
  ASSERT_NEAR(interaction("central"         , 0)->energy(PROTON , Field::PAIRING ),     0.000000, err);
  ASSERT_NEAR(interaction("central"         , 1)->energy(PROTON , Field::PAIRING ),    -0.000000, err);
  ASSERT_NEAR(interaction("rearrangement"   , 0)->energy(PROTON , Field::DIRECT  ),    88.095730, err);
}

//==============================================================================

TEST(Meta, energy_D1S_2ct)
{
  double err = 1e-6;
  DataTree dataTree("examples/42Ca_deformed_2x9.msg.gz");
  State state(dataTree);

  Interaction interaction("D1S", &state);
  interaction.calcEnergies();

  // ALL_OUT;
  // INFO("%14.6f", interaction("kinetic"         , 0)->energy(NEUTRON, Field::DIRECT  ));
  // INFO("%14.6f", interaction("coulomb (Slater)", 0)->energy(NEUTRON, Field::DIRECT  ));
  // INFO("%14.6f", interaction("central"         , 0)->energy(NEUTRON, Field::DIRECT  ));
  // INFO("%14.6f", interaction("central"         , 1)->energy(NEUTRON, Field::DIRECT  ));
  // INFO("%14.6f", interaction("central"         , 0)->energy(NEUTRON, Field::EXCHANGE));
  // INFO("%14.6f", interaction("central"         , 1)->energy(NEUTRON, Field::EXCHANGE));
  // INFO("%14.6f", interaction("spin-orbit"      , 0)->energy(NEUTRON, Field::DIRECT  ));
  // INFO("%14.6f", interaction("density"         , 0)->energy(NEUTRON, Field::DIRECT  ));
  // INFO("%14.6f", interaction("2-body COM cor." , 0)->energy(NEUTRON, Field::DIRECT  ));
  // INFO("%14.6f", interaction("2-body COM cor." , 0)->energy(NEUTRON, Field::PAIRING ));
  // INFO("%14.6f", interaction("coulomb (Slater)", 0)->energy(NEUTRON, Field::EXCHANGE));
  // INFO("%14.6f", interaction("central"         , 0)->energy(NEUTRON, Field::PAIRING ));
  // INFO("%14.6f", interaction("central"         , 1)->energy(NEUTRON, Field::PAIRING ));
  // INFO("%14.6f", interaction("rearrangement"   , 0)->energy(NEUTRON, Field::DIRECT  ));
  // INFO("%14.6f", interaction("kinetic"         , 0)->energy(PROTON , Field::DIRECT  ));
  // INFO("%14.6f", interaction("coulomb (Slater)", 0)->energy(PROTON , Field::DIRECT  ));
  // INFO("%14.6f", interaction("central"         , 0)->energy(PROTON , Field::DIRECT  ));
  // INFO("%14.6f", interaction("central"         , 1)->energy(PROTON , Field::DIRECT  ));
  // INFO("%14.6f", interaction("central"         , 0)->energy(PROTON , Field::EXCHANGE));
  // INFO("%14.6f", interaction("central"         , 1)->energy(PROTON , Field::EXCHANGE));
  // INFO("%14.6f", interaction("spin-orbit"      , 0)->energy(PROTON , Field::DIRECT  ));
  // INFO("%14.6f", interaction("density"         , 0)->energy(PROTON , Field::DIRECT  ));
  // INFO("%14.6f", interaction("2-body COM cor." , 0)->energy(PROTON , Field::DIRECT  ));
  // INFO("%14.6f", interaction("2-body COM cor." , 0)->energy(PROTON , Field::PAIRING ));
  // INFO("%14.6f", interaction("coulomb (Slater)", 0)->energy(PROTON , Field::EXCHANGE));
  // INFO("%14.6f", interaction("central"         , 0)->energy(PROTON , Field::PAIRING ));
  // INFO("%14.6f", interaction("central"         , 1)->energy(PROTON , Field::PAIRING ));
  // INFO("%14.6f", interaction("rearrangement"   , 0)->energy(PROTON , Field::DIRECT  ));

  ASSERT_NEAR(interaction("kinetic"         , 0)->energy(NEUTRON, Field::DIRECT  ),   370.706378, err);
  ASSERT_NEAR(interaction("coulomb (Slater)", 0)->energy(NEUTRON, Field::DIRECT  ),     0.000000, err);
  ASSERT_NEAR(interaction("central"         , 0)->energy(NEUTRON, Field::DIRECT  ), -1007.124174, err);
  ASSERT_NEAR(interaction("central"         , 1)->energy(NEUTRON, Field::DIRECT  ),   -44.992954, err);
  ASSERT_NEAR(interaction("central"         , 0)->energy(NEUTRON, Field::EXCHANGE),   458.164338, err);
  ASSERT_NEAR(interaction("central"         , 1)->energy(NEUTRON, Field::EXCHANGE),  -529.498430, err);
  ASSERT_NEAR(interaction("spin-orbit"      , 0)->energy(NEUTRON, Field::DIRECT  ),    -7.154445, err);
  ASSERT_NEAR(interaction("density"         , 0)->energy(NEUTRON, Field::DIRECT  ),   557.056358, err);
  ASSERT_NEAR(interaction("2-body COM cor." , 0)->energy(NEUTRON, Field::DIRECT  ),     4.629348, err);
  ASSERT_NEAR(interaction("2-body COM cor." , 0)->energy(NEUTRON, Field::PAIRING ),     0.201743, err);
  ASSERT_NEAR(interaction("coulomb (Slater)", 0)->energy(NEUTRON, Field::EXCHANGE),     0.000000, err);
  ASSERT_NEAR(interaction("central"         , 0)->energy(NEUTRON, Field::PAIRING ),     4.329008, err);
  ASSERT_NEAR(interaction("central"         , 1)->energy(NEUTRON, Field::PAIRING ),    -8.604896, err);
  ASSERT_NEAR(interaction("rearrangement"   , 0)->energy(NEUTRON, Field::DIRECT  ),    96.940644, err);
  ASSERT_NEAR(interaction("kinetic"         , 0)->energy(PROTON , Field::DIRECT  ),   310.321578, err);
  ASSERT_NEAR(interaction("coulomb (Slater)", 0)->energy(PROTON , Field::DIRECT  ),    79.170771, err);
  ASSERT_NEAR(interaction("central"         , 0)->energy(PROTON , Field::DIRECT  ), -1014.983218, err);
  ASSERT_NEAR(interaction("central"         , 1)->energy(PROTON , Field::DIRECT  ),   -21.069868, err);
  ASSERT_NEAR(interaction("central"         , 0)->energy(PROTON , Field::EXCHANGE),   451.122647, err);
  ASSERT_NEAR(interaction("central"         , 1)->energy(PROTON , Field::EXCHANGE),  -516.202991, err);
  ASSERT_NEAR(interaction("spin-orbit"      , 0)->energy(PROTON , Field::DIRECT  ),    -1.685199, err);
  ASSERT_NEAR(interaction("density"         , 0)->energy(PROTON , Field::DIRECT  ),   557.056360, err);
  ASSERT_NEAR(interaction("2-body COM cor." , 0)->energy(PROTON , Field::DIRECT  ),     3.678276, err);
  ASSERT_NEAR(interaction("2-body COM cor." , 0)->energy(PROTON , Field::PAIRING ),     0.000000, err);
  ASSERT_NEAR(interaction("coulomb (Slater)", 0)->energy(PROTON , Field::EXCHANGE),    -7.449355, err);
  ASSERT_NEAR(interaction("central"         , 0)->energy(PROTON , Field::PAIRING ),     0.000000, err);
  ASSERT_NEAR(interaction("central"         , 1)->energy(PROTON , Field::PAIRING ),    -0.000000, err);
  ASSERT_NEAR(interaction("rearrangement"   , 0)->energy(PROTON , Field::DIRECT  ),    88.744058, err);
}

//==============================================================================

TEST(Meta, energy_D12S_1ct)
{
  double err = 1e-6;
  DataTree dataTree("examples/42Ca_deformed_1x11.msg.gz");
  State state(dataTree);

  Interaction interaction("D12S", &state);
  interaction.calcEnergies();

  // ALL_OUT;
  // INFO("%14.6f", interaction("kinetic"         , 0)->energy(NEUTRON, Field::DIRECT  ));
  // INFO("%14.6f", interaction("coulomb (Slater)", 0)->energy(NEUTRON, Field::DIRECT  ));
  // INFO("%14.6f", interaction("central"         , 0)->energy(NEUTRON, Field::DIRECT  ));
  // INFO("%14.6f", interaction("central"         , 1)->energy(NEUTRON, Field::DIRECT  ));
  // INFO("%14.6f", interaction("central"         , 0)->energy(NEUTRON, Field::EXCHANGE));
  // INFO("%14.6f", interaction("central"         , 1)->energy(NEUTRON, Field::EXCHANGE));
  // INFO("%14.6f", interaction("spin-orbit"      , 0)->energy(NEUTRON, Field::DIRECT  ));
  // INFO("%14.6f", interaction("2-body COM cor." , 0)->energy(NEUTRON, Field::DIRECT  ));
  // INFO("%14.6f", interaction("2-body COM cor." , 0)->energy(NEUTRON, Field::PAIRING ));
  // INFO("%14.6f", interaction("coulomb (Slater)", 0)->energy(NEUTRON, Field::EXCHANGE));
  // INFO("%14.6f", interaction("central"         , 0)->energy(NEUTRON, Field::PAIRING ));
  // INFO("%14.6f", interaction("central"         , 1)->energy(NEUTRON, Field::PAIRING ));
  // INFO("%14.6f", interaction("density FR"      , 0)->energy(NEUTRON, Field::DIRECT  ));
  // INFO("%14.6f", interaction("density FR"      , 0)->energy(NEUTRON, Field::EXCHANGE));
  // INFO("%14.6f", interaction("density FR"      , 0)->energy(NEUTRON, Field::PAIRING ));
  // INFO("%14.6f", interaction("rearrangement FR", 0)->energy(NEUTRON, Field::DIRECT  ));
  // INFO("%14.6f", interaction("rearrangement FR", 0)->energy(NEUTRON, Field::EXCHANGE));
  // INFO("%14.6f", interaction("rearrangement FR", 0)->energy(NEUTRON, Field::PAIRING ));
  // INFO("%14.6f", interaction("kinetic"         , 0)->energy(PROTON , Field::DIRECT  ));
  // INFO("%14.6f", interaction("coulomb (Slater)", 0)->energy(PROTON , Field::DIRECT  ));
  // INFO("%14.6f", interaction("central"         , 0)->energy(PROTON , Field::DIRECT  ));
  // INFO("%14.6f", interaction("central"         , 1)->energy(PROTON , Field::DIRECT  ));
  // INFO("%14.6f", interaction("central"         , 0)->energy(PROTON , Field::EXCHANGE));
  // INFO("%14.6f", interaction("central"         , 1)->energy(PROTON , Field::EXCHANGE));
  // INFO("%14.6f", interaction("spin-orbit"      , 0)->energy(PROTON , Field::DIRECT  ));
  // INFO("%14.6f", interaction("2-body COM cor." , 0)->energy(PROTON , Field::DIRECT  ));
  // INFO("%14.6f", interaction("2-body COM cor." , 0)->energy(PROTON , Field::PAIRING ));
  // INFO("%14.6f", interaction("coulomb (Slater)", 0)->energy(PROTON , Field::EXCHANGE));
  // INFO("%14.6f", interaction("central"         , 0)->energy(PROTON , Field::PAIRING ));
  // INFO("%14.6f", interaction("central"         , 1)->energy(PROTON , Field::PAIRING ));
  // INFO("%14.6f", interaction("density FR"      , 0)->energy(PROTON , Field::DIRECT  ));
  // INFO("%14.6f", interaction("density FR"      , 0)->energy(PROTON , Field::EXCHANGE));
  // INFO("%14.6f", interaction("density FR"      , 0)->energy(PROTON , Field::PAIRING ));
  // INFO("%14.6f", interaction("rearrangement FR", 0)->energy(PROTON , Field::DIRECT  ));
  // INFO("%14.6f", interaction("rearrangement FR", 0)->energy(PROTON , Field::EXCHANGE));
  // INFO("%14.6f", interaction("rearrangement FR", 0)->energy(PROTON , Field::PAIRING ));

  ASSERT_NEAR(interaction("kinetic"         , 0)->energy(NEUTRON, Field::DIRECT  ),   369.860559, err);
  ASSERT_NEAR(interaction("coulomb (Slater)", 0)->energy(NEUTRON, Field::DIRECT  ),     0.000000, err);
  ASSERT_NEAR(interaction("central"         , 0)->energy(NEUTRON, Field::DIRECT  ), -1001.729780, err);
  ASSERT_NEAR(interaction("central"         , 1)->energy(NEUTRON, Field::DIRECT  ),   -44.636063, err);
  ASSERT_NEAR(interaction("central"         , 0)->energy(NEUTRON, Field::EXCHANGE),   456.097494, err);
  ASSERT_NEAR(interaction("central"         , 1)->energy(NEUTRON, Field::EXCHANGE),  -527.734175, err);
  ASSERT_NEAR(interaction("spin-orbit"      , 0)->energy(NEUTRON, Field::DIRECT  ),    -7.100354, err);
  ASSERT_NEAR(interaction("2-body COM cor." , 0)->energy(NEUTRON, Field::DIRECT  ),     4.621814, err);
  ASSERT_NEAR(interaction("2-body COM cor." , 0)->energy(NEUTRON, Field::PAIRING ),     0.205243, err);
  ASSERT_NEAR(interaction("coulomb (Slater)", 0)->energy(NEUTRON, Field::EXCHANGE),     0.000000, err);
  ASSERT_NEAR(interaction("central"         , 0)->energy(NEUTRON, Field::PAIRING ),     4.492858, err);
  ASSERT_NEAR(interaction("central"         , 1)->energy(NEUTRON, Field::PAIRING ),    -8.902637, err);
  ASSERT_NEAR(interaction("density FR"      , 0)->energy(NEUTRON, Field::DIRECT  ),  1155.987374, err);
  ASSERT_NEAR(interaction("density FR"      , 0)->energy(NEUTRON, Field::EXCHANGE),  -603.282136, err);
  ASSERT_NEAR(interaction("density FR"      , 0)->energy(NEUTRON, Field::PAIRING ),     0.000000, err);
  ASSERT_NEAR(interaction("rearrangement FR", 0)->energy(NEUTRON, Field::DIRECT  ),   192.664562, err);
  ASSERT_NEAR(interaction("rearrangement FR", 0)->energy(NEUTRON, Field::EXCHANGE),   -96.524952, err);
  ASSERT_NEAR(interaction("rearrangement FR", 0)->energy(NEUTRON, Field::PAIRING ),     0.000000, err);
  ASSERT_NEAR(interaction("kinetic"         , 0)->energy(PROTON , Field::DIRECT  ),   309.733474, err);
  ASSERT_NEAR(interaction("coulomb (Slater)", 0)->energy(PROTON , Field::DIRECT  ),    79.044427, err);
  ASSERT_NEAR(interaction("central"         , 0)->energy(PROTON , Field::DIRECT  ), -1009.472105, err);
  ASSERT_NEAR(interaction("central"         , 1)->energy(PROTON , Field::DIRECT  ),   -21.056042, err);
  ASSERT_NEAR(interaction("central"         , 0)->energy(PROTON , Field::EXCHANGE),   449.152348, err);
  ASSERT_NEAR(interaction("central"         , 1)->energy(PROTON , Field::EXCHANGE),  -514.596069, err);
  ASSERT_NEAR(interaction("spin-orbit"      , 0)->energy(PROTON , Field::DIRECT  ),    -1.675028, err);
  ASSERT_NEAR(interaction("2-body COM cor." , 0)->energy(PROTON , Field::DIRECT  ),     3.673729, err);
  ASSERT_NEAR(interaction("2-body COM cor." , 0)->energy(PROTON , Field::PAIRING ),     0.000000, err);
  ASSERT_NEAR(interaction("coulomb (Slater)", 0)->energy(PROTON , Field::EXCHANGE),    -7.440938, err);
  ASSERT_NEAR(interaction("central"         , 0)->energy(PROTON , Field::PAIRING ),     0.000000, err);
  ASSERT_NEAR(interaction("central"         , 1)->energy(PROTON , Field::PAIRING ),    -0.000000, err);
  ASSERT_NEAR(interaction("density FR"      , 0)->energy(PROTON , Field::DIRECT  ),  1059.246056, err);
  ASSERT_NEAR(interaction("density FR"      , 0)->energy(PROTON , Field::EXCHANGE),  -506.540819, err);
  ASSERT_NEAR(interaction("density FR"      , 0)->energy(PROTON , Field::PAIRING ),     0.000000, err);
  ASSERT_NEAR(interaction("rearrangement FR", 0)->energy(PROTON , Field::DIRECT  ),   176.541009, err);
  ASSERT_NEAR(interaction("rearrangement FR", 0)->energy(PROTON , Field::EXCHANGE),   -88.445540, err);
  ASSERT_NEAR(interaction("rearrangement FR", 0)->energy(PROTON , Field::PAIRING ),     0.000000, err);
}

//==============================================================================

TEST(Meta, energy_D12S_2ct)
{
  double err = 1e-6;
  DataTree dataTree("examples/42Ca_deformed_2x9.msg.gz");
  State state(dataTree);

  Interaction interaction("D12S", &state);
  interaction.calcEnergies();

  // ALL_OUT;
  // INFO("%14.6f", interaction("kinetic"         , 0)->energy(NEUTRON, Field::DIRECT  ));
  // INFO("%14.6f", interaction("coulomb (Slater)", 0)->energy(NEUTRON, Field::DIRECT  ));
  // INFO("%14.6f", interaction("central"         , 0)->energy(NEUTRON, Field::DIRECT  ));
  // INFO("%14.6f", interaction("central"         , 1)->energy(NEUTRON, Field::DIRECT  ));
  // INFO("%14.6f", interaction("central"         , 0)->energy(NEUTRON, Field::EXCHANGE));
  // INFO("%14.6f", interaction("central"         , 1)->energy(NEUTRON, Field::EXCHANGE));
  // INFO("%14.6f", interaction("spin-orbit"      , 0)->energy(NEUTRON, Field::DIRECT  ));
  // INFO("%14.6f", interaction("2-body COM cor." , 0)->energy(NEUTRON, Field::DIRECT  ));
  // INFO("%14.6f", interaction("2-body COM cor." , 0)->energy(NEUTRON, Field::PAIRING ));
  // INFO("%14.6f", interaction("coulomb (Slater)", 0)->energy(NEUTRON, Field::EXCHANGE));
  // INFO("%14.6f", interaction("central"         , 0)->energy(NEUTRON, Field::PAIRING ));
  // INFO("%14.6f", interaction("central"         , 1)->energy(NEUTRON, Field::PAIRING ));
  // INFO("%14.6f", interaction("density FR"      , 0)->energy(NEUTRON, Field::DIRECT  ));
  // INFO("%14.6f", interaction("density FR"      , 0)->energy(NEUTRON, Field::EXCHANGE));
  // INFO("%14.6f", interaction("density FR"      , 0)->energy(NEUTRON, Field::PAIRING ));
  // INFO("%14.6f", interaction("rearrangement FR", 0)->energy(NEUTRON, Field::DIRECT  ));
  // INFO("%14.6f", interaction("rearrangement FR", 0)->energy(NEUTRON, Field::EXCHANGE));
  // INFO("%14.6f", interaction("rearrangement FR", 0)->energy(NEUTRON, Field::PAIRING ));
  // INFO("%14.6f", interaction("kinetic"         , 0)->energy(PROTON , Field::DIRECT  ));
  // INFO("%14.6f", interaction("coulomb (Slater)", 0)->energy(PROTON , Field::DIRECT  ));
  // INFO("%14.6f", interaction("central"         , 0)->energy(PROTON , Field::DIRECT  ));
  // INFO("%14.6f", interaction("central"         , 1)->energy(PROTON , Field::DIRECT  ));
  // INFO("%14.6f", interaction("central"         , 0)->energy(PROTON , Field::EXCHANGE));
  // INFO("%14.6f", interaction("central"         , 1)->energy(PROTON , Field::EXCHANGE));
  // INFO("%14.6f", interaction("spin-orbit"      , 0)->energy(PROTON , Field::DIRECT  ));
  // INFO("%14.6f", interaction("2-body COM cor." , 0)->energy(PROTON , Field::DIRECT  ));
  // INFO("%14.6f", interaction("2-body COM cor." , 0)->energy(PROTON , Field::PAIRING ));
  // INFO("%14.6f", interaction("coulomb (Slater)", 0)->energy(PROTON , Field::EXCHANGE));
  // INFO("%14.6f", interaction("central"         , 0)->energy(PROTON , Field::PAIRING ));
  // INFO("%14.6f", interaction("central"         , 1)->energy(PROTON , Field::PAIRING ));
  // INFO("%14.6f", interaction("density FR"      , 0)->energy(PROTON , Field::DIRECT  ));
  // INFO("%14.6f", interaction("density FR"      , 0)->energy(PROTON , Field::EXCHANGE));
  // INFO("%14.6f", interaction("density FR"      , 0)->energy(PROTON , Field::PAIRING ));
  // INFO("%14.6f", interaction("rearrangement FR", 0)->energy(PROTON , Field::DIRECT  ));
  // INFO("%14.6f", interaction("rearrangement FR", 0)->energy(PROTON , Field::EXCHANGE));
  // INFO("%14.6f", interaction("rearrangement FR", 0)->energy(PROTON , Field::PAIRING ));

  ASSERT_NEAR(interaction("kinetic"         , 0)->energy(NEUTRON, Field::DIRECT  ),   370.706378, err);
  ASSERT_NEAR(interaction("coulomb (Slater)", 0)->energy(NEUTRON, Field::DIRECT  ),     0.000000, err);
  ASSERT_NEAR(interaction("central"         , 0)->energy(NEUTRON, Field::DIRECT  ), -1007.124174, err);
  ASSERT_NEAR(interaction("central"         , 1)->energy(NEUTRON, Field::DIRECT  ),   -44.992954, err);
  ASSERT_NEAR(interaction("central"         , 0)->energy(NEUTRON, Field::EXCHANGE),   458.164338, err);
  ASSERT_NEAR(interaction("central"         , 1)->energy(NEUTRON, Field::EXCHANGE),  -529.498430, err);
  ASSERT_NEAR(interaction("spin-orbit"      , 0)->energy(NEUTRON, Field::DIRECT  ),    -7.154445, err);
  ASSERT_NEAR(interaction("2-body COM cor." , 0)->energy(NEUTRON, Field::DIRECT  ),     4.629348, err);
  ASSERT_NEAR(interaction("2-body COM cor." , 0)->energy(NEUTRON, Field::PAIRING ),     0.201743, err);
  ASSERT_NEAR(interaction("coulomb (Slater)", 0)->energy(NEUTRON, Field::EXCHANGE),     0.000000, err);
  ASSERT_NEAR(interaction("central"         , 0)->energy(NEUTRON, Field::PAIRING ),     4.329008, err);
  ASSERT_NEAR(interaction("central"         , 1)->energy(NEUTRON, Field::PAIRING ),    -8.604896, err);
  ASSERT_NEAR(interaction("density FR"      , 0)->energy(NEUTRON, Field::DIRECT  ),  1165.672001, err);
  ASSERT_NEAR(interaction("density FR"      , 0)->energy(NEUTRON, Field::EXCHANGE),  -608.620103, err);
  ASSERT_NEAR(interaction("density FR"      , 0)->energy(NEUTRON, Field::PAIRING ),     0.000000, err);
  ASSERT_NEAR(interaction("rearrangement FR", 0)->energy(NEUTRON, Field::DIRECT  ),   194.278667, err);
  ASSERT_NEAR(interaction("rearrangement FR", 0)->energy(NEUTRON, Field::EXCHANGE),   -97.338078, err);
  ASSERT_NEAR(interaction("rearrangement FR", 0)->energy(NEUTRON, Field::PAIRING ),     0.000000, err);
  ASSERT_NEAR(interaction("kinetic"         , 0)->energy(PROTON , Field::DIRECT  ),   310.321578, err);
  ASSERT_NEAR(interaction("coulomb (Slater)", 0)->energy(PROTON , Field::DIRECT  ),    79.170771, err);
  ASSERT_NEAR(interaction("central"         , 0)->energy(PROTON , Field::DIRECT  ), -1014.983218, err);
  ASSERT_NEAR(interaction("central"         , 1)->energy(PROTON , Field::DIRECT  ),   -21.069868, err);
  ASSERT_NEAR(interaction("central"         , 0)->energy(PROTON , Field::EXCHANGE),   451.122647, err);
  ASSERT_NEAR(interaction("central"         , 1)->energy(PROTON , Field::EXCHANGE),  -516.202991, err);
  ASSERT_NEAR(interaction("spin-orbit"      , 0)->energy(PROTON , Field::DIRECT  ),    -1.685199, err);
  ASSERT_NEAR(interaction("2-body COM cor." , 0)->energy(PROTON , Field::DIRECT  ),     3.678276, err);
  ASSERT_NEAR(interaction("2-body COM cor." , 0)->energy(PROTON , Field::PAIRING ),     0.000000, err);
  ASSERT_NEAR(interaction("coulomb (Slater)", 0)->energy(PROTON , Field::EXCHANGE),    -7.449355, err);
  ASSERT_NEAR(interaction("central"         , 0)->energy(PROTON , Field::PAIRING ),     0.000000, err);
  ASSERT_NEAR(interaction("central"         , 1)->energy(PROTON , Field::PAIRING ),    -0.000000, err);
  ASSERT_NEAR(interaction("density FR"      , 0)->energy(PROTON , Field::DIRECT  ),  1067.084336, err);
  ASSERT_NEAR(interaction("density FR"      , 0)->energy(PROTON , Field::EXCHANGE),  -510.032438, err);
  ASSERT_NEAR(interaction("density FR"      , 0)->energy(PROTON , Field::PAIRING ),     0.000000, err);
  ASSERT_NEAR(interaction("rearrangement FR", 0)->energy(PROTON , Field::DIRECT  ),   177.847389, err);
  ASSERT_NEAR(interaction("rearrangement FR", 0)->energy(PROTON , Field::EXCHANGE),   -89.104013, err);
  ASSERT_NEAR(interaction("rearrangement FR", 0)->energy(PROTON , Field::PAIRING ),     0.000000, err);
}

//==============================================================================

TEST(Meta, energy_D2X_1ct)
{
  double err = 1e-6;

  DataTree dataTree("misc/data/rho_export_Fm_def.dat");
  State state(dataTree);

  Interaction interaction("D2X", &state);
  interaction.calcEnergies();

  // ALL_OUT;
  // INFO("%14.6f", interaction("kinetic"         , 0)->energy(NEUTRON, Field::DIRECT  ));
  // INFO("%14.6f", interaction("coulomb (exact)" , 0)->energy(NEUTRON, Field::DIRECT  ));
  // INFO("%14.6f", interaction("central"         , 0)->energy(NEUTRON, Field::DIRECT  ));
  // INFO("%14.6f", interaction("central"         , 1)->energy(NEUTRON, Field::DIRECT  ));
  // INFO("%14.6f", interaction("central"         , 0)->energy(NEUTRON, Field::EXCHANGE));
  // INFO("%14.6f", interaction("central"         , 1)->energy(NEUTRON, Field::EXCHANGE));
  // INFO("%14.6f", interaction("spin-orbit"      , 0)->energy(NEUTRON, Field::DIRECT  ));
  // INFO("%14.6f", interaction("2-body COM cor." , 0)->energy(NEUTRON, Field::DIRECT  ));
  // INFO("%14.6f", interaction("2-body COM cor." , 0)->energy(NEUTRON, Field::PAIRING ));
  // INFO("%14.6f", interaction("coulomb (exact)" , 0)->energy(NEUTRON, Field::EXCHANGE));
  // INFO("%14.6f", interaction("central"         , 0)->energy(NEUTRON, Field::PAIRING ));
  // INFO("%14.6f", interaction("central"         , 1)->energy(NEUTRON, Field::PAIRING ));
  // INFO("%14.6f", interaction("density FR"      , 0)->energy(NEUTRON, Field::DIRECT  ));
  // INFO("%14.6f", interaction("density FR"      , 0)->energy(NEUTRON, Field::EXCHANGE));
  // INFO("%14.6f", interaction("density FR"      , 0)->energy(NEUTRON, Field::PAIRING ));
  // INFO("%14.6f", interaction("rearrangement FR", 0)->energy(NEUTRON, Field::DIRECT  ));
  // INFO("%14.6f", interaction("rearrangement FR", 0)->energy(NEUTRON, Field::EXCHANGE));
  // INFO("%14.6f", interaction("rearrangement FR", 0)->energy(NEUTRON, Field::PAIRING ));
  // INFO("%14.6f", interaction("kinetic"         , 0)->energy(PROTON , Field::DIRECT  ));
  // INFO("%14.6f", interaction("coulomb (exact)" , 0)->energy(PROTON , Field::DIRECT  ));
  // INFO("%14.6f", interaction("central"         , 0)->energy(PROTON , Field::DIRECT  ));
  // INFO("%14.6f", interaction("central"         , 1)->energy(PROTON , Field::DIRECT  ));
  // INFO("%14.6f", interaction("central"         , 0)->energy(PROTON , Field::EXCHANGE));
  // INFO("%14.6f", interaction("central"         , 1)->energy(PROTON , Field::EXCHANGE));
  // INFO("%14.6f", interaction("spin-orbit"      , 0)->energy(PROTON , Field::DIRECT  ));
  // INFO("%14.6f", interaction("2-body COM cor." , 0)->energy(PROTON , Field::DIRECT  ));
  // INFO("%14.6f", interaction("2-body COM cor." , 0)->energy(PROTON , Field::PAIRING ));
  // INFO("%14.6f", interaction("coulomb (exact)" , 0)->energy(PROTON , Field::EXCHANGE));
  // INFO("%14.6f", interaction("central"         , 0)->energy(PROTON , Field::PAIRING ));
  // INFO("%14.6f", interaction("central"         , 1)->energy(PROTON , Field::PAIRING ));
  // INFO("%14.6f", interaction("density FR"      , 0)->energy(PROTON , Field::DIRECT  ));
  // INFO("%14.6f", interaction("density FR"      , 0)->energy(PROTON , Field::EXCHANGE));
  // INFO("%14.6f", interaction("density FR"      , 0)->energy(PROTON , Field::PAIRING ));
  // INFO("%14.6f", interaction("rearrangement FR", 0)->energy(PROTON , Field::DIRECT  ));
  // INFO("%14.6f", interaction("rearrangement FR", 0)->energy(PROTON , Field::EXCHANGE));
  // INFO("%14.6f", interaction("rearrangement FR", 0)->energy(PROTON , Field::PAIRING ));

  ASSERT_NEAR(interaction("kinetic"         , 0)->energy(NEUTRON, Field::DIRECT  ),   3274.078325, err);
  ASSERT_NEAR(interaction("coulomb (exact)" , 0)->energy(NEUTRON, Field::DIRECT  ),      0.000000, err);
  ASSERT_NEAR(interaction("central"         , 0)->energy(NEUTRON, Field::DIRECT  ), -14511.681319, err);
  ASSERT_NEAR(interaction("central"         , 1)->energy(NEUTRON, Field::DIRECT  ),    685.806531, err);
  ASSERT_NEAR(interaction("central"         , 0)->energy(NEUTRON, Field::EXCHANGE),  10437.859688, err);
  ASSERT_NEAR(interaction("central"         , 1)->energy(NEUTRON, Field::EXCHANGE),  -5364.201496, err);
  ASSERT_NEAR(interaction("spin-orbit"      , 0)->energy(NEUTRON, Field::DIRECT  ),    -66.138714, err);
  ASSERT_NEAR(interaction("2-body COM cor." , 0)->energy(NEUTRON, Field::DIRECT  ),      9.098913, err);
  ASSERT_NEAR(interaction("2-body COM cor." , 0)->energy(NEUTRON, Field::PAIRING ),      0.086315, err);
  ASSERT_NEAR(interaction("coulomb (exact)" , 0)->energy(NEUTRON, Field::EXCHANGE),      0.000000, err);
  ASSERT_NEAR(interaction("central"         , 0)->energy(NEUTRON, Field::PAIRING ),      3.034850, err);
  ASSERT_NEAR(interaction("central"         , 1)->energy(NEUTRON, Field::PAIRING ),    -10.467761, err);
  ASSERT_NEAR(interaction("density FR"      , 0)->energy(NEUTRON, Field::DIRECT  ),   9652.809572, err);
  ASSERT_NEAR(interaction("density FR"      , 0)->energy(NEUTRON, Field::EXCHANGE),  -5214.718808, err);
  ASSERT_NEAR(interaction("density FR"      , 0)->energy(NEUTRON, Field::PAIRING ),      1.824975, err);
  ASSERT_NEAR(interaction("rearrangement FR", 0)->energy(NEUTRON, Field::DIRECT  ),   1615.079466, err);
  ASSERT_NEAR(interaction("rearrangement FR", 0)->energy(NEUTRON, Field::EXCHANGE),   -796.159209, err);
  ASSERT_NEAR(interaction("rearrangement FR", 0)->energy(NEUTRON, Field::PAIRING ),      0.000000, err);
  ASSERT_NEAR(interaction("kinetic"         , 0)->energy(PROTON , Field::DIRECT  ),   1604.955641, err);
  ASSERT_NEAR(interaction("coulomb (exact)" , 0)->energy(PROTON , Field::DIRECT  ),   1131.795543, err);
  ASSERT_NEAR(interaction("central"         , 0)->energy(PROTON , Field::DIRECT  ), -10769.270880, err);
  ASSERT_NEAR(interaction("central"         , 1)->energy(PROTON , Field::DIRECT  ),    602.476271, err);
  ASSERT_NEAR(interaction("central"         , 0)->energy(PROTON , Field::EXCHANGE),   7583.481392, err);
  ASSERT_NEAR(interaction("central"         , 1)->energy(PROTON , Field::EXCHANGE),  -4560.317027, err);
  ASSERT_NEAR(interaction("spin-orbit"      , 0)->energy(PROTON , Field::DIRECT  ),    -35.968086, err);
  ASSERT_NEAR(interaction("2-body COM cor." , 0)->energy(PROTON , Field::DIRECT  ),      4.028044, err);
  ASSERT_NEAR(interaction("2-body COM cor." , 0)->energy(PROTON , Field::PAIRING ),      0.099620, err);
  ASSERT_NEAR(interaction("coulomb (exact)" , 0)->energy(PROTON , Field::EXCHANGE),    -38.339404, err);
  ASSERT_NEAR(interaction("central"         , 0)->energy(PROTON , Field::PAIRING ),      4.138216, err);
  ASSERT_NEAR(interaction("central"         , 1)->energy(PROTON , Field::PAIRING ),    -15.277117, err);
  ASSERT_NEAR(interaction("density FR"      , 0)->energy(PROTON , Field::DIRECT  ),   6372.610204, err);
  ASSERT_NEAR(interaction("density FR"      , 0)->energy(PROTON , Field::EXCHANGE),  -2685.505957, err);
  ASSERT_NEAR(interaction("density FR"      , 0)->energy(PROTON , Field::PAIRING ),      2.611214, err);
  ASSERT_NEAR(interaction("rearrangement FR", 0)->energy(PROTON , Field::DIRECT  ),   1055.823830, err);
  ASSERT_NEAR(interaction("rearrangement FR", 0)->energy(PROTON , Field::EXCHANGE),   -519.805885, err);
  ASSERT_NEAR(interaction("rearrangement FR", 0)->energy(PROTON , Field::PAIRING ),      0.000000, err);
}
