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

#include "system.h"

std::list<KeyStruct > System::validKeys =
  {
    { "system/nNeut", "Neutron number", "94"   , "I"  },
    { "system/nProt", "Proton number" , "146"  , "I"  }
  };

//==============================================================================
//==============================================================================
//==============================================================================

/** The constructor from a DataTree instance.
 */

System::System(const DataTree &dataTree)
{
  dataTree.get(nProt, "system/nProt", true);
  dataTree.get(nNeut, "system/nNeut", true);
}

//==============================================================================
//==============================================================================
//==============================================================================

/** The constructor from a number of protons and a number of neutrons.
 */

System::System(INT _nProt, INT _nNeut) : nProt(_nProt), nNeut(_nNeut)
{
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Return a DataTree instance.
 */

DataTree System::getDataTree(void) const
{
  DBG_ENTER;

  DataTree result;

  result.set("system/nProt", nProt);
  result.set("system/nNeut", nNeut);

  DBG_RETURN(result);
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Return a nice info string.
 */

const std::string System::info(bool isShort) const
{
  DBG_ENTER;

  std::string result = "";
  result += PF("%d", nProt + nNeut) + getElementSymbol(nProt);

  if (!isShort)
  {
    result += " (" + PF_YELLOW("Z") + ":" + PF_BLUE("%d", nProt) + ",";
    result +=        PF_YELLOW("N") + ":" + PF_BLUE("%d", nNeut) + ")";
  }

  DBG_RETURN(result);
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Return the symbol of an element.
 */

const std::string System::getElementSymbol(INT Z)
{
  DBG_ENTER;

  std::string elemName[] = {"H",                                                                                                 "He",
                            "Li", "Be",                                                             "B", "C", "N", "O", "F", "Ne",
                            "Na", "Mg",                                                             "Al", "Si", "P", "S", "Cl", "Ar",
                            "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr",
                            "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe",
                            "Cs", "Ba", "La",
                            "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu",
                            "Hf", "Ta", "W", "Re", "Os", "Il", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn",
                            "Fr", "Ra", "Ac",
                            "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr",
                            "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og"
                           };
  std::string name = "X";

  if ((Z >= 1) && (Z <= 118)) name = elemName[Z - 1];

  DBG_RETURN(name);
}
