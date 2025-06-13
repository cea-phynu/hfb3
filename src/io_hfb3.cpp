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

#include "io_hfb3.h"
#include "tools.h"
#include <iomanip>

/** \file
 *  \brief Methods of the IOhfb3 class.
 */

//==============================================================================
//==============================================================================
//==============================================================================

/** The constructor with a file to load.
 *
 */

IOhfb3::IOhfb3(void)
{
  DBG_ENTER;
  DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Return a DataTree instance.
 */

DataTree IOhfb3::fromContent(const std::string &content) const
{
  DBG_ENTER;
  DataTree result;
  std::stringstream stream;
  stream << content;
  std::string section = "";
  // magic number
  std::string inputStr;
  std::getline(stream, inputStr);

  if (inputStr != "#!HFB3!#")
  {
    std::stringstream ss;
    ss << "not an HFB3 format in IOhfb3::fromContent()";
    Tools::debug(ss.str());
    DBG_RETURN(DataTree());
  }

  while (!stream.eof())
  {
    std::getline(stream, inputStr);

    if (inputStr.length() == 0) continue;

    if (inputStr[0] == '#') continue;

    if (inputStr[0] == '[')
    {
      UINT vlast = inputStr.find_first_of("]");

      if (vlast > 1)
      {
        section = inputStr.substr(1, vlast - 1) + "/";
      }

      continue;
    }

    INT id = 0;
    UINT vfirst = 0;

    for (UINT i = 0; i < inputStr.size(); i++)
    {
      if ((id == 0) && ((inputStr[i] == ' ') || (inputStr[i] == '=') || (inputStr[i] == ':'))) id = 1;

      if ((id == 1) && ((inputStr[i] != ' ') && (inputStr[i] != '=') && (inputStr[i] != ':')))
      {
        vfirst = i;
        break;
      }
    }

    UINT vlast = inputStr.find_first_of(" =:");
    std::string key = "";
    std::string val = "";
    std::string type = "";

    if ((vfirst == 0) || (vlast == -1))
    {
      key = section + Tools::trim(inputStr.substr(0, vlast));
      val = "";
      type = "empty";
    }
    else
    {
      key = section + Tools::trim(inputStr.substr(0, vlast));
      val = Tools::trim(inputStr.substr(vfirst));
      type = "";

      if (val[0] == '[')
      {
        UINT tlast = val.find_first_of("]");
        type = Tools::trim(val.substr(1, tlast - 1));
        val = Tools::trim(val.substr(tlast + 1));
      }
      if (val[0] == '(')
      {
        UINT tlast = val.find_first_of(")");
        type = "VD";
        val = Tools::trim(val.substr(1, tlast - 1));
      }
    }

    updateDataTree(result, key, val, type);
  }

  DBG_RETURN(result);
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Create or update a (key, value) couple in a DataTree from a std::string representation.
 *
 *  \param key The key.
 *  \param val The value.
 *  \param type The type of the value (can be deduced from the value).
 */

void IOhfb3::updateDataTree(DataTree &dataTree, const std::string key, const std::string val, const std::string type)
{
  DBG_ENTER;

  if (type == "")
  {
    INT ival;
    double dval;
    std::stringstream st(val);
    st >> std::ws;
    st >> ival;

    if (!st.fail() && st.eof())
    {
      dataTree.set(key, ival);
      goto end_label;
    }
    else
    {
      std::stringstream st2(val);
      st2 >> std::ws;
      st2 >> std::setprecision(16) >> std::scientific >> dval;

      if (!st2.fail() && st2.eof())
      {
        dataTree.set(key, dval);
        goto end_label;
      }
      else
      {
        dataTree.set(key, val);
        goto end_label;
      }
    }
  }

  if (type == "empty")
  {
    dataTree.set(key);
  }
  else if (type == "I")
  {
    std::stringstream st(val);
    INT ival;
    st >> ival;
    st >> std::ws;

    if (!st.fail() && st.eof())
    {
      dataTree.set(key, ival);
    }
    else
    {
      std::stringstream ss;
      ss << "wrong INT value in '" << val << "'";
      Tools::warning(ss.str());
    }
  }
  else if (type == "D")
  {
    std::stringstream st(val);
    double dval;
    st >> std::setprecision(16) >> dval;
    st >> std::ws;

    if (!st.fail() && st.eof())
    {
      dataTree.set(key, dval);
    }
    else
    {
      std::stringstream ss;
      ss << "wrong double value in '" << val << "'";
      Tools::warning(ss.str());
    }
  }
  else if (type == "S")
  {
    dataTree.set(key, val);
  }
  else if (type == "VI")
  {
    // count commas
    INT nbCommas = 0;
    for (auto &c: val) if (c == ',') nbCommas++;

    INT nbValues = nbCommas + 1;
    IVEC tempIVec(nbValues);

    if (nbValues != 0)
    {
      INT lastPos = 0;
      INT idVal = 0;

      std::string val2 = val + ",";

      for (INT ic = 0; ic < val2.size(); ic++)
      {
        if (val2[ic] == ',')
        {
          std::string subVal = val2.substr(lastPos, ic - lastPos);
          std::stringstream st(subVal);
          INT ival;
          st >> ival;

          if (st.fail()) Tools::warning("wrong int value in '" + val2 + "' (" + subVal + ")");

          tempIVec[idVal] = ival;

          lastPos = ic + 1;
          idVal++;
        }
      }

      dataTree.set(key, tempIVec);
    }
  }
  else if (type == "VD")
  {
    // count commas
    INT nbCommas = 0;
    for (auto &c: val) if (c == ',') nbCommas++;

    INT nbValues = nbCommas + 1;
    arma::vec tempVec(nbValues);

    if (nbValues != 0)
    {
      INT lastPos = 0;
      INT idVal = 0;

      std::string val2 = val + ",";

      for (INT ic = 0; ic < val2.size(); ic++)
      {
        if (val2[ic] == ',')
        {
          std::string subVal = val2.substr(lastPos, ic - lastPos);
          std::stringstream st(subVal);
          double dval;
          st >> dval;

          if (st.fail()) Tools::warning("wrong double value in '" + val2 + "' (" + subVal + ")");

          tempVec[idVal] = dval;

          lastPos = ic + 1;
          idVal++;
        }
      }

      dataTree.set(key, tempVec);
    }
  }
  // else if (type == "VI")
  // {
  //   std::string decoded = Base64::decode(val);
  //   UINT nb = decoded.size() / sizeof(INT);
  //   IVEC tempIVec(nb);
  //
  //   for (INT i = 0; i < nb; i++)
  //     tempIVec(i) = ((INT *)(decoded.data()) )[i];
  //
  //   dataTree.set(key, tempIVec);
  // }
  // else if (type == "VD")
  // {
  //   std::string decoded = Base64::decode(val);
  //   UINT nb = decoded.size() / sizeof(double);
  //   arma::vec tempVec(nb);
  //
  //   for (INT i = 0; i < nb; i++)
  //     tempVec(i) = ((double *)(decoded.data()) )[i];
  //
  //   dataTree.set(key, tempVec);
  // }
  else if (type == "V")
  {
    std::stringstream s;
    s << Base64::decode(val);
    arma::vec temp;
    temp.load(s);
    dataTree.set(key, temp);
  }
  else if ((type == "M") || (type == "AM"))
  {
    std::stringstream s;
    s << Base64::decode(val);
    arma::mat temp;
    temp.load(s);
    dataTree.set(key, temp);
  }
  else if (type == "C")
  {
    std::stringstream s;
    s << Base64::decode(val);
    arma::cube temp;
    temp.load(s);
    dataTree.set(key, temp);
  }
  else if (type == "IV")
  {
    std::stringstream s;
    s << Base64::decode(val);
    IVEC temp;
    temp.load(s);
    dataTree.set(key, temp);
  }
  else if (type == "IM")
  {
    std::stringstream s;
    s << Base64::decode(val);
    IMAT temp;
    temp.load(s);
    dataTree.set(key, temp);
  }
  else if (type == "IC")
  {
    std::stringstream s;
    s << Base64::decode(val);
    ICUBE temp;
    temp.load(s);
    dataTree.set(key, temp);
  }

end_label:
  DBG_LEAVE;
}
