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

#include <iomanip>
#include "io_json.h"
#include "tools.h"
#include "datatree.h"
#include "base64.h"

/** \file
 *  \brief Methods of the IOjson class.
 */

//==============================================================================
//==============================================================================
//==============================================================================

/** The empty constructor.
 *
 */

IOjson::IOjson(void)
{
  DBG_ENTER;
  DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Return a DataTree instance.
 */

DataTree IOjson::fromContent(const std::string &_content) const
{
  DBG_ENTER;

  DataTree result;

  if ((_content[0] != '{') or (_content[_content.size() - 1] != '}'))
  {
    Tools::debug("invalid JSON format: missing starting '{' or ending '}' in IOjson::loadDataTree()");
    DBG_RETURN(DataTree());
  }

  std::string content = _content.substr(1, _content.size() - 2);
  std::vector<std::string> elems = Tools::stringSplit(content, ',');

  for (auto &e : elems)
  {
    std::string e2 = Tools::trim(e);
    std::vector<std::string> words = Tools::stringSplit(e2, ':');

    if (words.size() == 1)
    {
      std::string key = Tools::trim(words[0]);
      updateDataTree(result, key, "", "empty");
    }
    else if (words.size() == 2)
    {
      std::string key = Tools::trim(words[0]);
      std::string val = Tools::trim(words[1]);
      updateDataTree(result, key, val, "");
    }
    else
    {
      std::stringstream ss;
      ss << "invalid JSON format: missing or multiple ':' character(s) in '" << e2 << "' in IOjson::loadDataTree()";
      Tools::warning(ss.str());
      DBG_RETURN(DataTree());
    }
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

void IOjson::updateDataTree(DataTree &dataTree, const std::string key, const std::string val, const std::string type) const
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
    PF("set empty key '%s'", key.c_str());
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
    std::string decoded = Base64::decode(val);
    UINT nb = decoded.size() / sizeof(INT);
    //    std::vector<INT> newTab;
    //    for (INT i = 0; i < nb; i++) newTab.push_back( ((INT *)(decoded.data()) )[i]);
    //    vimap[key] = newTab;
    IVEC tempIVec(nb);

    for (INT i = 0; i < nb; i++)
      tempIVec(i) = ((INT *)(decoded.data()))[i];

    dataTree.set(key, tempIVec);
  }
  else if (type == "VD")
  {
    std::string decoded = Base64::decode(val);
    UINT nb = decoded.size() / sizeof(double);
    //    std::vector<double> newTab;
    //    for (INT i = 0; i < nb; i++) newTab.push_back(((double *)decoded.data())[i]);
    //    vdmap[key] = newTab;
    arma::vec tempVec(nb);

    for (INT i = 0; i < nb; i++)
      tempVec(i) = ((double *)(decoded.data()) )[i];

    dataTree.set(key, tempVec);
  }
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

