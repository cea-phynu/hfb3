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

#include "qnumbers.h"
#include "tools.h"

/** \file
 *  \brief Methods of the Qnumbers class.
 */

//==============================================================================
//==============================================================================
//==============================================================================

/** The constructor.
 *
 * \param _n Number of quantum numbers.
 */

Qnumbers::Qnumbers(UINT _n)
{
  nbq = _n > 0 ? _n : 0;
  nb = 0;
  name = std::vector<std::string>(nbq);
}

//==============================================================================
//==============================================================================
//==============================================================================

/** The constructor for a single linear index.
 *
 * \param _name Name of the single quantum number.
 * \param _nb Number of values of the single quantum number.
 */

Qnumbers::Qnumbers(std::string _name, UINT _nb)
{
  nbq = 0;
  nb = _nb;
  mat = IMAT(1, _nb, arma::fill::zeros);
  appendName(_name);

  for (INT i = 0; i < _nb; i++)
  {
    mat(0, i) = i;
  }
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Set the names of the quantum numbers.
 *
 *  \param l List of names.
 */

void Qnumbers::setNames(const std::vector<std::string> &l)
{
  if (nbq == 0)
  {
    nbq = l.size();
  }

  ASSERT(l.size() == nbq, "bad dimension");

  name = l;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Find a state by its quantum numbers.
 *
 * \param i First quantum number value.
 * \param j Second quantum number value.
 * \param k Third quantum number value.
 * \param l Fourth quantum number value.
 * \param m Fifth quantum number value.
 */

INT Qnumbers::find(const IVEC &m) const
{
  try
  {
    return INT(index.at(m));
  }
  catch (const std::out_of_range &)
  {
    return -1;
  }
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Append a set of quantum numbers.
 *
 * \param i First quantum number value.
 * \param j Second quantum number value.
 * \param k Third quantum number value.
 * \param l Fourth quantum number value.
 * \param m Fifth quantum number value.
 */

UINT Qnumbers::append(const IVEC &v)
{
  mat.resize(nbq, nb + 1);
  mat.col(nb) = v.head(nbq);
  index[v] = nb;
  nb++;
  return nb;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Reorder the vectors.
 *
 *  \param r A transformation matrix.
 */

void Qnumbers::reorder(IMAT r)
{
  IMAT result = mat;

  for (UINT i = 0; i < nbq; i++)
  {
    mat.row(i) = r * result.row(i);
  }

  index.clear();

  for (UINT i2 = 0; i2 < nb; i2++)
  {
    index[mat.col(i2)] = i2;
  }
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Quantum number value access.
 *
 * \param i Quantum number id.
 * \param j State id;
 */

INT Qnumbers::operator()(UINT i, UINT j) const
{
  return INT(mat(i, j));
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Quantum number value access.
 *
 * \param i State id;
 */

IVEC Qnumbers::operator()(UINT i) const
{
  return mat.col(i);
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Append a quantum number name.
 *
 * \param s Name;
 */

void Qnumbers::appendName(const std::string &s)
{
  name.push_back(s);
  nbq++;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Return a subset of quantum numbers.
 *
 * \param ids State ids;
 */

const Qnumbers Qnumbers::sub(UVEC ids) const
{
  Qnumbers toto(nbq);
  toto.name = name;

  for (UINT j = 0; j < ids.n_rows; j++)
  {
    long id = ids(j);
    toto.append(mat.col(id));
  }

  return toto;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Return a single-valued Qnumbers from a multi-valued Qnumbers object.
 *
 * \param i State id;
 */

const Qnumbers Qnumbers::extract(UINT i) const
{
  Qnumbers toto(nbq);
  toto.name = name;
  toto.nb = 1;
  toto.mat = mat.col(i);
  return toto;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Split quantum numbers in blocks.
 *
 * \param rows Rows to be grouped.
 */

UINT Qnumbers::calcBlocks(const UVEC &rows)
{
  DBG_ENTER;

  blocks.clear();
  IMAT mat2 = mat.rows(rows);
  INT nbBlocks = 0;
  IVEC blockId = arma::zeros<IVEC >(nb);
  IVEC blockIndex = arma::zeros<IVEC >(nb);
  IVEC blockSize = arma::zeros<IVEC >(nb);

  for (INT i = 0; i < (INT)nb; i++)
  {
    bool first = true;

    for (INT j = 0; j < i; j++)
    {
      if (arma::all(mat2.col(i) == mat2.col(j)))
      {
        blockId(i) = blockId(j);
        blockIndex(i) = blockSize(blockId(j));
        blockSize(blockId(j))++;
        first = false;
        break;
      }
    }

    if (first)
    {
      blockId(i) = nbBlocks;
      blockIndex(i) = 0;
      blockSize(nbBlocks) = 1;
      nbBlocks++;
    }
  }

  arma::field<UVEC> filters(nbBlocks);

  for (INT i = 0; i < nbBlocks; i++)
  {
    filters(i) = UVEC(blockSize(i));
  }

  for (INT i = 0; i < (INT)nb; i++)
  {
    filters(blockId(i))(blockIndex(i)) = i;
  }

  for (UINT i = 0; i < filters.n_rows; i++)
  {
    Qnumbers newQN;
    newQN.name = name;
    newQN.mat = mat.cols(filters(i));
    newQN.nb = newQN.mat.n_cols;
    newQN.nbq = nbq;
    newQN.filter = filters(i);

    for (UINT j = 0; j < newQN.nb; j++)
      newQN.index[newQN(j)] = j;

    blocks.push_back(newQN);
  }

  DBG_RETURN(blocks.size());
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Check if a matrix preserves Omega-symmetry.
 */

void Qnumbers::checkOmegaSym(const arma::mat &mat, const std::string &name)
{
  DBG_ENTER;

  UINT nbBlocks = calcBlocks({0, 4}); // get blocks of same (m, s) QN

  for (UINT msa = 0; msa < nbBlocks; msa++)
  {
    Qnumbers &qna = blocks[msa];
    INT ma = qna(0, 0);
    INT sa = qna(4, 0);
    INT oa = ma - sa;

    for (UINT msb = 0; msb < nbBlocks; msb++)
    {
      Qnumbers &qnb = blocks[msb];
      INT mb = qnb(0, 0);
      INT sb = qnb(4, 0);
      INT ob = mb - sb;

      if (oa != ob)
      {
        if (arma::norm(mat.submat(qna.filter, qnb.filter), "inf") > 1e-12) ERROR("broken Omega-symmetry in " + name + ".");
      }
    }
  }

  Tools::debug("Omega-symmetry in matrix " + name + ": OK");

  DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Return a nice std::string representation of a state.
 */

const std::string Qnumbers::niceStr(UINT id) const
{
  std::string result = "{";

  for (UINT i = 0; i < nbq; i++)
  {
    result += PF_BLUE(name[i]) + PF(": ") + PF("%d", mat(i, id));

    if (i != nbq - 1)
    {
      result += ", ";
    }
  }

  result += "}";

  return result;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Clear object.
 */

void Qnumbers::clear(void)
{
  nbq = 0;
  nb = 0;
  name.clear();
  blocks.clear();
  filter.clear();
  index.clear();
  mat.clear();
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Return some info.
 */

const std::string Qnumbers::info(bool isShort) const
{
  DBG_ENTER;

  std::string result = "";

  if (isShort)
  {
    result += Tools::treeStr(
    {
      {"nbq", Tools::infoStr(nbq)},
      {"nb", Tools::infoStr(nb)},
      {"filter", filter.is_empty() ? "no" : "yes"},
      {"blocks", Tools::infoStr((INT)blocks.size())},
    }, true);
  }
  else
  {
    std::vector<std::pair<std::string, std::string> > values;
    values.push_back({"Qnumbers", ""});
    values.push_back({"#qnumb", Tools::infoStr(nbq)});
    values.push_back({"qnames", Tools::infoStr(name)});
    values.push_back({"#state", Tools::infoStr(nb)});
    values.push_back({"filter", Tools::infoStr(filter.is_empty())});
    values.push_back({"blocks", Tools::infoStr((INT)blocks.size())});

    for (UINT i = 0; i < nb; i++)
    {
      values.push_back({PF("state %d", i), niceStr(i)});
    }

    result += Tools::treeStr(values, false);
  }

  DBG_RETURN(result);
}
