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
 *  \brief Test suite for the Base64 namespace.
 */

#include "gtest/gtest.h"
#include "base64.h"

//==============================================================================

TEST(Base64, encode)
{
  std::string input = "long live HFB3 !";
  std::string result = Base64::encode(input.c_str(), input.size());
  std::string expected = "bG9uZyBsaXZlIEhGQjMgIQ==";
  ASSERT_EQ(result, expected);
}

//==============================================================================

TEST(Base64, encode_)
{
  std::string input = "long live HFB3 !";
  std::string result = Base64::encode(input);
  std::string expected = "bG9uZyBsaXZlIEhGQjMgIQ==";
  ASSERT_EQ(result, expected);
}

//==============================================================================

TEST(Base64, decode)
{
  std::string input = "bG9uZyBsaXZlIEhGQjMgIQ==";
  std::string result = Base64::decode(input);
  std::string expected = "long live HFB3 !";
  ASSERT_EQ(result, expected);
}

