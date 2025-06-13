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

#ifndef LOGGER_H
#define LOGGER_H

/** \file
 *  \brief Headers for the Logger class.
 */

#include <iostream>
#include <string>

/** \brief A pure virtual class used in a callback technique.
 *
 * This class is used to receive a callback function from Python.
 */

class Callback
{
public:
  virtual ~Callback();
  virtual void run(const std::string &mesg);
};

//==============================================================================
//==============================================================================
//==============================================================================

/** \brief A logging manager.
 *
 * This class implements a logging manager that can be overloaded by Python.
 */

class Logger
{
public :
  Logger(void);
  ~Logger();

  void log(const std::string &mesg = "") const;
  void setCallback(Callback *cb);
  void delCallback();

private:
  Callback *_callback;
};

#endif // LOGGER_H
