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

#include "logger.h"

/** \file
 *  \brief Methods of the Logger class.
 */

//==============================================================================
//==============================================================================
//==============================================================================

/** Dummy constructor.
 */

Logger::Logger(void) : _callback(new Callback)
{
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Print a log message.
 */

void Logger::log(const std::string &mesg) const
{
  if (_callback) _callback->run(mesg);
  else std::cout << "no callback function defined in Logger::log()" << std::endl;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Set a CallBack function.
 */

void Logger::setCallback(Callback *cb)
{
  delCallback();
  _callback = cb;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Set a CallBack function.
 */

void Logger::delCallback(void)
{
  delete _callback;
  _callback = 0;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Destructor.
 */

Logger::~Logger()
{
  // FIXME: commenting the following line causes a (small) memory leak,
  // but avoids a Python bug. This should be fixed...

  // delCallback();
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Callback destructor.
 */

Callback::~Callback()
{
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Default Callback action.
 */

void Callback::run(const std::string &mesg)
{
  std::cout << mesg << std::endl;
}

//==============================================================================
//==============================================================================
//==============================================================================


