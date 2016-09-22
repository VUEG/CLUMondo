/*****************************************************************
 * This file is part of CLUMondo.
 *
 * CLUMondo is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * any later version.
 *
 * CLUMondois distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with CLUMondo.  If not, see <http://www.gnu.org/licenses/>.
 *
 * CLUMondo utilities
 * CLUMondo S.1.1
 *
 * Development:
 *	Peter Verburg
 *	31.3.2015
 *
 * All rights: Peter Verburg <peter.verburg@vu.nl>
 *
 ********************************************************************/
#include "../include/util.h"

#include <iostream>

void show_error()
{
	std::cout << "Error - see LOG.FIL for details" << std::endl;
}

void show_finished()
{
	std::cout << "Finished" << std::endl;
}

void show_usage()
{
	std::cout << "Usage:" << std::endl;
	std::cout << "CLUMondo.exe demand_file region_file [model_mode]" << std::endl;
	std::cout << "  Valid values for the optional parameter [model mode] are:" << std::endl;
	std::cout << "    - 1: use meta model mode" << std::endl;
	std::cout << "    - 2: create all probability maps in advance" << std::endl;
	std::cout << "    - 3: write yearly probability map" << std::endl;
	std::cout << "    - 4: only write initial probability map" << std::endl;
}