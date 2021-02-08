/** @file Apegrunt.h
	Top-level include file for Apegrunt.

	Copyright (c) 2016-2017 Santeri Puranen.

	This program is free software: you can redistribute it and/or modify
	it under the terms of the GNU Affero General Public License as
	published by the Free Software Foundation, either version 3 of the
	License, or (at your option) any later version.

	This program is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
	GNU Affero General Public License for more details.

	You should have received a copy of the GNU Affero General Public License
	along with this program. If not, see <http://www.gnu.org/licenses/>.

	@author Santeri Puranen
	$Id: $
*/

#ifndef APEGRUNT_H
#define APEGRUNT_H

#ifndef APEGRUNT_NO_TBB // Threading with Threading Building Blocks
#include "tbb/tbb_stddef.h"
#include "tbb/task_scheduler_init.h"
#endif // #ifndef APEGRUNT_NO_TBB

#include "Apegrunt_options.h"
#include "Apegrunt_utility.hpp"

#include "Alignment.h"
#include "Alignment_parsers.hpp"
#include "Alignment_generators.hpp"
#include "Alignment_impl_block_compressed_storage.hpp"
#include "StateVector_impl_block_compressed_alignment_storage.hpp"
#include "Alignment_utility.hpp"
#include "Apegrunt_IO_misc.hpp"
#include "ValueVector_parser.hpp"

#endif // APEGRUNT_H
