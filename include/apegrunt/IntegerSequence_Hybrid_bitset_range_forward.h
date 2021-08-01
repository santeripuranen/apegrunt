/** @file IntegerSequence_Hybrid_bitset_range_forward.h

	Copyright (c) 2018-2021 Santeri Puranen.

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

#ifndef APEGRUNT_INTEGERSEQUENCE_HYBRID_BITSET_RANGE_FORWARD_H
#define APEGRUNT_INTEGERSEQUENCE_HYBRID_BITSET_RANGE_FORWARD_H

#include <ostream>
#include <vector>

#include "IntegerSequence_forward.h"

namespace apegrunt {

template< typename IndexT, bool Aligned=true >
struct Hybrid_bitset_range_element;

template< typename IndexT, bool Aligned=true >
using Apegrunt_bitset = IntegerSequence< Hybrid_bitset_range_element<IndexT,Aligned> >;

template< typename IndexT, bool Aligned > struct hash< Hybrid_bitset_range_element<IndexT,Aligned> >;
template< typename IndexT, bool Aligned > struct equal_to< Hybrid_bitset_range_element<IndexT,Aligned> >;

} // namespace apegrunt

#endif // APEGRUNT_INTEGERSEQUENCE_HYBRID_BITSET_RANGE_FORWARD_H
