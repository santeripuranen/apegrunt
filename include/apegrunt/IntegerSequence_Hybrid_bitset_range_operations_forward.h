/** @file IntegerSequence_Hybrid_bitset_range_operations_forward.h

	Copyright (c) 2018-2020 Santeri Puranen.

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

#ifndef APEGRUNT_INTEGERSEQUENCE_HYBRID_BITSET_RANGE_OPERATIONS_FORWARD_H
#define APEGRUNT_INTEGERSEQUENCE_HYBRID_BITSET_RANGE_OPERATIONS_FORWARD_H

#include <ostream>
#include <vector>
#include <functional> // for std::reference_wrapper

#include "IntegerSequence_forward.h"
#include "IntegerSequence_operations_forward.h"
#include "IntegerSequence_Hybrid_bitset_range_forward.h"

namespace apegrunt {

template< typename IndexT, bool Aligned >
std::ostream& operator<< ( std::ostream& os, const IntegerSequence< Hybrid_bitset_range_element<IndexT,Aligned> >& container );

template< typename IndexT >
std::ostream& operator<< ( std::ostream& os, const IntegerSequence< Hybrid_bitset_range_element<IndexT,true> >& container );

template< typename IndexT, bool Aligned=true >
Hybrid_bitset_sequence<IndexT,Aligned> set_intersection( const Hybrid_bitset_sequence<IndexT,Aligned>& a, const Hybrid_bitset_sequence<IndexT,Aligned>& b );

template< typename IndexT >
Hybrid_bitset_sequence<IndexT,true> set_intersection( const Hybrid_bitset_sequence<IndexT,true>& a, const Hybrid_bitset_sequence<IndexT,true>& b );

template< typename IndexT, typename RealT, bool Aligned=true >
RealT intersect_and_gather( const Hybrid_bitset_sequence<IndexT,Aligned>& a, const Hybrid_bitset_sequence<IndexT,Aligned>& b, const std::vector<RealT>& weights );

template< typename IndexT, typename RealT >
RealT intersect_and_gather( const Hybrid_bitset_sequence<IndexT,true>& a, const Hybrid_bitset_sequence<IndexT,true>& b, const std::vector<RealT>& weights );

template< typename IndexT, bool Aligned=true >
Hybrid_bitset_sequence<IndexT,Aligned> set_union( const Hybrid_bitset_sequence<IndexT,Aligned>& a, const Hybrid_bitset_sequence<IndexT,Aligned>& b );

template< typename IndexT >
Hybrid_bitset_sequence<IndexT,true> set_union( const Hybrid_bitset_sequence<IndexT,true>& a, const Hybrid_bitset_sequence<IndexT,true>& b );

template< typename IndexT, bool Aligned=true >
Hybrid_bitset_sequence<IndexT,Aligned> set_union( const std::vector< std::reference_wrapper< const Hybrid_bitset_sequence<IndexT,Aligned> > >& sets );

template< typename IndexT >
Hybrid_bitset_sequence<IndexT,true> set_union( const std::vector< std::reference_wrapper< const Hybrid_bitset_sequence<IndexT,true> > >& sets );

template< typename IndexT, bool Aligned >
struct hash< Hybrid_bitset_sequence<IndexT,Aligned> >;

template< typename IndexT, bool Aligned >
struct equal_to< Hybrid_bitset_sequence<IndexT,Aligned> >;

} // namespace apegrunt

#endif // APEGRUNT_INTEGERSEQUENCE_HYBRID_BITSET_RANGE_OPERATIONS_FORWARD_H
