/** @file IntegerSequence_Hybrid_bitset_range_operations_forward.h

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

#ifndef APEGRUNT_INTEGERSEQUENCE_HYBRID_BITSET_RANGE_OPERATIONS_FORWARD_H
#define APEGRUNT_INTEGERSEQUENCE_HYBRID_BITSET_RANGE_OPERATIONS_FORWARD_H

#include <ostream>
#include <vector>

#include "IntegerSequence_forward.h"
#include "IntegerSequence_operations_forward.h"
#include "IntegerSequence_Hybrid_bitset_range_forward.h"

namespace apegrunt {

template< typename IndexT, bool Aligned=true >
IntegerSequence< Hybrid_bitset_range_element<IndexT,Aligned> > set_union( const IntegerSequence< Hybrid_bitset_range_element<IndexT,Aligned> >& a, const IntegerSequence< Hybrid_bitset_range_element<IndexT,Aligned> >& b );

template< typename IndexT >
IntegerSequence< Hybrid_bitset_range_element<IndexT,true> > set_union( const IntegerSequence< Hybrid_bitset_range_element<IndexT,true> >& a, const IntegerSequence< Hybrid_bitset_range_element<IndexT,true> >& b );

template< typename IndexT, bool Aligned=true >
IntegerSequence< Hybrid_bitset_range_element<IndexT,Aligned> > set_union( const std::vector< std::reference_wrapper< const IntegerSequence< Hybrid_bitset_range_element<IndexT,Aligned> > > >& sets );

template< typename IndexT >
IntegerSequence< Hybrid_bitset_range_element<IndexT,true> > set_union( const std::vector< std::reference_wrapper< const IntegerSequence< Hybrid_bitset_range_element<IndexT,true> > > >& sets );

template< typename IndexT, typename UnaryFunction, bool Aligned=true >
UnaryFunction inplace_set_differences_and_for_each_in_intersection( IntegerSequence< Hybrid_bitset_range_element<IndexT,Aligned> >& a, IntegerSequence< Hybrid_bitset_range_element<IndexT,Aligned> >& b, UnaryFunction f );

template< typename IndexT, typename UnaryFunction >
UnaryFunction inplace_set_differences_and_for_each_in_intersection( IntegerSequence< Hybrid_bitset_range_element<IndexT,true> >& a, IntegerSequence< Hybrid_bitset_range_element<IndexT,true> >& b, UnaryFunction f );

template< typename IndexT, bool Aligned=true >
struct hash< IntegerSequence< Hybrid_bitset_range_element<IndexT,Aligned> > >;

template< typename IndexT, bool Aligned=true >
struct equal_to< IntegerSequence< Hybrid_bitset_range_element<IndexT,Aligned> > >;

} // namespace apegrunt

#endif // APEGRUNT_INTEGERSEQUENCE_HYBRID_BITSET_RANGE_OPERATIONS_FORWARD_H
