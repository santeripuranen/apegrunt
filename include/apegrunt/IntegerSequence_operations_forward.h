/** @file IntegerSequence_operations_forward.h

	Copyright (c) 2016-2020 Santeri Puranen.

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

#ifndef APEGRUNT_INTEGERSEQUENCE_OPERATIONS_FORWARD_H
#define APEGRUNT_INTEGERSEQUENCE_OPERATIONS_FORWARD_H

#include <vector>
#include <ostream>
#include <functional> // for std::reference_wrapper

#include "IntegerSequence_forward.h"

namespace apegrunt {

// forward declarations

template< typename SequenceT >
std::ostream& operator<< ( std::ostream& os, const IntegerSequence< SequenceT >& container );

template< typename ContainerT >
ContainerT set_intersection( const ContainerT& a, const ContainerT& b );

// sugar: shorthand for set_intersection; however, we shouldn't declare anything with such a generic signature here, or else there will be trouble elsewhere
//template< typename ContainerT >
//ContainerT operator&( const ContainerT& a, const ContainerT& b );

template< typename SequenceT, typename RealT=double >
RealT intersect_and_gather( const IntegerSequence< SequenceT >& a, const IntegerSequence< SequenceT >& b, const std::vector<RealT>& weights );

template< typename ContainerT >
ContainerT set_union( const ContainerT& a, const ContainerT& b );

// sugar: shorthand for set_union; however, we shouldn't declare anything with such a generic signature here, or else there will be trouble elsewhere
//template< typename ContainerT >
//ContainerT operator|( const ContainerT& a, const ContainerT& b );

template< typename ContainerT >
ContainerT set_union( const std::vector< std::reference_wrapper<const ContainerT> > & sets );

} // namespace apegrunt

#endif // APEGRUNT_INTEGERSEQUENCE_OPERATIONS_FORWARD_H
