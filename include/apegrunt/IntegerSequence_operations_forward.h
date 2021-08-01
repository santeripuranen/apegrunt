/** @file IntegerSequence_operations_forward.h

	Copyright (c) 2016-2021 Santeri Puranen.

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
template< typename ContainerT >
ContainerT set_union( const ContainerT& a, const ContainerT& b );

template< typename ContainerT >
ContainerT set_union( const std::vector< std::reference_wrapper<const ContainerT> > & sets );

template< typename ContainerT, typename UnaryFunction >
UnaryFunction inplace_set_differences_and_for_each_in_intersection( ContainerT& a, ContainerT& b, UnaryFunction f );

} // namespace apegrunt

#endif // APEGRUNT_INTEGERSEQUENCE_OPERATIONS_FORWARD_H
