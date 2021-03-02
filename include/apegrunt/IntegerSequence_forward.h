/** @file IntegerSequence_forward.h

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

#ifndef APEGRUNT_INTEGERSEQUENCE_FORWARD_H
#define APEGRUNT_INTEGERSEQUENCE_FORWARD_H

#include <memory> // for std::shared_ptr and std::make_shared

#include "Apegrunt_functional.h"

namespace apegrunt {

// forward declarations
template< typename SequenceT >
class IntegerSequence;

template< typename SequenceT >
using IntegerSequence_ptr = std::shared_ptr< IntegerSequence< SequenceT > >;

template< typename IntegerSequenceT, typename... Args >
IntegerSequence_ptr< typename IntegerSequenceT::value_type > make_IntegerSequence_ptr( Args&... args )
{
	return std::make_shared< IntegerSequenceT >( args... );
}

template< typename SequenceT > struct hash< IntegerSequence< SequenceT > >;

template< typename SequenceT > struct equal_to< IntegerSequence< SequenceT > >;

} // namespace apegrunt

#endif // APEGRUNT_INTEGERSEQUENCE_FORWARD_H
