/** @file Value_vector_forward.h

	Copyright (c) 2018 Santeri Puranen.

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

#ifndef APEGRUNT_VALUE_VECTOR_FORWARD_H
#define APEGRUNT_VALUE_VECTOR_FORWARD_H

#include <memory> // std::shared_ptr, std::make_shared

namespace apegrunt {

// forward declarations
template< typename ValueT >
using ValueVector = std::vector<ValueT>;

template< typename ValueT >
using ValueVector_ptr = std::shared_ptr<ValueVector<ValueT> >;

template< typename ValueT, typename... Args >
ValueVector_ptr<ValueT> make_ValueVector_ptr( Args&&... args )
{
	return std::make_shared< ValueVector<ValueT> >( std::forward<Args>(args)... );
}

template< typename ValueT >
ValueVector_ptr<ValueT> make_ValueVector_ptr( const ValueVector<ValueT>& value_vector )
{
	return std::make_shared< ValueVector<ValueT> >( value_vector );
}

template< typename ValueT >
typename apegrunt::ValueVector<ValueT>::iterator begin( apegrunt::ValueVector_ptr<ValueT>& ptr ) { using std::begin; return begin( *ptr ); }

template< typename ValueT >
typename apegrunt::ValueVector<ValueT>::iterator end( apegrunt::ValueVector_ptr<ValueT>& ptr ) { using std::end; return end( *ptr ); }

template< typename ValueT >
typename apegrunt::ValueVector<ValueT>::const_iterator cbegin( const apegrunt::ValueVector_ptr<ValueT>& ptr ) { using std::cbegin; return cbegin( *ptr ); }

template< typename ValueT >
typename apegrunt::ValueVector<ValueT>::const_iterator cend( const apegrunt::ValueVector_ptr<ValueT>& ptr ) { using std::cend; return cend( *ptr ); }

} // namespace apegrunt

#endif // APEGRUNT_VALUE_VECTOR_FORWARD_H

