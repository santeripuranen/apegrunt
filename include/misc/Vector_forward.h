/** @file Vector_forward.h

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

#ifndef APEGRUNT_VECTOR_FORWARD_H
#define APEGRUNT_VECTOR_FORWARD_H

namespace apegrunt {

// forward declarations
template< typename RealT, std::size_t N, bool View >
struct Vector;

template< typename RealT, std::size_t Capacity>
inline constexpr Vector<RealT,Capacity,false> make_Vector( RealT* const p, std::size_t stride=1 )
{
	return Vector<RealT,Capacity,false>(p,stride);
}

template< typename RealT, std::size_t Capacity>
inline constexpr Vector<RealT,Capacity,false> make_Vector( const RealT& e )
{
	return Vector<RealT,Capacity,false>(e);
}

template< typename RealT, std::size_t Capacity>
inline Vector<RealT,Capacity,true> make_Vector_view( RealT* const p, std::size_t stride=1 )
{
	return Vector<RealT,Capacity,true>(p,stride);
}

} // namespace apegrunt

#endif // APEGRUNT_VECTOR_FORWARD_H

