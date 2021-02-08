/** @file IntegerSequence_interface.hpp
 
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

#ifndef APEGRUNT_INTEGERSEQUENCE_OPERATIONS_HPP
#define APEGRUNT_INTEGERSEQUENCE_OPERATIONS_HPP

#include <vector>
#include <iterator> // for std::back_inserter
#include <ostream>
#include <limits>
#include <algorithm> // for std::set_intersection

#include "IntegerSequence_forward.h"
#include "IntegerSequence_interface.hpp"

namespace apegrunt {

template< typename RangeT >
typename IntegerSequence<RangeT>::const_iterator cbegin( const IntegerSequence<RangeT>& container ) { return container.cbegin(); }

template< typename RangeT >
typename IntegerSequence<RangeT>::const_iterator cend( const IntegerSequence<RangeT>& container ) { return container.cend(); }

template< typename RangeT >
typename IntegerSequence<RangeT>::iterator begin( const IntegerSequence<RangeT>& container ) { return container.begin(); }

template< typename RangeT >
typename IntegerSequence<RangeT>::iterator end( const IntegerSequence<RangeT>& container ) { return container.end(); }

template< typename T >
std::size_t bytesize( const IntegerSequence<T>& v ) { return v.bytesize(); }

template< typename T >
inline IntegerSequence< T >& operator<< ( IntegerSequence< T >& container, typename T::value_type value )
{
	container.push_back(value); return container;
}

template< typename T, typename ImplicitT >
inline std::vector< T >& operator<< ( std::vector< T >& container, ImplicitT value )
{
	container.push_back( T(value) ); return container;
}

template< typename T >
const std::vector<T>& get_dense( const std::vector<T>& container ) { return container; }

template< typename T >
std::vector<typename T::value_type> get_dense( const IntegerSequence< T >& container )
{
	std::vector<typename IntegerSequence<T>::value_type> dense;
	dense.reserve( container.size() );
	for( auto i: container ) { dense.push_back(i); }
	return dense;
}

template< typename IterableT >
std::ostream& encode_range_unaligned(std::ostream& os, const IterableT& container )
{
	using std::cbegin; using std::cend;
	using apegrunt::cbegin; using apegrunt::cend;

	auto currpos = cbegin( container );
	auto begpos = currpos;
	auto oldpos = currpos;
	const auto endpos = cend( container );

	os << std::hex << *currpos;
	while( currpos != endpos )
	{
		++currpos;
		if( *currpos > (*oldpos)+1 )
		{
			if( oldpos != begpos )
			{
				os << "-" << *oldpos;
			}
			os << "," << *currpos;
			begpos = currpos;
		}
		oldpos = currpos;
	}
	return os;
}

template< typename IndexT >
std::ostream& operator<< ( std::ostream& os, const std::vector< IndexT >& container )
{
	return encode_range_unaligned( os, container );
}

template< typename RangeT >
std::ostream& operator<< ( std::ostream& os, const IntegerSequence< RangeT >& container )
{
	return encode_range_unaligned( os, container );
}

template< typename RangeT >
struct unaligned_encoding_wrap
{
	unaligned_encoding_wrap() = delete;
	unaligned_encoding_wrap( const IntegerSequence< RangeT >& container ) : m_container(container) { }
	~unaligned_encoding_wrap() = default;

	const IntegerSequence< RangeT >& m_container;
};

template< typename RangeT >
unaligned_encoding_wrap<RangeT> unaenc( const IntegerSequence< RangeT >& container ) { return unaligned_encoding_wrap<RangeT>(container); }

template< typename RangeT >
std::ostream& operator<< ( std::ostream& os, const unaligned_encoding_wrap<RangeT>& wrapper )
{
	return encode_range_unaligned( os, wrapper.m_container );
}

// gather
template< typename IndexT, typename RealT >
RealT gather( const std::vector<IndexT>& a, const RealT *weights )
{
	using std::cbegin; using std::cend;

	return std::accumulate( cbegin(a), cend(a), RealT(0), [&weights]( RealT sum, auto offset ) { return sum + *(weights+offset); } );
}

// intersection
template< typename IndexT >
std::vector<IndexT> set_intersection( const std::vector<IndexT>& a, const std::vector<IndexT>& b )
{
	using index_t = IndexT;
	std::vector<index_t> isect;
	isect.reserve( std::min(a.size(), b.size()) ); // can't need more than this, but could be less

	const auto isect_end = std::set_intersection(
			cbegin(a), cend(a),
			cbegin(b), cend(b),
			std::back_inserter(isect)
	);

	isect.shrink_to_fit();

	return isect;
}

template< typename IndexT >
std::vector<IndexT> operator&( const std::vector<IndexT>& a, const std::vector<IndexT>& b ) { return set_intersection(a,b); }

// union
template< typename IndexT >
std::vector<IndexT> set_union( const std::vector<IndexT>& a, const std::vector<IndexT>& b )
{
	using index_t = IndexT;
	std::vector<index_t> theunion;
	theunion.reserve( std::max(a.size(), b.size()) ); // will need at least this much, maybe more

	const auto theunion_end = std::set_union(
			cbegin(a), cend(a),
			cbegin(b), cend(b),
			std::back_inserter(theunion)
	);

	theunion.shrink_to_fit();

	return theunion;
}

template< typename IndexT >
std::vector<IndexT> operator|( const std::vector<IndexT>& a, const std::vector<IndexT>& b ) { return set_union(a,b); }

template< typename IntegerT >
constexpr IntegerT create_range_flag_mask() { return IntegerT(1) << std::numeric_limits<IntegerT>::digits-1; } // set the most significant bit

// complement

// split into intersection and complement

} // namespace apegrunt

#endif // APEGRUNT_INTEGERSEQUENCE_OPERATIONS_HPP
