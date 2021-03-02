/** @file IntegerSequence_interface.hpp
 
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

#ifndef APEGRUNT_INTEGERSEQUENCE_OPERATIONS_HPP
#define APEGRUNT_INTEGERSEQUENCE_OPERATIONS_HPP

#include <vector>
#include <iterator> // for std::back_inserter
#include <ostream>
#include <limits>
#include <algorithm> // for std::set_intersection
#include <functional> // for std::reference_wrapper

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

template< typename IndexT >
std::ostream& operator<< ( std::ostream& os, const std::vector< IndexT >& container )
{
	using std::cbegin; using std::cend;

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

template< typename RangeT >
std::ostream& operator<< ( std::ostream& os, const IntegerSequence< RangeT >& container )
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

template< typename RangeT >
struct unaligned_encoding
{
	unaligned_encoding() = delete;
	unaligned_encoding( const IntegerSequence< RangeT >& container ) : m_container(container) { }
	~unaligned_encoding() = default;

	const IntegerSequence< RangeT >& m_container;
};

template< typename RangeT >
unaligned_encoding<RangeT> unaenc( const IntegerSequence< RangeT >& container ) { return unaligned_encoding<RangeT>(container); }

template< typename RangeT >
std::ostream& operator<< ( std::ostream& os, const unaligned_encoding<RangeT>& wrapper )
{
	using std::cbegin; using std::cend;
	using apegrunt::cbegin; using apegrunt::cend;

	auto currpos = cbegin( wrapper.m_container );
	auto begpos = currpos;
	auto oldpos = currpos;
	const auto endpos = cend( wrapper.m_container );

	os << std::hex << *currpos;
	while( currpos != endpos )
	{
		++currpos;
		if( *currpos > (*oldpos)+1 )
		{
			if( oldpos != begpos )
			{
				os << "-" << *oldpos; // os.flush();
			}
			os << "," << *currpos; // os.flush();
			begpos = currpos;
		}
		oldpos = currpos;
	}
	return os;
}

// gather
template< typename IndexT, typename RealT >
RealT gather( const std::vector<IndexT>& a, const RealT *weights )
{
	using std::cbegin; using std::cend;

	return std::accumulate( cbegin(a), cend(a), RealT(0), [&weights]( RealT sum, auto offset ) { return sum + *(weights+offset); } );
}

// gather
template< typename RealT >
struct gatherer
{
	using real_t = RealT;

	gatherer( const real_t *weights ) : w(weights), sum(0) { }

	template< typename ElementT >
	void operator()( const ElementT& element ) { sum += ElementT::gather(w,element); }

	const real_t *w;
	real_t sum;
};

// intersection
template< typename IndexT >
std::vector<IndexT> set_intersection( const std::vector<IndexT>& a, const std::vector<IndexT>& b )
{
	using index_t = IndexT;
	std::vector<index_t> isect;
	isect.reserve( std::min(a.size(), b.size()) ); // can't need more than this, but might be less

	const auto isect_end = std::set_intersection(
			cbegin(a), cend(a),
			cbegin(b), cend(b),
			std::back_inserter(isect)
	);

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

	return theunion;
}

template< typename IndexT >
std::vector<IndexT> operator|( const std::vector<IndexT>& a, const std::vector<IndexT>& b ) { return set_union(a,b); }

template< typename IndexT >
std::vector<IndexT> set_union(
		const std::vector< std::reference_wrapper< const std::vector<IndexT> > >& sets
	)
{
	using container_t = std::vector<IndexT>;
	using std::cbegin; using std::cend;
	using std::begin; using std::end;

	std::size_t reserve = std::accumulate( cbegin(sets), cend(sets), 0, []( auto value, const auto& set ){ return value+set.get().size(); } );
	container_t the_union; the_union.reserve(reserve);
	for( const auto& set: sets )
	{
		the_union.insert( the_union.end(), cbegin(set.get()), cend(set.get()) );
	}
	std::sort( the_union.begin(), the_union.end() );
	auto last = std::unique( the_union.begin(), the_union.end());
	the_union.erase(last, the_union.end());

	return the_union;
}

template< typename SetT >
struct lazy_set_union
{
	using my_type = lazy_set_union<SetT>;
	using wrap_t = std::reference_wrapper< const SetT >;
	lazy_set_union() : m_sets() { }
	~lazy_set_union() = default;

	inline my_type& add( const SetT& set ) {  m_sets.emplace_back( wrap_t(set) ); return *this; }
	inline my_type& operator|=( const SetT& set ) { return this->add(set); }
	inline operator SetT() const { return apegrunt::set_union(m_sets); }

	std::vector< std::reference_wrapper< const SetT > > m_sets;
};

// complement

// symmetric difference

// relative complements

// split into intersection and complement


} // namespace apegrunt

#endif // APEGRUNT_INTEGERSEQUENCE_OPERATIONS_HPP
