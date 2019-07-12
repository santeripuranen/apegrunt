/** @file IntegerSequence_interface.hpp
 
	Copyright (c) 2016-2019 Santeri Puranen. All rights reserved.
 
	By installing, copying or otherwise using the attached
	material ("product" or "software") you acknowledge and
	agree that the attached	material contains proprietary
	information of the copyright holder(s). Any use of the
	material is prohibited except as expressly agreed between
	the copyright holder(s) and the recipient.
 
	THIS PRODUCT ("SOFTWARE") IS PROVIDED "AS IS", WITHOUT WARRANTY
	OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO
	THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
	PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT SHALL THE
	COPYRIGHT HOLDER(S) BE LIABLE FOR ANY DAMAGES OR OTHER LIABILITY,
	WHETHER IN CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR
	IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
	THE SOFTWARE.

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

// Generic outstream operator for IntegerSequences
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

// outstream operator for specialization for vectors of integers
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

// set intersection for vectors of integers
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

template< typename IntegerT >
constexpr IntegerT create_range_flag_mask() { return IntegerT(1) << std::numeric_limits<IntegerT>::digits-1; } // set the most significant bit

// complement

// union

// split into intersection and complement

} // namespace apegrunt

#endif // APEGRUNT_INTEGERSEQUENCE_OPERATIONS_HPP
