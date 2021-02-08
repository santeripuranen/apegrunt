/** @file Apegrunt_utility.hpp

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
#ifndef APEGRUNT_UTILITY_HPP
#define APEGRUNT_UTILITY_HPP

#include <iostream>
#include <iomanip>
#include <cmath>

#include <boost/tuple/tuple.hpp>
#include <boost/range/iterator_range_core.hpp>
#include <boost/iterator/zip_iterator.hpp>

namespace apegrunt {

// Helpers for iterables that are passed in shared_ptrs
template< typename T > auto begin( std::shared_ptr<T>& sptr ) { return begin(*sptr); }
template< typename T > auto end( std::shared_ptr<T>& sptr ) { return end(*sptr); }

template< typename T > auto begin( const std::shared_ptr<T>& sptr ) { return begin(*sptr); }
template< typename T > auto end( const std::shared_ptr<T>& sptr ) { return end(*sptr); }

template< typename T > auto cbegin( const std::shared_ptr<T>& sptr ) { return cbegin(*sptr); }
template< typename T > auto cend( const std::shared_ptr<T>& sptr ) { return cend(*sptr); }

template< typename... Containers >
auto zip_range( Containers&&... containers ) // universal reference; will bind to anything
	-> decltype(
		boost::make_iterator_range(
			boost::make_zip_iterator( boost::make_tuple( begin(containers) ... ) ),
			boost::make_zip_iterator( boost::make_tuple( end(containers) ... ) )
		)
	)
{
	return { boost::make_zip_iterator( boost::make_tuple( begin(containers) ... ) ),
		     boost::make_zip_iterator( boost::make_tuple( end(containers) ... ) )
	};
}

// To use, inherit from extend_comparison_operators and
// define these public operators in your struct/class:
// inline bool T::operator==( const T& rhs ) const
// inline bool operator<( const my_type& rhs ) const
struct extend_comparison_operators { using enable_extended_comparison_operators = std::true_type; };

template< typename T, typename std::enable_if< std::is_same< typename T::enable_extended_comparison_operators,std::true_type >::type >::type >
inline bool operator!=( const T& lhs, const T& rhs ) { return !(lhs == rhs); }

template< typename T, typename std::enable_if< std::is_same< typename T::enable_extended_comparison_operators,std::true_type >::type >::type >
inline bool operator> ( const T& lhs, const T& rhs ) { return (rhs < lhs); }

template< typename T, typename std::enable_if< std::is_same< typename T::enable_extended_comparison_operators,std::true_type >::type >::type >
inline bool operator<=( const T& lhs, const T& rhs ) { return !(lhs > rhs); }

template< typename T, typename std::enable_if< std::is_same< typename T::enable_extended_comparison_operators,std::true_type >::type >::type >
inline bool operator>=( const T& lhs, const T& rhs ) { return !(lhs < rhs); }

// Binary search in an iterator-defined range; returns an iterator to the matching element, else last.
template< typename IteratorT >
inline IteratorT binary_search( IteratorT first, IteratorT last, const typename IteratorT::value_type& key )
{
    for( auto a = first, b = last, mid = first + std::distance(first,last) / 2; a != b; mid = a + std::distance(a,b) / 2 )
    {
    	if( *mid == key ) { return mid; }
        else if( *mid < key ) { a = mid + 1; }
        else { b = mid; }
    }
    return last;
}

struct my_div_t { uint64_t quot; uint64_t rem; };
my_div_t my_div( uint64_t n, uint64_t div ) { my_div_t result{}; result.quot=n/div; result.rem=n%div; return result; }

struct memory_string
{
	memory_string( uint64_t mem_in_bytes ) : m_mem_in_bytes(mem_in_bytes) { }
	~memory_string() { }

	std::ostream& operator()( std::ostream& os ) const
	{
		my_div_t result{};
		result.quot = m_mem_in_bytes;
		std::size_t n = 0;

		//std::cout << m_elapsed_time << std::endl;

		while( result.quot > 1024 && n < 4 )
		{
			++n;
			result = my_div(result.quot,1000);
			//std::cout << "rem=" << result.rem << " quot=" << result.quot << " n=" << n << "\n";
		};

		std::string unit;
		switch(n)
		{
			case 0: unit = "B"; break;
			case 1: unit = "KiB"; break;
			case 2: unit = "MiB"; break;
			case 3: unit = "GiB"; break;
			case 4: unit = "TiB"; break;
		}

		os << std::fixed << std::setprecision(2) << double(m_mem_in_bytes)/double(std::pow(1024,n)) << unit;

		return os;
	}

	uint64_t m_mem_in_bytes;
};

std::ostream& operator<< ( std::ostream& os, const memory_string& mem )
{
	return mem(os);
}

} // namespace apegrunt

#endif // APEGRUNT_UTILITY_HPP
