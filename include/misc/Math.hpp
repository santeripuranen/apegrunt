/** @file Math.hpp
 
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
#ifndef APEGRUNT_MATH_HPP
#define APEGRUNT_MATH_HPP

#include <vector>
#include <cmath>
#include <numeric> // for std::accumulate

namespace apegrunt {

/*
// integer pow
int ipow(int base, int exp)
{
    int result = 1;
    while (exp)
    {
        if (exp & 1) { result *= base; }
        exp >>= 1;
        base *= base;
    }

    return result;
}
*/
template< typename IntegerT >
IntegerT ipow(IntegerT base, int exp)
{
	IntegerT result = 1;
    while (exp)
    {
        if (exp & 1) { result *= base; }
        exp >>= 1;
        base *= base;
    }

    return result;
}

template< typename RealT >
std::size_t norm2( const std::vector< RealT >& v )
{
	using real_t = RealT;
	using std::pow;
	using std::cbegin;
	using std::cend;

	return std::accumulate( cbegin(v), cend(v), real_t(0), [=]( real_t sum, real_t x ) { return sum += x*x; } );
}

template< typename RealT >
RealT norm( const std::vector< std::size_t >& v )
{
	using real_t = RealT;
	using std::sqrt;
	auto n2 = norm2(v);
	return ( n2 != 0.0 ? sqrt( real_t(n2) ) : 0 );
}

} // namespace apegrunt

#endif // APEGRUNT_MATH_HPP
