/** @file accumulators_utility.hpp
 
	Copyright (c) 2016-2019 Santeri Puranen.

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
#ifndef APEGRUNT_ACCUMUALTORS_UTILITY_HPP
#define APEGRUNT_ACCUMUALTORS_UTILITY_HPP

#include <string>
#include <functional> // for std::hash
#include <memory> // for std::unique_ptr and std::make_unique
#include <algorithm> // for std::min, std::max
#include <array>
#include <cmath>

#include <boost/accumulators/accumulators.hpp>

namespace apegrunt {

template< typename RealT, std::size_t N, typename ...Args >
boost::accumulators::accumulator_set< RealT, Args... >& operator<< ( boost::accumulators::accumulator_set< RealT, Args... >& acc, const std::array< std::array<RealT,N>, N >& a )
{
	for( std::size_t i = 0; i < N; ++i )
	{
		for( std::size_t j = 0; j < N; ++j )
		{
			acc( a[i][j] );
		}
	}
	return acc;
}

template< typename RealT, std::size_t N, typename ...Args >
boost::accumulators::accumulator_set< RealT, Args... >& operator<< ( boost::accumulators::accumulator_set< RealT, Args... >& acc, const std::array< RealT, N >& a )
{
	for( std::size_t i = 0; i < N; ++i )
	{
		if( 0.0 != a[i] ) { acc( a[i] ); } // this might not be a good thing in general
	}
	return acc;
}

template< typename RealT, std::size_t N, typename ...Args >
boost::accumulators::accumulator_set< RealT, Args... >& operator<< ( boost::accumulators::accumulator_set< RealT, Args... >& acc, const std::vector< std::array<RealT,N> >& v )
{
	for( const auto& a: v )
	{
		acc << a;
	}
	return acc;
}
} // namespace apegrunt

#endif // APEGRUNT_ACCUMUALTORS_UTILITY_HPP
