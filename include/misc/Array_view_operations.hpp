/** @file Array_view_operations.hpp
 
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
#ifndef APEGRUNT_ARRAY_VIEW_OPERATIONS_HPP
#define APEGRUNT_ARRAY_VIEW_OPERATIONS_HPP

#include <cstring> // for std::memcpy
#include <array>
#include <cmath>

#include "misc/Array_view_forward.h"
#include "misc/Array_view.hpp"

namespace apegrunt {

template< typename RealT, std::size_t Extent >
void copy( const Array_view< Array_view< RealT, Extent >, Extent >& source, Array_view< Array_view< RealT, Extent >, Extent >& dest, bool transpose=false )
{
	if( transpose )
	{
		for( std::size_t i=0; i < Extent; ++i )
		{
			for( std::size_t j=0; j < Extent; ++j )
			{
				dest[i][j] = source[j][i];
			}
		}
	}
	else
	{
		std::memcpy( dest.data(), source.data(), Extent*Extent );
/*
		for( std::size_t i=0; i < Extent; ++i )
		{
			for( std::size_t j=0; j < Extent; ++j )
			{
				dest[i][j] = source[i][j];
			}
		}
*/
	}
}

template< typename RealT, std::size_t N >
std::array< std::array<RealT,N>, N > transpose( const Array_view< Array_view<RealT,N>, N >& a )
{
	std::array< std::array<RealT,N>, N > b;
	for( std::size_t i = 0; i < N; ++i )
	{
		for( std::size_t j = 0; j < N; ++j )
		{
			b[j][i] = a[i][j];
		}
	}
	return b;
}

template< typename RealT, std::size_t N >
RealT frobenius_norm( const Array_view< Array_view<RealT,N>, N >& a )
{
	using std::pow; using std::sqrt;

	RealT sum{0.0};
	for( std::size_t i = 0; i < N; ++i )
	{
		for( std::size_t j = 0; j < N; ++j )
		{
			sum += pow( a[i][j], 2 );
		}
	}
	return sqrt(sum);
}

template< typename RealT, std::size_t N >
RealT frobenius_norm( const Array_view< Array_view<RealT,N>, N >& aview, std::size_t exclude )
{
	using std::pow; using std::sqrt;

	RealT sum{0.0};

	if( exclude == N-1 )
	{
		for( std::size_t i=0; i<N-1; ++i )
		{
			for( std::size_t j=0; j<N-1; ++j ) { sum += pow( aview[i][j], 2 ); }
		}
	}
	else if( exclude == 0 )
	{
		for( std::size_t i = 1; i < N; ++i )
		{
			for( std::size_t j = 1; j < N; ++j ) { sum += pow( aview[i][j], 2 ); }
		}
	}

	else
	{
		for( std::size_t i = 0; i < N; ++i )
		{
			if( exclude==i ) { continue; }
			for( std::size_t j = 0; j < N; ++j )
			{
				if( exclude==j ) { continue; }
				sum += pow( aview[i][j], 2 );
			}
		}
	}
	return sqrt(sum);
}

template< typename RealT, std::size_t N >
std::array< std::array<RealT,N>, N > operator+( const Array_view< Array_view<RealT,N>, N >& a, const Array_view< Array_view<RealT,N>, N >& b )
{
	std::array< std::array<RealT,N>, N > c;

	for( std::size_t i = 0; i < N; ++i )
	{
		for( std::size_t j = 0; j < N; ++j )
		{
			c[i][j] = a[i][j] + b[i][j];
		}
	}
	return c;
}

template< typename RealT, std::size_t N >
std::array< std::array<RealT,N>, N > operator-( const std::array< std::array<RealT,N>, N >& a, const Array_view< Array_view<RealT,N>, N >& b )
{
	std::array< std::array<RealT,N>, N > c;

	for( std::size_t i = 0; i < N; ++i )
	{
		for( std::size_t j = 0; j < N; ++j )
		{
			c[i][j] = a[i][j] - b[i][j];
		}
	}
	return c;
}

template< typename RealT, std::size_t N >
std::array< std::array<RealT,N>, N > operator-( const Array_view< Array_view<RealT,N>, N >& a, const std::array< std::array<RealT,N>, N >& b )
{
	std::array< std::array<RealT,N>, N > c;

	for( std::size_t i = 0; i < N; ++i )
	{
		for( std::size_t j = 0; j < N; ++j )
		{
			c[i][j] = a[i][j] - b[i][j];
		}
	}
	return c;
}

template< typename RealT, std::size_t N >
std::array<RealT,N> operator+( std::array<RealT,N>& a, Array_view<RealT,N>&& aview )
{
	std::array<RealT,N> aret{0.0};
	for( std::size_t i=0; i<N; ++i ) { aret[i] = a[i] + aview[i]; }
	return aret;
}

template< typename RealT, std::size_t N >
RealT mean( const Array_view<RealT,N>& a )
{
	RealT am{0.0};
	for( std::size_t i=0; i<N; ++i ) { am += a[i]; }
	return am/RealT(N);
}

template< typename RealT, std::size_t N >
std::array<RealT,N> mean( const Array_view< Array_view<RealT,N>, N >& aview )
{
	std::array<RealT,N> am{0.0};
	for( std::size_t i=0; i<N; ++i )
	{
		for( std::size_t j=0; j<N; ++j ) { am[i] = am[i] + aview[i][j]; }
		am[i] = am[i] / RealT(N);
	}
	return am;
}

/*
template< typename RealT, std::size_t N >
RealT min( const Array_view<RealT,N>& a )
{
	using std::min;
	RealT am{a[0]};
	for( std::size_t i=1; i<N; ++i ) { am = min(am,a[i]); }
	return am;
}

template< typename RealT, std::size_t N >
std::array<RealT,N> min( const Array_view< Array_view<RealT,N>, N >& a )
{
	using std::min;
	std::array<RealT,N> am{0.0};
	for( std::size_t i=0; i<N; ++i )
	{
		am[i] = min(a[i]);
	}
	return am;
}

template< typename RealT, std::size_t N >
RealT max( const Array_view<RealT,N>& a )
{
	using std::max;
	RealT am{a[0]};
	for( std::size_t i=1; i<N; ++i ) { am = max(am,a[i]); }
	return am;
}

template< typename RealT, std::size_t N >
std::array<RealT,N> abs( const Array_view<RealT,N>& a )
{
	using std::abs;
	std::array<RealT,N> am{0};
	for( std::size_t i=1; i<N; ++i ) { am[i] = abs(a[i]); }
	return am;
}

template< typename RealT, std::size_t N >
std::array<RealT,N> max( const Array_view< Array_view<RealT,N>, N >& a )
{
	using std::max;
	using apegrunt::max;
	std::array<RealT,N> am{0.0};
	for( std::size_t i=0; i<N; ++i )
	{
		am[i] = max(a[i]);
	}
	return am;
}

template< typename RealT, std::size_t N >
RealT max_element( const Array_view< Array_view<RealT,N>, N >& a )
{
	using std::max;
	using apegrunt::max;
	RealT amax = a[0][0]; // first element
	for( std::size_t i=0; i<N; ++i )
	{
		amax = max( max(a[i]), amax );
	}
	return amax;
}

template< typename RealT, std::size_t N >
RealT max_abs_element( const Array_view< Array_view<RealT,N>, N >& a )
{
	using std::abs;
	using std::max;
	using apegrunt::max;
	RealT amax = abs( a[0][0] ); // first element
	for( std::size_t i=0; i<N; ++i )
	{
		amax = max( max( abs(a[i]) ), amax );
	}
	return amax;
}
*/

/*
template<>
void copy( const Array_view< Array_view< double, 4 >, 4 >& source, Array_view< Array_view< double, 4 >, 4 >& dest, bool transpose=false )
{
	if( transpose )
	{
		for( std::size_t i=0; i < 4; ++i )
		{
			for( std::size_t j=0; j < 4; ++j )
			{
				dest[i][j] = source[j][i];
			}
		}
	}
	else
	{

	}
}
*/
} // namespace apegrunt

#endif // APEGRUNT_ARRAY_VIEW_OPERATIONS_HPP
