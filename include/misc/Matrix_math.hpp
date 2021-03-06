/** @file Matrix_math.hpp
 
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
#ifndef APEGRUNT_MATRIX_MATH_HPP
#define APEGRUNT_MATRIX_MATH_HPP

#include <string>
#include <functional> // for std::hash
#include <memory> // for std::unique_ptr and std::make_unique
#include <algorithm> // for std::min, std::max
#include <array>
#include <cmath>

#include "misc/Array_view.h"

namespace apegrunt {

template< typename RealT, std::size_t N >
std::ostream& operator<< ( std::ostream& os, const std::array< std::array<RealT,N>, N >& a )
{
	for( std::size_t i = 0; i < N; ++i )
	{
		for( std::size_t j = 0; j < N; ++j )
		{
			if( i == 0 && j == 0 ) { os << a[i][j]; }
			else
			{
				os << " " << a[i][j];
			}
		}
	}
	return os;
}

template< typename RealT, std::size_t N >
std::array< std::array<RealT,N>, N > repmat( const std::array<RealT,N>& a, bool is_column_vector=false )
{
	std::array< std::array<RealT,N>, N > b;
	if( is_column_vector )
	{
		for( std::size_t i = 0; i < N; ++i )
		{
			for( std::size_t j = 0; j < N; ++j )
			{
				b[i][j] = a[i];
			}
		}
	}
	else
	{
		for( std::size_t i = 0; i < N; ++i )
		{
			for( std::size_t j = 0; j < N; ++j )
			{
				b[i][j] = a[j];
			}
		}
	}
	return b;
}

template< typename RealT, std::size_t N >
std::array< std::array<RealT,N>, N > transpose( const std::array< std::array<RealT,N>, N >& a )
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
RealT frobenius_norm( const std::array< std::array<RealT,N>, N >& a )
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
RealT frobenius_norm( const std::array< std::array<RealT,N>, N >& a, std::size_t exclude )
{
	using std::pow;

	RealT sum{0.0};

	if( exclude == N-1 )
	{
		for( std::size_t i = 0; i < N-1; ++i )
		{
			for( std::size_t j = 0; j < N-1; ++j )
			{
				sum += pow( a[i][j], 2 );
				//if( i != j ) { std::cout << std::scientific << a[i][j] << std::endl; }
			}
		}
	}

	else if( exclude == 0 )
	{
		for( std::size_t i = 1; i < N; ++i )
		{
			for( std::size_t j = 1; j < N; ++j )
			{
				sum += pow( a[i][j], 2 );
				//if( i != j ) { std::cout << std::scientific << a[i][j] << std::endl; }
			}
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
				sum += pow( a[i][j], 2 );
				//if( i != j ) { std::cout << std::scientific << a[i][j] << std::endl; }
			}
		}
	}

	return std::sqrt(sum);
}

template< typename RealT, typename RealT2, std::size_t N >
std::array< std::array<RealT,N>, N > operator*( const std::array< std::array<RealT,N>, N >& a, RealT2 b )
{
	std::array< std::array<RealT,N>, N > c;

	for( std::size_t i = 0; i < N; ++i )
	{
		for( std::size_t j = 0; j < N; ++j )
		{
			c[i][j] = a[i][j] * RealT(b);
		}
	}
	return c;
}

template< typename RealT, std::size_t N >
std::array< std::array<RealT,N>, N > operator+( const std::array< std::array<RealT,N>, N >& a, RealT b )
{
	std::array< std::array<RealT,N>, N > c;

	for( std::size_t i = 0; i < N; ++i )
	{
		for( std::size_t j = 0; j < N; ++j )
		{
			c[i][j] = a[i][j] + b;
		}
	}
	return c;
}

template< typename RealT, std::size_t N >
std::array< std::array<RealT,N>, N > operator+( const std::array< std::array<RealT,N>, N >& a, const std::array< std::array<RealT,N>, N >& b )
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
std::array< std::array<RealT,N>, N > operator-( const std::array< std::array<RealT,N>, N >& a, const std::array< std::array<RealT,N>, N >& b )
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
std::array<RealT,N> operator/( std::array<RealT,N>& a, RealT D )
{
	std::array<RealT,N> aret{0.0};
	for( std::size_t i=0; i<N; ++i ) { aret[i] = a[i] / D; }
	return aret;
}

template< typename RealT, std::size_t N >
RealT mean( const std::array<RealT,N>& a )
{
	RealT am{0.0};
	for( std::size_t i=0; i<N; ++i ) { am += a[i]; }
	return am/RealT(N);
}

template< typename RealT, std::size_t N >
std::array<RealT,N> mean( const std::array< std::array<RealT,N>, N >& a )
{
	std::array<RealT,N> am{0.0};
	for( std::size_t i=0; i<N; ++i )
	{
		for( std::size_t j=0; j<N; ++j ) { am[i] = am[i] + a[i][j]; }
		am[i] = am[i] / RealT(N);
	}
	return am;
}

/*
template< typename RealT, std::size_t N >
RealT min( const std::array<RealT,N>& a )
{
	using std::min;
	RealT am{a[0]};
	for( std::size_t i=1; i<N; ++i ) { am = min(am,a[i]); }
	return am;
}

template< typename RealT, std::size_t N >
RealT max( const std::array<RealT,N>& a )
{
	using std::max;
	RealT am{a[0]};
	for( std::size_t i=1; i<N; ++i ) { am = max(am,a[i]); }
	return am;
}
*/
} // namespace apegrunt

#endif // APEGRUNT_MATRIX_MATH_HPP
