/** @file Vector_operations.hpp
 
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
#ifndef APEGRUNT_VECTOR_OPERATIONS_HPP
#define APEGRUNT_VECTOR_OPERATIONS_HPP

#include <type_traits> // for std::enable_if
#include <cmath>

#include "Vector_forward.h"
#include "Vector_interface.hpp"

namespace apegrunt {

template< typename RealT, std::size_t Capacity, bool View >
inline RealT sum( const Vector<RealT,Capacity,View>& v )
{
	RealT s = 0;
	for( std::size_t i=0; i < Capacity; ++i ) { s += v[i]; }
	return s;
}

template< typename RealT, std::size_t Capacity, bool View, typename IntegerT >
inline RealT mask_sum( const Vector<RealT,Capacity,View>& v, IntegerT mask )
{
	RealT s = 0;
	//const std::bitset<std::numeric_limits<IntegerT>::digits> m(mask);
	//for( std::size_t i=0; i < Capacity; ++i ) { if( m[i] ) { s += v[i]; } }
	for( std::size_t i=0; i < Capacity; ++i ) { if( mask & (1<<i) ) { s += v[i]; } }
	return s;
}


template< typename RealT, std::size_t Capacity, bool View >
inline Vector<RealT,Capacity> log( const Vector<RealT,Capacity,View>& x )
{
	using std::log;
	Vector<RealT,Capacity> v;
	for( std::size_t i=0; i < Capacity; ++i ) { v[i] = log( x[i] ); }
	return v;
}

template< typename RealT, std::size_t Capacity, bool View >
inline Vector<RealT,Capacity> xaddlogx( const Vector<RealT,Capacity,View>& x )
{
	return x + log(x);
}

template< typename RealT >
inline RealT xlogx( RealT x )
{
	using std::log;
	return x * log(x);
}

template< typename RealT, std::size_t Capacity, bool View >
inline Vector<RealT,Capacity> xlogx( const Vector<RealT,Capacity,View>& x )
{
	return x * log(x);
}

template< typename RealT, std::size_t Capacity, bool View, typename IntegerT >
inline Vector<RealT,Capacity> mask_xaddlogx( const Vector<RealT,Capacity,View>& x, IntegerT mask )
{
	using std::log;
	Vector<RealT,Capacity> v;
	//const std::bitset<std::numeric_limits<IntegerT>::digits> m(mask);
	//for( std::size_t i=0; i < Capacity; ++i ) { v[i] = m[i] ? x[i] + log( x[i] ) : RealT(0); }
	for( std::size_t i=0; i < Capacity; ++i ) { v[i] = ( mask & (1<<i) ) ? x[i] + log( x[i] ) : RealT(0); }
	return v;
}

template< typename RealT, std::size_t Capacity, bool View, typename IntegerT >
inline Vector<RealT,Capacity> mask_xlogx( const Vector<RealT,Capacity,View>& x, IntegerT mask )
{
	using std::log;
	Vector<RealT,Capacity> v;
	//const std::bitset<std::numeric_limits<IntegerT>::digits> m(mask);
	//for( std::size_t i=0; i < Capacity; ++i ) { v[i] = m[i] ? x[i] * log( x[i] ) : RealT(0); }
	for( std::size_t i=0; i < Capacity; ++i ) { v[i] = ( mask & (1<<i) ) ? x[i] * log( x[i] ) : RealT(0); }
	return v;
}

template< typename RealT, std::size_t Capacity, bool View, typename IntegerT >
inline RealT mask_sum_xlogx( const Vector<RealT,Capacity,View>& x, IntegerT mask )
{
	using std::log;
	RealT s(0);
	for( std::size_t i=0; i < Capacity; ++i ) { if( mask & (1<<i) ) { s += x[i] * log( x[i] ); } }
	return s;
}

template< typename RealT, std::size_t Capacity, bool View >
inline RealT guard_sum_xlogx( const Vector<RealT,Capacity,View>& x )
{
	using std::log;
	RealT s(0);
	for( std::size_t i=0; i < Capacity; ++i ) { if( x[i] != 0 ) { s += x[i] * log( x[i] ); } }
	return s;
}

// basic binary arithmetic operators
template< typename RealT, std::size_t Capacity, bool aView, bool bView >
inline Vector<RealT,Capacity> operator+( const Vector<RealT,Capacity,aView>& lhs, const Vector<RealT,Capacity,bView>& rhs )
{
	Vector<RealT,Capacity> v;
	for( std::size_t i=0; i < Capacity; ++i ) { v[i] = lhs[i] + rhs[i]; }
	return v;
}

template< typename RealT, std::size_t Capacity, bool aView, bool bView >
inline Vector<RealT,Capacity> operator-( const Vector<RealT,Capacity,aView>& lhs, const Vector<RealT,Capacity,bView>& rhs )
{
	Vector<RealT,Capacity> v;
	for( std::size_t i=0; i < Capacity; ++i ) { v[i] = lhs[i] - rhs[i]; }
	return v;
}

template< typename RealT, std::size_t Capacity, bool aView, bool bView >
inline Vector<RealT,Capacity> operator*( const Vector<RealT,Capacity,aView>& lhs, const Vector<RealT,Capacity,bView>& rhs )
{
	Vector<RealT,Capacity> v;
	for( std::size_t i=0; i < Capacity; ++i ) { v[i] = lhs[i] * rhs[i]; }
	return v;
}

template< typename RealT, std::size_t Capacity, bool aView, bool bView >
inline Vector<RealT,Capacity> operator/( const Vector<RealT,Capacity,aView>& lhs, const Vector<RealT,Capacity,bView>& rhs )
{
	Vector<RealT,Capacity> v;
	for( std::size_t i=0; i < Capacity; ++i ) { v[i] = lhs[i] / rhs[i]; }
	return v;
}

// bitshift integer fields in lhs by the amount of steps specified by fields in rhs
template< typename RealT, std::size_t Capacity, bool View, typename StateT > // should use std::enable_if to drop from overload resolution if not std::is_integral<RealT>::value
inline Vector<RealT,Capacity> operator<< ( const Vector<RealT,Capacity,View>& lhs, const State_block<StateT,Capacity>& rhs )
{
	Vector<RealT,Capacity> v;
	for( std::size_t i = 0; i < Capacity; ++i )
	{
		v[i] = lhs[i] << std::size_t( rhs[i] );
	}
	return v;
}

// bitwise AND for array elements
template< typename RealT, std::size_t Capacity, bool View, bool View2 >
inline Vector<RealT,Capacity> operator& ( const Vector<RealT,Capacity,View>& lhs, const Vector<RealT,Capacity,View2>& rhs )
{
	Vector<RealT,Capacity> v;
	for( std::size_t i = 0; i < Capacity; ++i )
	{
		v[i] = lhs[i] & rhs[i];
	}
	return v;
}

// bitwise OR for array elements
template< typename RealT, std::size_t Capacity, bool View, bool View2 >
inline Vector<RealT,Capacity> operator| ( const Vector<RealT,Capacity,View>& lhs, const Vector<RealT,Capacity,View2>& rhs )
{
	Vector<RealT,Capacity> v;
	for( std::size_t i = 0; i < Capacity; ++i )
	{
		v[i] = lhs[i] | rhs[i];
	}
	return v;
}

template< typename RealT, std::size_t Capacity, bool View >
inline Vector<RealT,Capacity> popcnt_per_element( const Vector<RealT,Capacity,View>& a )
{
	Vector<RealT,Capacity> v;
	for( std::size_t i = 0; i < Capacity; ++i )
	{
		// OK, this is a lazy solution, but it works
		v[i] = std::bitset<std::numeric_limits<RealT>::digits>(a[i]).count();
	}
	return v;
}

template< typename RealT, std::size_t Capacity, bool View >
inline std::size_t popcnt( const Vector<RealT,Capacity,View>& a )
{
	std::size_t s;
	for( std::size_t i = 0; i < Capacity; ++i )
	{
		// OK, this is a lazy _and_ slow solution, but it works
		s += std::bitset<std::numeric_limits<RealT>::digits>(a[i]).count();
	}
	return s;
}

template< typename RealT >
inline std::size_t popcnt( RealT a )
{
	// OK, this is a lazy _and_ slow solution, but it works
	return std::bitset<std::numeric_limits<RealT>::digits>(a).count();
}

template< typename RealT, std::size_t Capacity, bool View >
inline RealT distance2( const Vector<RealT,Capacity,View>& a, const Vector<RealT,Capacity,View>& b )
{
	RealT d2(0);
	//std::cout << "distance2[";
	for( std::size_t i=0; i < Capacity; ++i )
	{
		//std::cout << " " << std::pow( a[i]-b[i], 2 );
		d2 += std::pow( a[i]-b[i], 2 );
	}
	//std::cout << " ]=" << d2 << std::endl;
	return d2;
}

template< typename RealT, std::size_t Capacity, bool View >
inline RealT distance( const Vector<RealT,Capacity,View>& a, const Vector<RealT,Capacity,View>& b )
{
	return std::sqrt( distance2(a,b) );
}

#ifndef NO_INTRINSICS

#ifdef __SSE2__
inline Vector<uint8_t,16,false> operator& ( Vector<uint8_t,16,false> lhs, Vector<uint8_t,16,false> rhs )
{
	return Vector<uint8_t,16,false>( _mm_and_si128( lhs(), rhs() ) );
}

inline Vector<uint8_t,16,false> operator| ( Vector<uint8_t,16,false> lhs, Vector<uint8_t,16,false> rhs )
{
	return Vector<uint8_t,16,false>( _mm_or_si128( lhs(), rhs() ) );
}
#endif // __SSE2__


#ifdef __AVX__

inline double sum( __m256d a )
{
	using vec_t = Vector<double,4,false>;
	using simd_t = __m256d;
	const simd_t b( _mm256_permute2f128_pd( a, a, 0b0001 ) );
	const simd_t c( _mm256_add_pd( a, b ) );
	const vec_t d( _mm256_hadd_pd( c, c ) );

	return typename vec_t::element_t( d.m_elem[0] );
}

inline double sum( const Vector<double,4,false>& v )
{
	return sum( v() ); // forward to sum( __m256d a )
}
// /*
inline double sum( const Vector<double,5,false>& v )
{
	return sum( Vector<double,4,false>( v.data() ) ) + v[4];
}
// */
template< uint Exponent >
inline __m256d pow( __m256d x ) { __m256d result = x; for( std::size_t i=1; i<Exponent; ++i ) { result = result * x; } return result; }
template<>
inline __m256d pow<2>( __m256d x ) { return x*x; }

#endif // __AVX__

#endif // #ifndef NO_INTRINSICS

} // namespace apegrunt

#endif // APEGRUNT_VECTOR_OPERATIONS_HPP
