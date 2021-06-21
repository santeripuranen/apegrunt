/** @file Vector_interface.hpp
 
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
#ifndef APEGRUNT_VECTOR_INTERFACE_HPP
#define APEGRUNT_VECTOR_INTERFACE_HPP

#include "Vector_forward.h" // consistency check

namespace apegrunt {

template< typename RealT, std::size_t Capacity, bool View=false >
struct Vector
{
	enum { N=Capacity };
	using element_t = RealT;
	using my_type = Vector<element_t,N,false>;

	inline Vector() { for( std::size_t i=0; i < N; ++i ) { m_elem[i] = element_t(0.0); } }
 	inline Vector( Vector<RealT,Capacity>&& v ) noexcept { for( std::size_t i=0; i < N; ++i ) { m_elem[i] = v.m_elem[i]; } }
	inline Vector( const Vector<RealT,Capacity>& v ) { for( std::size_t i=0; i < N; ++i ) { m_elem[i] = v.m_elem[i]; } }
	inline Vector( element_t e ) { for( std::size_t i=0; i < N; ++i ) { m_elem[i] = e; } }
	inline Vector( element_t e, uint8_t mask ) { for( std::size_t i=0; i < N; ++i ) { m_elem[i] = bool( mask & (1<<i) ) ? e : 0; } }
	inline Vector( const element_t* const p, std::size_t stride=1 ) { for( std::size_t i=0; i < N; ++i ) { m_elem[i] = *(p+i*stride); } }

	inline my_type& operator=( Vector<RealT,Capacity>&& v ) { for( std::size_t i=0; i < N; ++i ) { m_elem[i] = v.m_elem[i]; } return *this; }
	inline my_type& operator=( const Vector<RealT,Capacity>& v ) { for( std::size_t i=0; i < N; ++i ) { m_elem[i] = v.m_elem[i]; } return *this; }

	inline element_t* data() { return m_elem; }
	inline const element_t* data() const { return m_elem; }

	inline my_type& operator()() { return *this; }
	inline const my_type& operator()() const { return *this; }

	inline constexpr element_t& operator[]( uint i ) { return m_elem[i]; }
	inline constexpr const element_t& operator[]( uint i ) const { return m_elem[i]; }

	template< bool ViewFlag >
	inline my_type& operator+=( const Vector<RealT,Capacity,ViewFlag>& v ) { for( std::size_t i=0; i < N; ++i ) { m_elem[i] += v[i]; } return *this; }
	template< bool ViewFlag >
	inline my_type& operator-=( const Vector<RealT,Capacity,ViewFlag>& v ) { for( std::size_t i=0; i < N; ++i ) { m_elem[i] -= v[i]; } return *this; }
	template< bool ViewFlag >
	inline my_type& operator*=( const Vector<RealT,Capacity,ViewFlag>& v ) { for( std::size_t i=0; i < N; ++i ) { m_elem[i] *= v[i]; } return *this; }
	template< bool ViewFlag >
	inline my_type& operator/=( const Vector<RealT,Capacity,ViewFlag>& v ) { for( std::size_t i=0; i < N; ++i ) { m_elem[i] /= v[i]; } return *this; }

	union
	{
		element_t m_elem[N];
	};
};

template< typename RealT, std::size_t Capacity, bool View >
static std::ostream& operator<< ( std::ostream& os, const Vector<RealT,Capacity,View>& v )
{
	os << "[";
	for( std::size_t i=0; i < Capacity; ++i )
	{
		os << " " << v[i];
	}
	os << " ]";
	return os;
}

template< std::size_t Capacity, bool View >
static std::ostream& operator<< ( std::ostream& os, const Vector<uint8_t,Capacity,View>& v )
{
	os << "[";
	for( std::size_t i=0; i < Capacity; ++i )
	{
		os << " " << std::size_t(v[i]); // uint8_t won't print properly without casting
	}
	os << " ]";
	return os;
}

template< std::size_t Capacity, bool View >
static std::ostream& operator<< ( std::ostream& os, const Vector<uint32_t,Capacity,View>& v )
{
	os << "[";
	for( std::size_t i=0; i < Capacity; ++i )
	{
		/*
		os << " ";
		for( std::size_t j=0; j < std::numeric_limits<uint32_t>::digits; ++j )
		{
			// let's do binary
			os << ( ( v[i] & ( uint32_t(1) << j ) ) ? "1" : "0" );
		}
		*/
		os << " " << std::size_t(v[i]);
	}
	os << " ]";
	return os;
}

template< typename RealT, std::size_t Capacity >
struct Vector<RealT,Capacity,true>
{
	enum { N=Capacity };
	using element_t = RealT;
	using my_type = Vector<element_t,N,true>;

	Vector() = delete; // Vector views always need pointer to storage
	inline constexpr Vector( element_t* const p, std::size_t stride=1 ) : m_p(p), m_stride(stride) { }
	inline constexpr Vector( my_type&& v ) : m_p( std::move(v.m_p) ), m_stride( std::move(v.m_stride) ) { }
	inline constexpr Vector( const my_type& v ) : m_p(v.m_p), m_stride(v.m_stride) { }
	// generic Vector
	inline my_type& operator=( Vector<RealT,Capacity>&& v ) { for( std::size_t i=0; i < N; ++i ) { (*this)[i] = v[i]; } return *this; }
	inline my_type& operator=( const Vector<RealT,Capacity>& v ) { for( std::size_t i=0; i < N; ++i ) { (*this)[i] = v[i]; } return *this; }
	// my_type
	inline my_type& operator=( my_type&& v ) { for( std::size_t i=0; i < N; ++i ) { (*this)[i] = v[i]; } return *this; }
	inline my_type& operator=( const my_type& v ) { for( std::size_t i=0; i < N; ++i ) { (*this)[i] = v[i]; } return *this; }

	inline element_t* data() { return m_p; }
	inline const element_t* data() const { return m_p; }

	inline my_type& operator()() { return *this; }
	inline const my_type& operator()() const { return *this; }

	inline constexpr element_t& operator[]( uint i ) { return *(m_p+i*m_stride); }
	inline constexpr const element_t& operator[]( uint i ) const { return *(m_p+i*m_stride); }

	template< bool ViewFlag >
	inline my_type& operator+=( const Vector<RealT,Capacity,ViewFlag>& v ) { for( std::size_t i=0; i < N; ++i ) { (*this)[i] += v[i]; } return *this; }
	template< bool ViewFlag >
	inline my_type& operator-=( const Vector<RealT,Capacity,ViewFlag>& v ) { for( std::size_t i=0; i < N; ++i ) { (*this)[i] -= v[i]; } return *this; }
	template< bool ViewFlag >
	inline my_type& operator*=( const Vector<RealT,Capacity,ViewFlag>& v ) { for( std::size_t i=0; i < N; ++i ) { (*this)[i] *= v[i]; } return *this; }
	template< bool ViewFlag >
	inline my_type& operator/=( const Vector<RealT,Capacity,ViewFlag>& v ) { for( std::size_t i=0; i < N; ++i ) { (*this)[i] /= v[i]; } return *this; }

	inline my_type& operator+=( RealT e ) { for( std::size_t i=0; i < N; ++i ) { (*this)[i] += e; } return *this; }
	inline my_type& operator-=( RealT e ) { for( std::size_t i=0; i < N; ++i ) { (*this)[i] -= e; } return *this; }
	inline my_type& operator*=( RealT e ) { for( std::size_t i=0; i < N; ++i ) { (*this)[i] *= e; } return *this; }
	inline my_type& operator/=( RealT e ) { for( std::size_t i=0; i < N; ++i ) { (*this)[i] /= e; } return *this; }

	element_t* const m_p;
	const std::size_t m_stride;
};

} // namespace apegrunt

#include "SIMD_intrinsics.h"

namespace apegrunt {

#ifndef NO_INTRINSICS

#ifdef __SSE2__

template<>
struct alignas(16) Vector<uint8_t,16,false>
{
	enum { N=16 };
	using element_t = uint8_t;
	using simd_t = __m128i;
	using my_type = Vector<element_t,N,false>;

	inline Vector() : m_vec( _mm_setzero_si128() ) { }
 	inline Vector( Vector<element_t,N>&& v ) noexcept : m_vec(v()) { }
	inline Vector( const Vector<element_t,N>& v ) : m_vec(v()) { }
	inline Vector( simd_t v ) : m_vec(v) { }
	inline Vector( element_t e ) : m_vec( _mm_set1_epi8(e) ) { }
	inline Vector( const element_t* const e, std::size_t stride=1 )
	{
		switch(stride)
		{
			case 1: m_vec = _mm_loadu_si128( (__m128i const*)e ); break;
			default: for( std::size_t i=0; i < N; ++i ) { m_elem[i] = *(e+i*stride); }
		}
	}

	inline my_type& operator=( Vector<element_t,N>&& v ) { m_vec = v(); return *this; }
	inline my_type& operator=( const Vector<element_t,N>& v ) { m_vec = v(); return *this; }

	inline my_type& operator=( simd_t v ) { m_vec = v; return *this; }

	inline element_t* data() { return m_elem; }
	inline const element_t* data() const { return m_elem; }

	inline simd_t operator()() { return m_vec; }
	inline simd_t operator()() const { return m_vec; }

	inline element_t& operator[]( uint i ) { return m_elem[i]; }
	inline const element_t& operator[]( uint i ) const { return m_elem[i]; }

	template< bool ViewFlag >
	inline my_type& operator+=( const Vector<element_t,N,ViewFlag>& v ) { m_vec = _mm_add_epi8( m_vec, v() ); return *this; }
	template< bool ViewFlag >
	inline my_type& operator-=( const Vector<element_t,N,ViewFlag>& v ) { m_vec = _mm_sub_epi8( m_vec, v() ); return *this; }
	template< bool ViewFlag >
	inline my_type& operator*=( const Vector<element_t,N,ViewFlag>& v ) { for( std::size_t i=0; i < N; ++i ) { (*this)[i] *= v[i]; }; return *this; }
	template< bool ViewFlag >
	inline my_type& operator/=( const Vector<element_t,N,ViewFlag>& v ) { for( std::size_t i=0; i < N; ++i ) { (*this)[i] /= v[i]; }; return *this; }
/*
	inline my_type& operator+=( const simd_t v ) { m_vec += v; return *this; }
	inline my_type& operator-=( const simd_t v ) { m_vec -= v; return *this; }
	inline my_type& operator*=( const simd_t v ) { m_vec *= v; return *this; }
	inline my_type& operator/=( const simd_t v ) { m_vec /= v; return *this; }
*/
	inline my_type& operator+=( element_t e ) { m_vec += my_type(e)(); return *this; }
	inline my_type& operator-=( element_t e ) { m_vec -= my_type(e)(); return *this; }
	inline my_type& operator*=( element_t e ) { m_vec *= my_type(e)(); return *this; }
	inline my_type& operator/=( element_t e ) { m_vec /= my_type(e)(); return *this; }

	union
	{
		element_t m_elem[N];
		simd_t m_vec;
	};
};

#endif // __SSE2__

#ifdef __AVX__

template<>
struct alignas(32) Vector<double,4,false>
{
	enum { N=4 };
	using element_t = double;
	using simd_t = __m256d;
	using my_type = Vector<element_t,N,false>;

	inline Vector() : m_vec( _mm256_setzero_pd() ) { }
 	inline Vector( Vector<element_t,N>&& v ) noexcept : m_vec(v()) { }
	inline Vector( const Vector<element_t,N>& v ) : m_vec(v()) { }
	inline Vector( simd_t v ) : m_vec(v) { }
	inline Vector( double e ) : m_vec( _mm256_set1_pd(e) ) { }
	inline Vector( const double* const e, std::size_t stride=1 )
	{
		switch(stride)
		{
			case 1: m_vec = _mm256_loadu_pd(e); break;
			default: for( std::size_t i=0; i < N; ++i ) { m_elem[i] = *(e+i*stride); }
		}
	}

	inline my_type& operator=( Vector<element_t,N>&& v ) { m_vec = v(); return *this; }
	inline my_type& operator=( const Vector<element_t,N>& v ) { m_vec = v(); return *this; }

	inline my_type& operator=( simd_t v ) { m_vec = v; return *this; }

	inline element_t* data() { return m_elem; }
	inline const element_t* data() const { return m_elem; }

	inline simd_t operator()() { return m_vec; }
	inline simd_t operator()() const { return m_vec; }

	inline element_t& operator[]( uint i ) { return m_elem[i]; }
	inline const element_t& operator[]( uint i ) const { return m_elem[i]; }

	template< bool ViewFlag >
	inline my_type& operator+=( const Vector<element_t,N,ViewFlag>& v ) { m_vec += v(); return *this; }
	template< bool ViewFlag >
	inline my_type& operator-=( const Vector<element_t,N,ViewFlag>& v ) { m_vec -= v(); return *this; }
	template< bool ViewFlag >
	inline my_type& operator*=( const Vector<element_t,N,ViewFlag>& v ) { m_vec *= v(); return *this; }
	template< bool ViewFlag >
	inline my_type& operator/=( const Vector<element_t,N,ViewFlag>& v ) { m_vec /= v(); return *this; }

	inline my_type& operator+=( const simd_t v ) { m_vec += v; return *this; }
	inline my_type& operator-=( const simd_t v ) { m_vec -= v; return *this; }
	inline my_type& operator*=( const simd_t v ) { m_vec *= v; return *this; }
	inline my_type& operator/=( const simd_t v ) { m_vec /= v; return *this; }

	inline my_type& operator+=( element_t e ) { m_vec += _mm256_set1_pd(e); return *this; }
	inline my_type& operator-=( element_t e ) { m_vec -= _mm256_set1_pd(e); return *this; }
	inline my_type& operator*=( element_t e ) { m_vec *= _mm256_set1_pd(e); return *this; }
	inline my_type& operator/=( element_t e ) { m_vec /= _mm256_set1_pd(e); return *this; }

	union
	{
		element_t m_elem[N];
		simd_t m_vec;
	};
};

template<>
struct Vector<double,4,true>
{
	enum { N=4 };
	using element_t = double;
	using simd_t = __m256d;
	using my_type = Vector<element_t,N,true>;

	inline Vector( element_t* const p ) : m_p(p) { }
	// generic
	inline my_type& operator=( Vector<element_t,N>&& v ) { this->store( v() ); return *this; }
	inline my_type& operator=( const Vector<element_t,N>& v ) { this->store( v() ); return *this; }
	// my_type
	inline my_type& operator=( my_type&& v ) { this->store( v() ); return *this; }
	inline my_type& operator=( const my_type& v ) { this->store( v() ); return *this; }

	inline my_type& operator=( simd_t v ) { this->store( v ); return *this; }

	inline element_t* data() { return m_p; }
	inline const element_t* data() const { return m_p; }

	inline simd_t operator()() { return _mm256_load_pd(m_p); }
	inline simd_t operator()() const { return _mm256_load_pd(m_p); }

	inline element_t& operator[]( uint i ) { return *(m_p+i); }
	inline const element_t& operator[]( uint i ) const { return *(m_p+i); }

	inline void store( simd_t vec ) { _mm256_store_pd( m_p, vec ); }

	template< bool ViewFlag >
	inline my_type& operator+=( const Vector<element_t,N,ViewFlag>& v ) { this->store( _mm256_add_pd( (*this)(), v() ) ); return *this; }
	template< bool ViewFlag >
	inline my_type& operator-=( const Vector<element_t,N,ViewFlag>& v ) { this->store( _mm256_sub_pd( (*this)(), v() ) ); return *this; }
	template< bool ViewFlag >
	inline my_type& operator*=( const Vector<element_t,N,ViewFlag>& v ) { this->store( _mm256_mul_pd( (*this)(), v() ) ); return *this; }
	template< bool ViewFlag >
	inline my_type& operator/=( const Vector<element_t,N,ViewFlag>& v ) { this->store( _mm256_div_pd( (*this)(), v() ) ); return *this; }

	inline my_type& operator+=( const simd_t v ) { this->store( _mm256_add_pd( (*this)(), v ) ); return *this; }
	inline my_type& operator-=( const simd_t v ) { this->store( _mm256_sub_pd( (*this)(), v ) ); return *this; }
	inline my_type& operator*=( const simd_t v ) { this->store( _mm256_mul_pd( (*this)(), v ) ); return *this; }
	inline my_type& operator/=( const simd_t v ) { this->store( _mm256_div_pd( (*this)(), v ) ); return *this; }

	inline my_type& operator+=( element_t e ) { this->store( _mm256_add_pd( (*this)(), _mm256_set1_pd(e) ) ); return *this; }
	inline my_type& operator-=( element_t e ) { this->store( _mm256_sub_pd( (*this)(), _mm256_set1_pd(e) ) ); return *this; }
	inline my_type& operator*=( element_t e ) { this->store( _mm256_mul_pd( (*this)(), _mm256_set1_pd(e) ) ); return *this; }
	inline my_type& operator/=( element_t e ) { this->store( _mm256_div_pd( (*this)(), _mm256_set1_pd(e) ) ); return *this; }

	element_t* const m_p;
};

#endif // __AVX__

#endif // #ifndef NO_INTRINSICS

} // namespace apegrunt

#endif // APEGRUNT_VECTOR_INTERFACE_HPP
