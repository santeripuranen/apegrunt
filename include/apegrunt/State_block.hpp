/** @file State_block.hpp

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

#ifndef APEGRUNT_STATE_BLOCK_HPP
#define APEGRUNT_STATE_BLOCK_HPP

#include <cctype>
#include <functional> // for std::hash
#include <cstring> // for std::memcpy

#include <bitset>
#include <limits>

#include <boost/serialization/serialization.hpp>
#include <boost/serialization/nvp.hpp>

#include "StateVector_state_types.hpp"
#include "aligned_allocator.hpp"

#include "misc/SIMD_intrinsics.h"

namespace apegrunt {

template< typename MaskT, std::size_t N, typename StencilT >
constexpr MaskT create_clear( StencilT stencil ) { MaskT mask=0; uint n=0; while( n < N ) { mask = mask | ( MaskT(stencil) << n*8*sizeof(StencilT) ); ++n; } return mask; }

// forward declaration
template< typename StateT, std::size_t Size > struct State_block;

// const_State_view (a const proxy) wraps a raw state type and enables state references.
// /*
template< typename StateT, std::size_t Size >
struct const_State_view
{
	enum { N=Size };
	using state_t = StateT;
	using const_state_view_t = const_State_view<state_t,N>;

	using state_block_t = apegrunt::State_block<state_t,Size>;

	using my_type = const_State_view<state_t,N>;

	const_State_view() = delete;
	const_State_view( const state_block_t& parent, std::size_t pos ) : m_parent(parent), m_pos(pos) { assert(m_pos<Size); }
	~const_State_view() = default;

	const state_block_t& m_parent;
	const std::size_t m_pos;

	//template< typename Q = state_t, std::enable_if_t<!std::is_same<Q, char>::value,bool> = true >
	inline operator state_t() const { return m_parent.m_states[m_pos]; }

	//template< typename IntegerT, typename = typename std::enable_if< std::is_integral<IntegerT>::value && !std::is_same<IntegerT, char>::value >::type >
	template< typename IntegerT, std::enable_if_t< std::is_integral<IntegerT>::value && !std::is_same<IntegerT, char>::value, bool > = true >
	inline operator IntegerT() const { return IntegerT( this->operator state_t() ); }

	//template< typename Q = state_t, typename = std::enable_if_t<!std::is_same<Q, char>::value> >
	template< typename Q = state_t, std::enable_if_t<!std::is_same<Q, char>::value,bool> = true >
	inline operator char() const { return to_char( state_t(*this) ); }

	inline bool operator==( const my_type& rhs ) const { return state_t(*this) == state_t(rhs); }
	// we choose to define state_t order based on their internal numerical value; an arbitrary choice, of course, since we deal with categorical variables
	inline bool operator<( const my_type& rhs ) const { return state_t(*this) < state_t(rhs); }
};
// */
/*
template< typename StateT, std::size_t Size >
struct const_State_view
{
	enum { N=Size };
	using state_t = StateT;
	using const_state_view_t = const_State_view<state_t,N>;

	using state_block_t = apegrunt::State_block<state_t,Size>;

	using my_type = const_State_view<state_t,N>;

	const_State_view() = delete;
	const_State_view( const state_t& state ) : m_state(state) { }
	~const_State_view() = default;

	const state_t& m_state;

	//template< typename Q = state_t, std::enable_if_t<!std::is_same<Q, char>::value,bool> = true >
	inline operator state_t() const { return m_state; }

	//template< typename IntegerT, typename = typename std::enable_if< std::is_integral<IntegerT>::value && !std::is_same<IntegerT, char>::value >::type >
	template< typename IntegerT, std::enable_if_t< std::is_integral<IntegerT>::value && !std::is_same<IntegerT, char>::value, bool > = true >
	//inline operator IntegerT() const { return IntegerT( this->operator state_t() ); }
	inline operator IntegerT() const { return IntegerT( m_state ); }

	//template< typename Q = state_t, typename = std::enable_if_t<!std::is_same<Q, char>::value> >
	template< typename Q = state_t, std::enable_if_t<!std::is_same<Q, char>::value,bool> = true >
	inline operator char() const { return to_char( state_t(*this) ); }

	inline bool operator==( const my_type& rhs ) const { return state_t(*this) == state_t(rhs); }
	// we choose to define state_t order based on their internal numerical value; an arbitrary choice, of course, since we deal with categorical variables
	inline bool operator<( const my_type& rhs ) const { return state_t(*this) < state_t(rhs); }
};
*/

template< typename StateT, std::size_t Size >
std::ostream& operator<< ( std::ostream& os, const const_State_view<StateT,Size>& state )
{
	os << to_char(StateT(state));
	return os;
}

template< typename StateT, std::size_t Size >
constexpr StateT to_state( const const_State_view<StateT,Size>& val ) { return val; }


// State_view (a mutator proxy) enables us to take a reference of, and modify, otherwise unreferenceable state
// types (specifically those that are stored using less that a byte worth of bits).
// prototype class template; please specialize
template< typename StateT, std::size_t Size >
struct State_view
{
	enum { N=Size };
	using state_t = StateT;
	using stateholder_t = apegrunt::State_holder<state_t>;
	using const_state_view_t = const_State_view<state_t,N>;

	using state_block_t = apegrunt::State_block<state_t,Size>;
	//using integer_t = typename state_block_t::internal_t;

	using my_type = State_view<state_t,N>;

	State_view() = delete;
	State_view( state_block_t& parent, std::size_t pos ) : m_parent(parent), m_pos(pos) { assert(m_pos<Size); }
	~State_view() = default;

	state_block_t& m_parent;
	const std::size_t m_pos;

	//template< typename Q = state_t, std::enable_if_t<!std::is_same<Q, char>::value,bool> = true >
	inline operator state_t() const { return m_parent.m_states[m_pos]; }

	//template< typename IntegerT, typename = typename std::enable_if< std::is_integral<IntegerT>::value && !std::is_same<IntegerT, char>::value >::type >
	template< typename IntegerT, std::enable_if_t< std::is_integral<IntegerT>::value && !std::is_same<IntegerT, char>::value, bool > = true >
	inline operator IntegerT() const { return IntegerT( this->operator state_t() ); }

	template< typename Q = state_t, std::enable_if_t<!std::is_same<Q, char>::value,bool> = true >
	inline operator char() const { return to_char( state_t(*this) ); }

	inline my_type& operator=( state_t state )
	{
		m_parent.m_states[m_pos] = state;
		return *this;
	}

	inline my_type& operator=( const_state_view_t state )
	{
		return *this = state_t(state); // get the concrete state from const_state_view_t and forward
	}

	inline my_type& operator=( stateholder_t state )
	{
		return *this = state_t(state); // peel off stateholder and forward
	}

	inline bool operator==( const my_type& rhs ) const { return state_t(*this) == state_t(rhs); }
	// we choose to define state_t order based on their internal numerical value; an arbitrary choice, of course, since we deal with categorical variables
	inline bool operator<( const my_type& rhs ) const { return state_t(*this) < state_t(rhs); }
};

template< typename StateT, std::size_t Size >
std::ostream& operator<< ( std::ostream& os, const State_view<StateT,Size>& state )
{
	os << to_char(StateT(state));
	return os;
}

template< typename StateT, std::size_t Size >
constexpr StateT to_state( const State_view<StateT,Size>& val ) { return val; }



template< typename StateT, std::size_t Size, std::size_t ParentSize >
struct State_block_view
{
	enum { N=Size };
	using state_t = StateT;
	using state_block_t = apegrunt::State_block<state_t,ParentSize>;

	using my_type = State_block_view<state_t,N,ParentSize>;

	State_block_view() = delete;
	State_block_view( const state_block_t& parent, std::size_t pos ) : m_parent(parent), m_pos(pos) { assert(m_pos+Size<=ParentSize); }
	~State_block_view() = default;

	const state_block_t& m_parent;
	const std::size_t m_pos;

	template< std::size_t BlockSize >
	inline operator State_block<state_t,BlockSize>() const { assert(BlockSize<=ParentSize-m_pos); return m_parent.m_states+m_pos; }

	template< std::size_t BlockSize >
	inline int cmp( const State_block<state_t,BlockSize>& rhs ) { return std::memcmp(m_parent.m_states+m_pos, rhs.m_states, sizeof(state_t)*BlockSize); }

	template< std::size_t BlockSize >
	inline bool operator==( const State_block<state_t,BlockSize>& rhs ) { return 0 == std::memcmp(m_parent.m_states+m_pos, rhs.m_states, sizeof(state_t)*BlockSize); }

	template< std::size_t BlockSize >
	inline bool operator<( const State_block<state_t,BlockSize>& rhs ) { return 0 > std::memcmp(m_parent.m_states+m_pos, rhs.m_states, sizeof(state_t)*BlockSize); }

	template< std::size_t BlockSize >
	inline bool operator>( const State_block<state_t,BlockSize>& rhs ) { return 0 < std::memcmp(m_parent.m_states+m_pos, rhs.m_states, sizeof(state_t)*BlockSize); }
};

template< typename StateT, std::size_t Size >
//struct alignas(sizeof(void*)) State_block
struct alignas(std::max(std::size_t(8),Size*sizeof(StateT))) State_block
{
	using state_t = StateT;
	enum { N=Size };
	enum : uint8_t { CLEAR=uint8_t(gap_state<state_t>::value) };
	//enum : uint8_t { CLEAR=uint8_t(0) };

	using my_type = State_block< state_t, N >;
	//using stateholder_t = apegrunt::State_holder<state_t>;
	using const_state_view_t = const_State_view<state_t,N>;
	using state_view_t = State_view<state_t,N>;

	template<std::size_t BlockSize> using sub_block_t = State_block<state_t,BlockSize>;
	template<std::size_t BlockSize> using sub_block_view_t = State_block_view<state_t,BlockSize,N>;

	state_t m_states[N];

	State_block() { this->clear(); }
	~State_block() = default;

	State_block( my_type&& other ) noexcept { this->set(other.m_states); }
	State_block( const my_type& other ) { this->set(other.m_states); }
	template< typename Q = state_t, typename = std::enable_if_t<!std::is_same<Q, char>::value> >
	State_block( const char* states, std::size_t n=N ) { this->set(states,n); }
	State_block( const state_t* states, std::size_t n=N ) { this->set(states,n); }

	template< typename Q = state_t, typename = std::enable_if_t<!std::is_same<Q, char>::value> >
	State_block( const State_block<char,N>& other ) { this->set(other.m_states); }

	my_type& operator=( my_type&& other ) noexcept { this->set(other.m_states); return *this; }
	my_type& operator=( const my_type& other ) { this->set(other.m_states); return *this; }

	// No range check
	template< typename Q = state_t, std::enable_if_t<!std::is_same<Q, char>::value,bool> = true >
	inline void set( const char* states, std::size_t n=N ) { assert(n<=N); for(uint i=0; i<n; ++i) { m_states[i] = apegrunt::char_to_state<state_t>(*(states+i)); } }
	inline void set( const state_t* states, std::size_t n=N ) { assert(n<=N); std::memcpy(m_states,states,sizeof(state_t)*n); }

	// No range check
	inline state_view_t operator[]( std::size_t pos ) { assert(pos<N); return state_view_t(*this,pos); }
	//inline const_state_view_t operator[]( std::size_t pos ) const { return m_states[pos]; }
	inline const_state_view_t operator[]( std::size_t pos ) const { assert(pos<N); return const_state_view_t(*this,pos); }

	// No range check
	template< std::size_t BlockSize >
	inline sub_block_t<BlockSize> block( std::size_t pos ) const { assert((pos+1)*BlockSize<=N); return m_states+pos*BlockSize; }

	// No range check
	template< std::size_t BlockSize >
	inline sub_block_view_t<BlockSize> blockview( std::size_t pos ) const { assert((pos+1)*BlockSize<=N); return sub_block_view_t<BlockSize>(*this,pos*BlockSize); }

	inline void clear() { std::memset(m_states,uint8_t(CLEAR),sizeof(state_t)*N); }

	inline int cmp( const my_type& rhs ) const { return std::memcmp(m_states, rhs.m_states, sizeof(state_t)*N); }

	inline bool operator==( const my_type& rhs ) const { return 0 == std::memcmp(m_states, rhs.m_states, sizeof(state_t)*N); }

	// this is rather arbitrary since state_t are categorical variables, but we define it nonetheless in order to support sorting.
	inline bool operator<( const my_type& rhs ) const { return 0 > std::memcmp(m_states, rhs.m_states, sizeof(state_t)*N); }
	inline bool operator>( const my_type& rhs ) const { return 0 < std::memcmp(m_states, rhs.m_states, sizeof(state_t)*N); }

	//inline my_type& operator|=( const my_type& rhs ) { for( auto i(0); i<N; ++i ) { m_ints[i] = ((m_ints[i] | rhs.m_ints[i]) & 0x3) & (((m_ints[i] & CLEAR)^(rhs.m_ints[i] & CLEAR)) ^ CLEAR); } return *this; }
	inline my_type& operator|=( const my_type& rhs )
	{
		for( auto i(0); i<N; ++i )
		{
			switch(m_states[i])
			{
			case gap_state<state_t>::value:
				switch(rhs.m_states[i])
				{
				case gap_state<state_t>::value: break;
				default: m_states[i] = rhs.m_states[i]; break;
				} break;
			default: break;
			}
		}
		return *this;
	}
};

template< typename BlockT >
using State_block_allocator_t = memory::AlignedAllocator<BlockT,alignof(BlockT)>;

template< typename StateT, std::size_t Size >
inline State_block<StateT,Size> operator|( const State_block<StateT,Size>& lhs, const State_block<StateT,Size>& rhs )
{
	State_block<StateT,Size> combined(lhs);
	combined |= rhs;
	return combined;
}

template< typename StateT, std::size_t Size >
constexpr std::size_t size( const State_block< StateT, Size >& block ) { return Size; }

template< typename StateT, std::size_t Size >
inline bool operator!=( const State_block<StateT,Size>& lhs, const State_block<StateT,Size>& rhs ) { return !(lhs == rhs); }

template< typename StateT, std::size_t Size >
inline bool operator> ( const State_block<StateT,Size>& lhs, const State_block<StateT,Size>& rhs ) { return (rhs < lhs); }

template< typename StateT, std::size_t Size >
inline bool operator<=( const State_block<StateT,Size>& lhs, const State_block<StateT,Size>& rhs ) { return !(lhs > rhs); }

template< typename StateT, std::size_t Size >
inline bool operator>=( const State_block<StateT,Size>& lhs, const State_block<StateT,Size>& rhs ) { return !(lhs < rhs); }

template< typename StateT, std::size_t Size >
std::ostream& operator<< ( std::ostream& os, const State_block<StateT,Size>& block )
{
	for( std::size_t i=0; i < Size; ++i )
	{
		os << char(block[i]);
	}
	return os;
}

template< typename StateT, std::size_t Size >
std::size_t count_identical( const State_block<StateT,Size>& lhs, const State_block<StateT,Size>& rhs )
{
	std::size_t n=0;
	for( std::size_t i=0; i<Size; ++i ) { n += (lhs[i] == rhs[i]); }
	return n;
}

#ifndef NO_INTRINSICS

#ifdef __POPCNT__

// 16 states per block
#ifdef __SSE2__
template< typename StateT >
inline std::size_t count_identical( State_block<StateT,16> lhs, State_block<StateT,16> rhs )
{
	auto load = [](const StateT *const p){ return _mm_load_si128( (__m128i*)p ); };
	return _mm_popcnt_u32( _mm_movemask_epi8( _mm_cmpeq_epi8( load(lhs.m_states), load(rhs.m_states) ) ) );
}
#endif // __SSE2__

// 32 states per block
#ifdef __AVX2__
template< typename StateT >
inline std::size_t count_identical( const State_block<StateT,32>& lhs, const State_block<StateT,32>& rhs )
{
	//std::cout << "count_identical<avx2>" << std::endl;
	auto load = [](const StateT *const p){ return _mm256_lddqu_si256( (__m256i*)p ); }; // may perform better than _mm256_loadu_si256, when loads cross cache-line boundary
	//auto load = [](const StateT *const p){ return _mm256_loadu_si256( (__m256i*)p ); }; // hmm.. _mm256_loadu_si256 // our data is 32B-aligned; loads will never cross cache-line boundary
	//auto load = [](const StateT *const p){ return _mm256_load_si256( (__m256i*)p ); }; // our data is 32B-aligned; loads will never cross cache-line boundary
	return _mm_popcnt_u32( _mm256_movemask_epi8( _mm256_cmpeq_epi8( load(lhs.m_states), load(rhs.m_states) ) ) );
}
#else if __SSE2__
template< typename StateT >
inline std::size_t count_identical( State_block<StateT,32> lhs, State_block<StateT,32> rhs )
{
	//std::cout << "count_identical<sse2>" << std::endl;
//	auto load = [](const StateT *const p){ return _mm_load_si128( (__m128i*)p ); };
//	return	_mm_popcnt_u32( _mm_movemask_epi8( _mm_cmpeq_epi8( load(lhs.m_states), load(rhs.m_states) ) ) ) +
//			_mm_popcnt_u32( _mm_movemask_epi8( _mm_cmpeq_epi8( load(lhs.m_states+16), load(rhs.m_states+16) ) ) );
	auto load = [](const StateT *const p, std::size_t offset){ return _mm_load_si128( (__m128i*)(p+offset) ); };
	return	_mm_popcnt_u32(
				uint32_t(_mm_movemask_epi8( _mm_cmpeq_epi8( load(lhs.m_states,0), load(rhs.m_states,0) ) ) ) |
				( uint32_t(_mm_movemask_epi8( _mm_cmpeq_epi8( load(lhs.m_states,16), load(rhs.m_states,16) ) ) ) << 16 )
			);
}
#endif // __SSE2__

// 64 states per block
#ifdef __AVX2__
template< typename StateT >
inline std::size_t count_identical( const State_block<StateT,64>& lhs, const State_block<StateT,64>& rhs )
{
	auto load = [](const StateT *const p, int offset){ return _mm256_lddqu_si256( (__m256i*)(p+offset) ); }; // may perform better than _mm256_loadu_si256, when loads cross cache-line boundary
	//auto load = [](const StateT *const p, int offset){ return _mm256_loadu_si256( (__m256i*)(p+offset) ); }; // hmm.. _mm256_loadu_si256 will produce the wrong result on zen2 // our data is 64B-aligned; loads will never cross cache-line boundary
	//auto load = [](const StateT *const p, int offset){ return _mm256_load_si256( (__m256i*)(p+offset) ); }; // our data is 64B-aligned; loads will never cross cache-line boundary
	return	_mm_popcnt_u32( _mm256_movemask_epi8( _mm256_cmpeq_epi8( load(lhs.m_states,0), load(rhs.m_states,0) ) ) ) +
			_mm_popcnt_u32( _mm256_movemask_epi8( _mm256_cmpeq_epi8( load(lhs.m_states,32), load(rhs.m_states,32) ) ) );
// this one is broken (doesn't crash, but gives the wrong result):
//	return	_mm_popcnt_u64(
//				uint64_t(_mm256_movemask_epi8( _mm256_cmpeq_epi8( load(lhs.m_states,0), load(rhs.m_states,0) ) ) ) |
//				( uint64_t(_mm256_movemask_epi8( _mm256_cmpeq_epi8( load(lhs.m_states,32), load(rhs.m_states,32) ) ) ) << 32 )
//			);
}
#else if __SSE2__
template< typename StateT >
inline std::size_t count_identical( State_block<StateT,64> lhs, State_block<StateT,64> rhs )
{
	auto load = [](const StateT *const p, int offset){ return _mm_load_si128( (__m128i*)(p+offset) ); }; // hmm.. _mm_loadu_si128 will produce the wrong result on zen2
// this version is correct, but perhaps a tiny bit slower than the alternative below
//	return	_mm_popcnt_u32( _mm_movemask_epi8( _mm_cmpeq_epi8( load(lhs.m_states,0), load(rhs.m_states,0) ) ) ) +
//			_mm_popcnt_u32( _mm_movemask_epi8( _mm_cmpeq_epi8( load(lhs.m_states,16), load(rhs.m_states,16) ) ) ) +
//			_mm_popcnt_u32( _mm_movemask_epi8( _mm_cmpeq_epi8( load(lhs.m_states,32), load(rhs.m_states,32) ) ) ) +
//			_mm_popcnt_u32( _mm_movemask_epi8( _mm_cmpeq_epi8( load(lhs.m_states,48), load(rhs.m_states,48) ) ) )
//			;
	return	_mm_popcnt_u64(
				(uint64_t)_mm_movemask_epi8( _mm_cmpeq_epi8( load(lhs.m_states,0), load(rhs.m_states,0) ) ) |
				( (uint64_t)_mm_movemask_epi8( _mm_cmpeq_epi8( load(lhs.m_states,16), load(rhs.m_states,16) ) ) << 16 ) |
				( (uint64_t)_mm_movemask_epi8( _mm_cmpeq_epi8( load(lhs.m_states,32), load(rhs.m_states,32) ) ) << 32 ) |
				( (uint64_t)_mm_movemask_epi8( _mm_cmpeq_epi8( load(lhs.m_states,48), load(rhs.m_states,48) ) ) << 48 )
			);
}
#endif // __SSE2__

#endif // __POPCNT__

#endif // #ifndef NO_INTRINSICS

} // namespace apegrunt

#endif // APEGRUNT_STATE_BLOCK_HPP

