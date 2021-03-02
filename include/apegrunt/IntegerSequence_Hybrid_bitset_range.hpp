/** @file IntegerSequence_Hybrid_bitset_range.hpp
 
	Copyright (c) 2018-2021 Santeri Puranen.
 
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

#ifndef APEGRUNT_INTEGERSEQUENCE_HYBRID_BITSET_RANGE_HPP
#define APEGRUNT_INTEGERSEQUENCE_HYBRID_BITSET_RANGE_HPP

#include <limits>
#include <ostream>
#include <type_traits> // for std::integral_constant
#include <numeric> // for std::accumulate

#include "IntegerSequence_Hybrid_bitset_range_forward.h"
#include "misc/Vector_operations.hpp" // basically only for apegrunt::popcnt
#include "misc/Math.hpp" // for apegrunt::ipow

namespace apegrunt {

template< typename ValueT >
struct Hybrid_bitset_range_element_internal;

template<> struct Hybrid_bitset_range_element_internal<uint16_t> { using type=uint8_t; static const uint8_t BLOCK_MASK=0x1; };
template<> struct Hybrid_bitset_range_element_internal<uint32_t> { using type=uint16_t; static const uint16_t BLOCK_MASK=0x1FF; };
template<> struct Hybrid_bitset_range_element_internal<uint64_t> { using type=uint32_t; static const uint32_t BLOCK_MASK=0x1FFFFFF; };

template< typename IndexT, bool Aligned >
struct alignas(sizeof(IndexT)) Hybrid_bitset_range_element
{
#pragma message("NOTE: internal range is 10 bits for value_type=uint16_t, 19 bits for value_type=uint32_t and 36 bits for value_type=uint64_t")
	using value_type = IndexT;
	using index_t = typename Hybrid_bitset_range_element_internal<value_type>::type;
	using my_type = Hybrid_bitset_range_element<value_type,Aligned>;
	using aligned = std::integral_constant<bool, Aligned>;
	using block_mask_type = uint64_t;

	static constexpr index_t create_range_flag_mask() { return index_t(1) << std::numeric_limits<index_t>::digits-1; } // set the most significant bit
	static constexpr index_t popcnt( index_t mask ) { index_t n(0); while( mask ) { ++n; mask >>= 1; } return n; }

	enum { ELEMENT_SIZE = std::numeric_limits<value_type>::digits }; // 16, 32, 64
	enum { BITSTRING_SIZE = std::numeric_limits<index_t>::digits }; // 8, 16, 32
	enum { RANGE_FLAG_MASK = my_type::create_range_flag_mask() };
	enum { RANGE_UNFLAG_MASK = (index_t)~my_type::create_range_flag_mask() };
	enum { MODULO_MASK=std::numeric_limits<index_t>::digits-1 }; // 0x7 for 8-bit, 0xF for 16-bit, 0x1F for 32-bit and 0x3F for 64-bit
	enum { MODULO_NSTEPS=my_type::popcnt(MODULO_MASK) }; // 3, 4, 5
	enum { RANGE = std::numeric_limits<index_t>::digits-1 + my_type::popcnt(MODULO_MASK) }; // numeric range in number of bits: 10 bits, 19 bits, 36 bits
	enum { BLOCK_MASK=Hybrid_bitset_range_element_internal<value_type>::BLOCK_MASK };
	enum { BLOCK_NSTEPS=my_type::popcnt(BLOCK_MASK) };
	enum { INTL_BLOCK_NSTEPS=my_type::popcnt(BLOCK_MASK)-MODULO_NSTEPS };

	Hybrid_bitset_range_element( ) : m_payload(0) { } // zero-initialize
	explicit Hybrid_bitset_range_element( value_type value ) : m_pos( value >> MODULO_NSTEPS ), m_bitfield(1 << (value & MODULO_MASK)) { }
	Hybrid_bitset_range_element( value_type value, index_t bitfield ) : m_pos(value >> MODULO_NSTEPS), m_bitfield(bitfield) { } // trust that 'value' is properly aligned
	Hybrid_bitset_range_element( const my_type& other, index_t bitfield ) noexcept : m_pos(other.m_pos & RANGE_UNFLAG_MASK), m_bitfield(bitfield) { } // other could be a range, so strip any range flag
	Hybrid_bitset_range_element( const my_type& other ) : m_payload(other.m_payload) { }
	Hybrid_bitset_range_element( my_type&& other ) : m_payload(std::move(other.m_payload)) { }

	inline my_type& operator=( const my_type& other ) { m_payload = other.m_payload; return *this; }
	inline my_type& operator=( my_type&& other ) { m_payload = std::move(other.m_payload); return *this; }
	inline my_type& operator=( value_type value ) { m_pos = (value >> MODULO_NSTEPS); m_bitfield = 1 << (value & MODULO_MASK); return *this; }

	inline value_type operator()() const noexcept { return value_type(m_pos & RANGE_UNFLAG_MASK) << MODULO_NSTEPS; } // range flag *really* needs to be removed

	// the bitsets are equal only if they exactly match each other, range property and all
	inline bool operator==( const my_type& other ) const { return m_payload == other.m_payload; }
	inline bool operator==( const value_type& value ) const { return this->range_begin() <= value && this->range_end() > value; }

	inline bool operator<( const my_type& other ) const { return ( (m_pos & RANGE_UNFLAG_MASK) < (other.m_pos & RANGE_UNFLAG_MASK) ); }
	inline bool operator<( const value_type& value ) const { return this->range_end() <= value; }

	inline operator bool() const { return bool(m_payload); } // the element is empty if the entire payload contains no set bits (note that a stored value of "zero" would have the least significant bit of m_bitfield set)

	inline my_type& set( value_type value, index_t bitfield ) noexcept { m_pos = value >> MODULO_NSTEPS; m_bitfield = bitfield; return *this; } // trust that 'value' is properly aligned
	inline my_type& set( const my_type& other, index_t bitfield ) noexcept { m_pos = other.m_pos & RANGE_UNFLAG_MASK; m_bitfield = bitfield; return *this; } // other could be a range, so strip any range flag

	// try to set bit at 'pos' and return whether successful or not: 'false' if out of range or if range flag is set.
	// can be safely called on range elements (never modifies, always returns 'false')
	inline bool set( value_type pos ) noexcept
	{
		// will return "false" if "this" is a range element -- it doesn't make sense to set individual bits in a range
		return ((pos >> MODULO_NSTEPS) == m_pos) && (m_bitfield |= (1 << (pos & MODULO_MASK))); // '(pos >> MODULO_NSTEPS) == m_pos' will be false for range elements
	}

	inline bool is_set( value_type pos ) const noexcept
	{
		// will return "false" if "this" is a range element
		return ((pos >> MODULO_NSTEPS) == m_pos) && ( m_bitfield & (1 << (pos & MODULO_MASK)) );
	}

	inline bool all_set() const noexcept { return !bool( (index_t)~m_bitfield ); } // cast needed when index_t == uint16_t (at least for gcc)

	inline bool holds( value_type pos ) const noexcept { return this->is_range() ? (pos>=this->range_begin() && pos<this->range_end()) : this->is_set(pos); }

	inline bool is_range() const noexcept { return bool( m_pos & RANGE_FLAG_MASK ); }
	inline my_type& set_range() noexcept { m_pos = index_t(m_pos) | index_t(RANGE_FLAG_MASK); return *this; }
	inline my_type& set_range( value_type pos, value_type end ) noexcept { if( pos < end ) { m_pos = (index_t(pos >> MODULO_NSTEPS) | index_t(RANGE_FLAG_MASK)); m_bitfield = index_t(end >> MODULO_NSTEPS)/*-(m_pos & RANGE_UNFLAG_MASK)*/; } return *this; }
	inline my_type& set_range_end( value_type end ) noexcept { m_bitfield = index_t(end >> MODULO_NSTEPS)/*-(m_pos & RANGE_UNFLAG_MASK)*/; return this->set_range(); } // m_bitfield is range size counter (0=no range, 1=BITSTRING_SIZE, 2=BITSTRING_SIZE*2, etc)

	inline value_type range_begin() const noexcept { return (*this)(); }
	inline value_type range_end() const noexcept { return value_type(this->is_range() ? m_bitfield : m_pos+1) << MODULO_NSTEPS; } // range_begin != range_end and past the end even for non-range elements

	inline index_t get_bitfield() const noexcept { return m_bitfield; }
	inline void set_bitfield( index_t bitfield ) { m_bitfield = bitfield; }
	inline void merge_bitfield( index_t bitfield ) { m_bitfield |= bitfield; }

	inline std::size_t size() const
	{
		return this->is_range() ? std::size_t( m_bitfield - (m_pos & RANGE_UNFLAG_MASK) ) << MODULO_NSTEPS : apegrunt::popcnt(m_bitfield);
	}

	friend inline void swap( my_type& lhs, my_type& rhs ) { using std::swap; swap( lhs.m_payload, rhs.m_payload ); }

	// Although every major compiler will accept this,
	// type punning as used here is not standards-conforming C++,
	// so potentially a trip to undefined behavior land..
	// ..and the memory debugger will complain
	//
	union {
		struct {
			index_t m_bitfield;
			index_t m_pos;
		};
		value_type m_payload;
	};

	template< typename RealT >
	static inline RealT gather( const RealT *weights, index_t selection_mask );

	template< typename RealT >
	static inline RealT gather( const RealT *weights, const my_type& selection_range );
};

template< typename IndexT >
inline Hybrid_bitset_range_element<IndexT,true> operator&( const Hybrid_bitset_range_element<IndexT,true>& lhs, const Hybrid_bitset_range_element<IndexT,true>& rhs )
{
	using element_t = Hybrid_bitset_range_element<IndexT,true>;
	element_t value; value.m_payload = lhs.m_payload & rhs.m_payload;
	return value;
}

template< typename IndexT, bool Aligned > template< typename RealT >
RealT Hybrid_bitset_range_element<IndexT,Aligned>::gather( const RealT *val_ptr, typename Hybrid_bitset_range_element<IndexT,Aligned>::index_t selection_mask )
{
	RealT sum(0);

	do { if(selection_mask & 1) { sum += *val_ptr; } } while( selection_mask >>= 1 && ++val_ptr );
	return sum;
}

template< typename IndexT, bool Aligned > template< typename RealT >
RealT Hybrid_bitset_range_element<IndexT,Aligned>::gather( const RealT *val_ptr, const Hybrid_bitset_range_element<IndexT,Aligned>& element )
{
	using element_t = Hybrid_bitset_range_element<IndexT,Aligned>;
	return element.is_range() ? std::accumulate( val_ptr+element.range_begin(), val_ptr+element.range_end(), RealT(0), []( auto sum, auto w ){ return sum+w; } ) : element_t::gather(val_ptr+element.range_begin(), element.get_bitfield() );
}

template< typename IndexT >
inline std::ostream& operator<< ( std::ostream& os, const Hybrid_bitset_range_element<IndexT,true>& element )
{
	os << element();

	if( element.is_range() ) { os << "-" << element.range_end(); }
	else { os << ":" << element.get_bitfield(); }

	return os;
}
// /*
template< typename IndexT, bool Aligned >
struct hash< IntegerSequence< Hybrid_bitset_range_element<IndexT,Aligned> > >
{
	std::size_t operator()( const IntegerSequence< Hybrid_bitset_range_element<IndexT,Aligned> >& s ) const { return s.hash(); }
};

template< typename IndexT, bool Aligned >
struct hash< std::reference_wrapper< IntegerSequence< Hybrid_bitset_range_element<IndexT,Aligned> > > >
{
	inline std::size_t operator()( const std::reference_wrapper< IntegerSequence< Hybrid_bitset_range_element<IndexT,Aligned> > >& wrapped ) const
	{
		return wrapped.get().hash();
	}
};

template< typename IndexT, bool Aligned >
struct equal_to< IntegerSequence< Hybrid_bitset_range_element<IndexT,Aligned> > >
{
	std::size_t operator()(
			const IntegerSequence< Hybrid_bitset_range_element<IndexT,Aligned> >& a,
			const IntegerSequence< Hybrid_bitset_range_element<IndexT,Aligned> >& b
		) const
	{ return a == b; }
};

template< typename IndexT, bool Aligned >
inline bool operator==(
		const std::reference_wrapper< IntegerSequence< Hybrid_bitset_range_element<IndexT,Aligned> > >& lhs,
		const std::reference_wrapper< IntegerSequence< Hybrid_bitset_range_element<IndexT,Aligned> > >& rhs
	)
{
	return lhs.get() == rhs.get();
}
// */

} // namespace apegrunt

#endif // APEGRUNT_INTEGERSEQUENCE_HYBRID_BITSET_RANGE_HPP
