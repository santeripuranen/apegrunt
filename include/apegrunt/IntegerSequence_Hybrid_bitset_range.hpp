/** @file IntegerSequence_Hybrid_bitset_range.hpp
 
	Copyright (c) 2018-2019 Santeri Puranen. All rights reserved.
 
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
/*
template< typename IndexT, bool Aligned >
struct Hybrid_bitset_range_element
{
	using index_t = IndexT;
	using value_type = index_t;
	using my_type = Hybrid_bitset_range_element<index_t,Aligned>;
	using aligned = std::integral_constant<bool, Aligned>;

	static constexpr index_t create_range_flag_mask() { return index_t(1) << std::numeric_limits<index_t>::digits-1; } // set the most significant bit

	enum { BITSTRING_SIZE = std::numeric_limits<index_t>::digits };
	//enum { BITSTRING_LAST_ELEMENT_INDEX = BITSTRING_SIZE-1 };
	enum { RANGE_FLAG_MASK = my_type::create_range_flag_mask() };
	enum { RANGE_UNFLAG_MASK = ~my_type::create_range_flag_mask() };
	enum { MODULO_MASK=std::numeric_limits<index_t>::digits-1 }; // 0x7 for 8-bit, 0xF for 16-bit, 0x1F for 32-bit and 0x3F for 64-bit

	Hybrid_bitset_range_element( ) : m_pos(0), m_bitfield(0) { } // zero-initialize
	Hybrid_bitset_range_element( index_t value ) : m_pos(value & ~MODULO_MASK), m_bitfield(RANGE_FLAG_MASK >> (value & MODULO_MASK)) { }
	Hybrid_bitset_range_element( const my_type& other ) : m_pos(other.m_pos), m_bitfield(other.m_bitfield) { }

	inline my_type& operator=( const my_type& other ) { m_pos = other.m_pos; m_bitfield = other.m_bitfield; return *this; }
	inline my_type& operator=( index_t value ) { m_pos = value & ~MODULO_MASK; m_bitfield = RANGE_FLAG_MASK >> (value & MODULO_MASK); return *this; }

	inline index_t operator()() const { return m_pos & RANGE_UNFLAG_MASK; }

	inline bool set( index_t pos ) noexcept
	{
		if( ( pos & ~MODULO_MASK ) == m_pos )
		{
			m_bitfield = m_bitfield | (RANGE_FLAG_MASK >> (pos & MODULO_MASK));
			return true;
		}
		return false;
	}
	inline bool is_set( index_t pos ) const noexcept
	{
		return (( pos & ~MODULO_MASK ) == m_pos) && ( m_bitfield & (RANGE_FLAG_MASK >> (pos & MODULO_MASK)) );
	}

	inline bool is_unset( index_t pos ) const { return !this->is_set(pos); }

	inline bool all_set() const { return !bool( ~m_bitfield ); }

	inline bool is_range() const { return bool( m_pos & RANGE_FLAG_MASK ); }
	inline void set_range() { m_pos = m_pos | RANGE_FLAG_MASK; }
	inline void set_range_end( index_t pos ) { this->set_range(); m_bitfield = pos; }

	inline index_t range_begin() const { return (*this)(); }
	inline index_t range_end() const { return this->is_range() ? m_bitfield : m_pos; }

	inline index_t get_bitfield() const { return m_bitfield; }

	index_t m_pos;
	index_t m_bitfield;
};
*/

namespace hybrid_bitset_range_detail {
template< typename ValueT >
struct Hybrid_bitset_range_element_internal;

template<> struct Hybrid_bitset_range_element_internal<uint16_t> { using type=uint8_t; static const uint8_t BLOCK_MASK=0x1; };
template<> struct Hybrid_bitset_range_element_internal<uint32_t> { using type=uint16_t; static const uint16_t BLOCK_MASK=0x1FF; };
template<> struct Hybrid_bitset_range_element_internal<uint64_t> { using type=uint32_t; static const uint8_t BLOCK_MASK=0x1FFFFFF; };
} // namespace hybrid_bitset_range_detail


template< typename IndexT, bool Aligned >
struct alignas(sizeof(IndexT)) Hybrid_bitset_range_element
{
#pragma message("NOTE: internal range is 10 bits for value_type=uint16_t, 19 bits for value_type=uint32_t and 36 bits for value_type=uint64_t")
	using value_type = IndexT;
	using index_t = typename hybrid_bitset_range_detail::Hybrid_bitset_range_element_internal<value_type>::type;
	using my_type = Hybrid_bitset_range_element<value_type,Aligned>;
	using aligned = std::integral_constant<bool, Aligned>;
	using block_mask_type = uint64_t;

	static constexpr index_t create_range_flag_mask() { return index_t(1) << std::numeric_limits<index_t>::digits-1; } // set the most significant bit
	static constexpr index_t count_set_bits( index_t mask ) { index_t n(0); while( mask ) { ++n; mask = mask >> 1; } return n; }

	enum { BITSTRING_SIZE = std::numeric_limits<index_t>::digits }; // 8, 16, 32
	enum { RANGE_FLAG_MASK = my_type::create_range_flag_mask() };
	enum { RANGE_UNFLAG_MASK = (index_t)~my_type::create_range_flag_mask() };
	enum { MODULO_MASK=std::numeric_limits<index_t>::digits-1 }; // 0x7 for 8-bit, 0xF for 16-bit, 0x1F for 32-bit and 0x3F for 64-bit
	enum { MODULO_NSTEPS=my_type::count_set_bits(MODULO_MASK) }; // 3, 4, 5
	enum { RANGE = std::numeric_limits<index_t>::digits-1 + my_type::count_set_bits(MODULO_MASK) }; // numeric range in number of bits; 10, 19, 36
	// 10-3-6=1, 19-4-6=9, 36-5-6=25
	enum { BLOCK_MASK=hybrid_bitset_range_detail::Hybrid_bitset_range_element_internal<value_type>::BLOCK_MASK };
	//enum { BLOCK_MASK=apegrunt::ipow(2,RANGE-MODULO_NSTEPS-my_type::count_set_bits(std::numeric_limits<block_mask_type>::digits-1))-1 }; // mask granularity is 1bit for uint16_t, 9bits for uint32_t and 25bits for uint64_t
	enum { BLOCK_NSTEPS=my_type::count_set_bits(BLOCK_MASK) };
	enum { INTL_BLOCK_NSTEPS=my_type::count_set_bits(BLOCK_MASK)-MODULO_NSTEPS };

	Hybrid_bitset_range_element( ) : m_payload(0) { } // zero-initialize
	Hybrid_bitset_range_element( value_type value ) : m_pos( value >> MODULO_NSTEPS ), m_bitfield(1 << (value & MODULO_MASK)) { }
	Hybrid_bitset_range_element( index_t value, index_t bitfield ) : m_pos(value >> MODULO_NSTEPS), m_bitfield(bitfield) { } // trust that 'value' is properly aligned

	Hybrid_bitset_range_element( const my_type& other ) : m_payload(other.m_payload) { }
	Hybrid_bitset_range_element( my_type&& other ) : m_payload(std::move(other.m_payload)) { }

	inline my_type& operator=( const my_type& other ) { m_payload = other.m_payload; return *this; }
	inline my_type& operator=( my_type&& other ) { m_payload = std::move(other.m_payload); return *this; }
	inline my_type& operator=( value_type value ) { m_pos = (value >> MODULO_NSTEPS); m_bitfield = 1 << (value & MODULO_MASK); return *this; }

	inline value_type operator()() const noexcept { return value_type(m_pos & RANGE_UNFLAG_MASK) << MODULO_NSTEPS; }

	// Hbres are equal only if they exactly match each other, range property and all
	inline bool operator==( const my_type& other ) const { return m_payload == other.m_payload; }

	inline bool operator<( const my_type& other ) const { return ( (m_pos & RANGE_UNFLAG_MASK) < (other.m_pos & RANGE_UNFLAG_MASK) ); }
	inline operator bool() const { return bool(m_payload); } // the element is empty if the entire payload contains no set bits; integer value "zero" would have the least significant bit of m_bitfield set

	// don't call on a range element; it's up to the caller to guard against this, for example by "if( !obj.is_range() ) { obj.set(val); }"
	inline bool set( value_type pos ) noexcept
	{
		// will return "false" if "this" is a range element
		return ((pos >> MODULO_NSTEPS) == m_pos) && (m_bitfield = (m_bitfield | (1 << (pos & MODULO_MASK))));
	}
	inline bool is_set( value_type pos ) const noexcept
	{
		// will return "false" if "this" is a range element
		return ((pos >> MODULO_NSTEPS) == m_pos) && ( m_bitfield & (1 << (pos & MODULO_MASK)) );
	}

	inline bool all_set() const noexcept { return !bool( (index_t)~m_bitfield ); }

	inline bool holds( value_type pos ) const noexcept { return this->is_range() ? (pos>=this->range_begin() && pos<this->range_end()) : this->is_set(pos); }

	inline bool is_range() const noexcept { return bool( m_pos & RANGE_FLAG_MASK ); }
	inline my_type& set_range() noexcept { m_pos = index_t(m_pos) | index_t(RANGE_FLAG_MASK); return *this; }
	inline my_type& set_range( value_type pos, value_type end ) noexcept { if( pos < end ) { m_pos = ((pos >> MODULO_NSTEPS) | RANGE_FLAG_MASK); m_bitfield = (end >> MODULO_NSTEPS); } return *this; }
	inline my_type& set_range_end( value_type pos ) noexcept { m_bitfield = (pos >> MODULO_NSTEPS); return this->set_range(); } // returns *this; end is past-the-end

	inline value_type range_begin() const noexcept { return (*this)(); }
	inline value_type range_end() const noexcept { return this->is_range() ? value_type(m_bitfield) << MODULO_NSTEPS : value_type(m_pos+1) << MODULO_NSTEPS; } // range_begin != range_end and past the end also for non-range elements

	inline index_t get_bitfield() const noexcept { return m_bitfield; }
	inline void set_bitfield( index_t bitfield ) { m_bitfield = bitfield; }
	inline void merge_bitfield( index_t bitfield ) { m_bitfield = m_bitfield | bitfield; }

	inline std::size_t size() const
	{
		return this->is_range() ? std::size_t( m_bitfield - (m_pos & RANGE_UNFLAG_MASK) ) << MODULO_NSTEPS : apegrunt::popcnt(m_bitfield);
	}

	friend inline void swap( my_type& lhs, my_type& rhs ) { using std::swap; swap( lhs.m_payload, rhs.m_payload ); }

	// Although every major compiler accept this,
	// type punning as used here is not standard-conforming C++,
	// so potentially a trip to undefined behavior land..
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

template< typename IndexT, bool Aligned > template< typename RealT >
RealT Hybrid_bitset_range_element<IndexT,Aligned>::gather( const RealT *weights, typename Hybrid_bitset_range_element<IndexT,Aligned>::index_t selection_mask )
{
	//return apegrunt::popcnt( selection_mask );
	RealT sum(0);

	while( selection_mask )
	{
		if(selection_mask & 1) { sum += *weights; }
		selection_mask = selection_mask >> 1;
		++weights;
	}

	return sum;
}

template< typename IndexT, bool Aligned > template< typename RealT >
RealT Hybrid_bitset_range_element<IndexT,Aligned>::gather( const RealT *weights, const Hybrid_bitset_range_element<IndexT,Aligned>& selection_range )
{
	//return selection_range.size();
	return std::accumulate( weights+selection_range.range_begin(), weights+selection_range.range_end(), RealT(0), []( auto sum, auto w ){ return sum+w; } );
}

template< typename IndexT >
inline std::ostream& operator<< ( std::ostream& os, const Hybrid_bitset_range_element<IndexT,true>& element )
{
	os << element();

	if( element.is_range() ) { os << "-" << element.range_end(); }
	else { os << ":" << element.get_bitfield(); }

	return os;
}

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

} // namespace apegrunt

#endif // APEGRUNT_INTEGERSEQUENCE_HYBRID_BITSET_RANGE_HPP
