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

#include "IntegerSequence_Hybrid_bitset_range_forward.h"

namespace apegrunt {

template< typename IndexT, bool Aligned >
struct alignas(sizeof(IndexT)*2) Hybrid_bitset_range_element
{
	using index_t = IndexT;
	using value_type = index_t;
	using my_type = Hybrid_bitset_range_element<index_t,Aligned>;
	using aligned = std::integral_constant<bool, Aligned>;

	static constexpr index_t create_range_flag_mask() { return index_t(1) << std::numeric_limits<index_t>::digits-1; } // set the most significant bit

	enum { BITSTRING_SIZE = std::numeric_limits<index_t>::digits };
	enum { RANGE_FLAG_MASK = my_type::create_range_flag_mask() };
	enum { RANGE_UNFLAG_MASK = ~my_type::create_range_flag_mask() };
	enum { MODULO_MASK=std::numeric_limits<index_t>::digits-1 }; // 0x7 for 8-bit, 0xF for 16-bit, 0x1F for 32-bit and 0x3F for 64-bit

	Hybrid_bitset_range_element( ) : m_pos(0), m_bitfield(0) { } // zero-initialize
	Hybrid_bitset_range_element( index_t value ) : m_pos(value & ~MODULO_MASK), m_bitfield(1 << (value & MODULO_MASK)) { }
	Hybrid_bitset_range_element( index_t value, index_t bitfield ) : m_pos(value), m_bitfield(bitfield) { } // trust that 'value' is properly aligned

	Hybrid_bitset_range_element( const my_type& other ) : m_pos(other.m_pos), m_bitfield(other.m_bitfield) { }

	inline my_type& operator=( const my_type& other ) { m_pos = other.m_pos; m_bitfield = other.m_bitfield; return *this; }
	inline my_type& operator=( index_t value ) { m_pos = value & ~MODULO_MASK; m_bitfield = 1 << (value & MODULO_MASK); return *this; }

	inline index_t operator()() const noexcept { return m_pos & RANGE_UNFLAG_MASK; }

	inline bool set( index_t pos ) noexcept
	{
		if( ( pos & ~MODULO_MASK ) == m_pos )
		{
			m_bitfield = m_bitfield | (1 << (pos & MODULO_MASK));
			return true;
		}
		return false;
	}
	inline bool is_set( index_t pos ) const noexcept
	{
		return (( pos & ~MODULO_MASK ) == m_pos) && ( m_bitfield & (1 << (pos & MODULO_MASK)) );
	}

	inline bool is_unset( index_t pos ) const noexcept { return !this->is_set(pos); }

	inline bool all_set() const noexcept { return !bool( ~m_bitfield ); }

	inline bool is_range() const noexcept { return bool( m_pos & RANGE_FLAG_MASK ); }
	inline void set_range() noexcept { m_pos = m_pos | RANGE_FLAG_MASK; }
	inline void set_range_end( index_t pos ) noexcept { this->set_range(); m_bitfield = pos; }

	inline index_t range_begin() const noexcept { return (*this)(); }
	inline index_t range_end() const noexcept { return this->is_range() ? m_bitfield : m_pos; }

	inline index_t get_bitfield() const noexcept { return m_bitfield; }

	index_t m_pos;
	index_t m_bitfield;

	template< typename RealT >
	static inline RealT gather( const RealT *weights, index_t selection_mask );
};

template< typename IndexT, bool Aligned > template< typename RealT >
RealT Hybrid_bitset_range_element<IndexT,Aligned>::gather( const RealT *weights, IndexT selection_mask )
{
	RealT sum(0);

	while( selection_mask )
	{
		if(selection_mask & 1) { sum += *weights; }
		selection_mask = selection_mask >> 1;
		++weights;
	}

	return sum;
}

template< typename IndexT >
inline std::ostream& operator<< ( std::ostream& os, const Hybrid_bitset_range_element<IndexT,true>& element )
{
	os << element();

	if( element.is_range() ) { os << "-" << element.range_end(); }
	else { os << ":" << element.get_bitfield(); }

	return os;
}

} // namespace apegrunt

#endif // APEGRUNT_INTEGERSEQUENCE_HYBRID_BITSET_RANGE_HPP
