/** @file IntegerSequence_Hybrid_bitset_range_iterator.hpp
 
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

#ifndef APEGRUNT_INTEGERSEQUENCE_HYBRID_BITSET_RANGE_ITERATOR_HPP
#define APEGRUNT_INTEGERSEQUENCE_HYBRID_BITSET_RANGE_ITERATOR_HPP

#include <vector>
#include <iterator>

#include "IntegerSequence_iterator.hpp"
#include "IntegerSequence_Hybrid_bitset_range.hpp"

namespace apegrunt {

namespace iterator {

// A reasonably fast version of the value_type::BITSTRING_SIZE aligned container iterator
template< typename IndexT >
struct IntegerSequence_const_iterator< Hybrid_bitset_range_element<IndexT,true> >
{
	using index_t = IndexT;
	using range_type = Hybrid_bitset_range_element<index_t,true>;

	// std::iterator requirements
	using value_type = typename range_type::index_t;
	using reference = const value_type&; // only const access, since value_type is a proxy
	using pointer = const value_type*; // only const access, since value_type is a proxy
	//using iterator_category = std::random_access_iterator_tag;
	using iterator_category = std::forward_iterator_tag;
	using difference_type = std::ptrdiff_t; //almost always ptrdiff_t

	using my_type = IntegerSequence_const_iterator<range_type>;

	IntegerSequence_const_iterator() = delete; // disallow default-constructed iterators
	IntegerSequence_const_iterator( const range_type *pos, const range_type * const end )
	: m_pos(pos), m_end(end), m_currpos(0), m_mode(0)
	{
		m_mode = (m_pos != m_end) ? ( m_pos->is_range() ? 1 : 2 ) : 0;
		if( m_mode )
		{
			m_currpos = (*m_pos)(); m_bitstring = m_pos->get_bitfield();
			if( m_mode == 2 && !(m_bitstring & 1) )
			{
				// fast-forward to first index in set
				while( m_bitstring && !this->advance_and_test() ) { }
			}
		}
	}
	~IntegerSequence_const_iterator() = default;

	inline reference operator*() { return m_currpos; }
    inline pointer operator->() { return &m_currpos; }
	inline void operator++()
	{
		switch(m_mode)
		{
		case 1: // interval
			++m_currpos == m_bitstring && (m_mode = 3);
			break;

		case 2: // bitstring/mask
			while( m_bitstring && !this->advance_and_test() ) { }
			if( m_bitstring ) { break; }

		case 3: // previous bitstring/interval has ended; start a new cycle
			if( ++m_pos != m_end )
			{
				m_currpos = (*m_pos)();
				m_bitstring = m_pos->get_bitfield();
				//m_mode = m_pos->is_range() ? 1 : 2;
				//if( m_mode == 2 && !(m_bitstring & 1) )
				if( (m_mode = m_pos->is_range() ? 1 : 2) == 2 && !(m_bitstring & 1) )
				{
					// fast-forward to first index in set
					while( m_bitstring && !this->advance_and_test() ) { }
				}
			}
			else
			{
				m_mode = 0;
				m_currpos = 0;
			}
			break;

		default: // end has been reached; do nothing
			(m_pos != m_end) && ( m_pos = m_end );
			break;
		}
	}

    inline bool operator==( const my_type& rhs ) const { return rhs.m_currpos == m_currpos && rhs.m_pos == m_pos; }
    inline bool operator!=( const my_type& rhs ) const { return !(rhs == *this); }
    inline bool operator<( const my_type& rhs ) const { return m_currpos < rhs.m_currpos; }
    inline bool operator>( const my_type& rhs ) const { return m_currpos > rhs.m_currpos; }

    inline bool advance_and_test() noexcept { ++m_currpos; return ( m_bitstring = m_bitstring >> 1 ) & 1; }

	index_t m_currpos;
	index_t m_bitstring;
	uint_fast8_t m_mode;
    const range_type *m_pos;
	const range_type * const m_end;
};

// specialize back_insert_iterator for Hybrid_bitset_range_element

// This version is aligned at value_type::BITSTRING_SIZE chunks, and does not
// have single index entries, only bitstrings and interval entries.
template< typename IndexT >
struct back_insert_iterator< std::vector< Hybrid_bitset_range_element<IndexT,true> > >
{
	using container_t = std::vector< Hybrid_bitset_range_element<IndexT,true> >;
	using value_type = typename container_t::value_type;
	using index_t = typename value_type::index_t;
	using my_type = back_insert_iterator<container_t>;

	back_insert_iterator() = delete;
	back_insert_iterator( container_t& container ) : m_container(container), m_bitfield_mode(false), m_is_contiguous(false) { }
	~back_insert_iterator() = default;

	constexpr back_insert_iterator( const my_type& other )
	: m_container( other.m_container ),
	  m_bitfield_mode( other.m_bitfield_mode ),
	  m_is_contiguous( other.m_is_contiguous )
	{ }

	constexpr back_insert_iterator( const my_type& other, container_t& new_container )
	: m_container( new_container ),
	  m_bitfield_mode( other.m_bitfield_mode ),
	  m_is_contiguous( other.m_is_contiguous )
	{ }

	constexpr back_insert_iterator( my_type&& other )
	: m_container( other.m_container ), // we don't want to move the container, since it's not owned by us
	  m_bitfield_mode( std::move(other.m_bitfield_mode) ),
	  m_is_contiguous( std::move(other.m_is_contiguous) )
	{ }

	inline my_type& operator=( const index_t& value )
	{
		switch(m_bitfield_mode)
		{
		case true: // fill bitset
			if( m_container.back().set(value) )
			{
				if( m_container.back().all_set() ) // is our bitset a contiguous range (from m_container.back()() to value)?
				{
					if( m_is_contiguous ) // were we already previously in a contiguous range?
					{
						// ..then remove the bitset element and extend the range (below)
						m_container.pop_back();
					}
					else // we were not previously in a contiguous range
					{
						// ..but we are now!
						m_is_contiguous = true;
					}
					m_container.back().set_range_end(value); // update range end value
					m_bitfield_mode = false; // switch back to new bitfield insertion mode
				}
				break;
			}
			m_is_contiguous = false; // we can no longer extend a contiguous range
			// fall through to default

		default: // insert new bitset range element
			m_container.push_back(value); // range begin value

			m_bitfield_mode = true; // switch to bitset fill mode
			break;
		}

		return *this;
	}

	container_t& m_container;
	bool m_bitfield_mode;
	bool m_is_contiguous;
};

} // namespace iterator

} // namespace apegrunt

#endif // APEGRUNT_INTEGERSEQUENCE_HYBRID_BITSET_RANGE_ITERATOR_HPP
