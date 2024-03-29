/** @file IntegerSequence_Hybrid_bitset_range_iterator.hpp
 
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
	using value_type = typename range_type::value_type;
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
		m_mode = (m_pos != m_end) ? m_pos->is_range() : 2;
		if( m_mode != 2 )
		{
			m_currpos = (*m_pos)(); m_bitstring = m_pos->get_bitfield();
			if( !m_mode && !(m_bitstring & 1) )
			{
				// fast-forward to first index in set
				//while( m_bitstring && !this->advance_and_test() ) { }
				this->find_next();
			}
		}
	}
	~IntegerSequence_const_iterator() = default;

	inline my_type& operator=( const my_type& other )
	{
		m_pos = other.m_pos;
		m_end = other.m_end;
		m_currpos = other.m_currpos;
		m_mode = other.m_mode;

		return *this;
	}

	inline reference operator*() const { return m_currpos; }
    inline pointer operator->() { return &m_currpos; }
	inline void operator++()
	{
		switch(m_mode)
		{
		case 0: // bitstring/mask
			if( !this->find_next() ) { this->next_block(); }
			break;

		case 1: // interval
			if( !((++m_currpos >> range_type::MODULO_NSTEPS) < m_bitstring) ) { this->next_block(); }
			break;

		default: // end has been reached; do nothing
			(m_pos != m_end) && (m_pos = m_end);
			break;
		}
	}

	inline void next_block()
	{
		if( ++m_pos != m_end )
		{
			m_currpos = (*m_pos)();
			m_bitstring = m_pos->get_bitfield();
			m_mode = m_pos->is_range();
			if( !m_mode && !(m_bitstring & 1) ) { this->find_next(); } // next_block() recursion not required, since find_next() should never fail here (unless a block is malformed)
		}
		else
		{
			m_mode = 2;
			m_currpos = 0;
			m_bitstring = 0;
		}
	}

    inline bool operator==( const my_type& rhs ) const { return rhs.m_currpos == m_currpos && rhs.m_pos == m_pos; }
    inline bool operator!=( const my_type& rhs ) const { return !(rhs == *this); }
    inline bool operator<( const my_type& rhs ) const { return m_currpos < rhs.m_currpos; }
    inline bool operator>( const my_type& rhs ) const { return m_currpos > rhs.m_currpos; }

    inline bool advance_and_test() noexcept { ++m_currpos; return (m_bitstring >>= 1) & 1; }
    inline bool find_next() noexcept { do { ++m_currpos; } while( m_bitstring && !( (m_bitstring >>= 1) & 1) ); return m_bitstring; }

    const range_type *m_pos;
	const range_type *m_end;
	value_type m_currpos;
	index_t m_bitstring;
	uint_fast8_t m_mode;
};

// specialize back_insert_iterator for Hybrid_bitset_range_element

// This version is aligned at value_type::BITSTRING_SIZE chunks, and does not
// have single index entries, only bitstrings and interval entries.
template< typename IndexT, typename AllocatorT >
struct back_insert_iterator< std::vector< Hybrid_bitset_range_element<IndexT,true>, AllocatorT > >
{
	using value_type = Hybrid_bitset_range_element<IndexT,true>;
	using allocator_t = AllocatorT;
	using container_t = std::vector< value_type, allocator_t >;
	using index_t = typename value_type::value_type;
	using my_type = back_insert_iterator<container_t>;

	back_insert_iterator() = delete;
	back_insert_iterator( container_t& container ) : m_container(container), m_bitfield_mode(false), m_is_contiguous(false) { }
	~back_insert_iterator() = default;

	constexpr back_insert_iterator( const my_type& other )
	: m_container( other.m_container ),
	  m_bitfield_mode( other.m_bitfield_mode ),
	  m_is_contiguous( other.m_is_contiguous )
	{ }

	constexpr back_insert_iterator( my_type&& other )
	: m_container( other.m_container ), // we don't want to move the container, since we only hold it by reference
	  m_bitfield_mode( std::move(other.m_bitfield_mode) ),
	  m_is_contiguous( std::move(other.m_is_contiguous) )
	{ }

	constexpr back_insert_iterator( const my_type& other, container_t& new_container )
	: m_container( new_container ),
	  m_bitfield_mode( other.m_bitfield_mode ),
	  m_is_contiguous( other.m_is_contiguous )
	{ }

	inline my_type& operator=( const my_type& other )
	{
		m_container = other.m_container;
		m_bitfield_mode = other.m_bitfield_mode;
		m_is_contiguous = other.m_is_contiguous;
		return *this;
	}

	inline my_type& operator=( my_type&& other )
	{
		m_container = other.m_container;
		m_bitfield_mode = std::move(other.m_bitfield_mode);
		m_is_contiguous = std::move(other.m_is_contiguous);
		return *this;
	}

	inline my_type& reset( const my_type& other, container_t& new_container )
	{
		m_container = new_container;
		m_bitfield_mode = other.m_bitfield_mode;
		m_is_contiguous = other.m_is_contiguous;
		return *this;
	}

	inline my_type& reset()
	{
		m_bitfield_mode = true;
		m_is_contiguous = m_is_contiguous && m_container.back().is_range();
		return *this;
	}

	inline my_type& operator=( index_t value )
	{
		switch(m_bitfield_mode)
		{
		case true: // fill bitset
			if( m_container.back().set(value) )
			{
				if( m_container.back().all_set() ) // is our bitset a contiguous range (from m_container.back()() to value)?
				{
					if( m_is_contiguous ) // are we already in a contiguous range?
					{
						// ..then remove the bitset element and extend the existing range (below)
						m_container.pop_back();
					}

					m_container.back().set_range_end(value+1); // end == past-the-end
					m_bitfield_mode = false; // switch back to new bitfield insertion mode
					m_is_contiguous = true; // all bits are set, so we're now in contiguous range mode
				}
				break;
			}
			m_is_contiguous = false; // we can no longer extend a contiguous range
			// no break == fall through to default

		default: // insert new bitset range element
			m_bitfield_mode = true; // switch to bitset fill mode
			m_is_contiguous = ( m_is_contiguous && (value == m_container.back().range_end()) ); // should we try to extend a previous range once the current bitfield is full?
			m_container.emplace_back(value); // range begin value

			break;
		}

		return *this;
	}

	inline my_type& operator=( const value_type& value )
	{
		if( m_is_contiguous && value.is_range() && m_container.back().range_end() >= value.range_begin() ) { m_container.back().set_range_end(value.range_end()); }
		else { m_container.emplace_back(value); }

		m_is_contiguous = value.is_range();
		m_bitfield_mode = m_is_contiguous; // put next single value through default pathway

		return *this;
	}

	inline my_type& operator=( value_type&& value )
	{
		if( m_is_contiguous && value.is_range() && m_container.back().range_end() >= value.range_begin() ) { m_container.back().set_range_end(value.range_end()); /*m_is_contiguous = true; // no need to set, as it already is true */ }
		else { m_is_contiguous = value.is_range(); m_container.emplace_back(std::move(value)); }

		m_bitfield_mode = m_is_contiguous; // put next single value through default pathway
		return *this;
	}

	container_t& m_container;
	bool m_bitfield_mode;
	bool m_is_contiguous;
};

} // namespace iterator

} // namespace apegrunt

#endif // APEGRUNT_INTEGERSEQUENCE_HYBRID_BITSET_RANGE_ITERATOR_HPP
