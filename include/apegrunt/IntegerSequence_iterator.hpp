/** @file IntegerSequence_iterator.hpp
 
	Copyright (c) 2016-2020 Santeri Puranen.

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

#ifndef APEGRUNT_INTEGERSEQUENCE_ITERATOR_HPP
#define APEGRUNT_INTEGERSEQUENCE_ITERATOR_HPP

#include <iterator>

#include "IntegerSequence_forward.h"

namespace apegrunt {

namespace iterator {

template< typename RangeT >
struct IntegerSequence_const_iterator
{
	using range_type = RangeT;

	// std::iterator requirements
	using value_type = typename range_type::index_t;
	using reference = const value_type&;
	using pointer = const value_type*;
	using iterator_category = std::random_access_iterator_tag;
	//using iterator_category = std::forward_iterator_tag;
	using difference_type = std::ptrdiff_t; //almost always ptrdiff_t

	using my_type = IntegerSequence_const_iterator<range_type>;

	IntegerSequence_const_iterator() = delete; // disallow default-constructed iterators
	IntegerSequence_const_iterator( const range_type *pos, const range_type * const end ) : m_pos(pos), m_currpos( pos != end ? (*pos)() : 0 ), m_end(end) { }
	~IntegerSequence_const_iterator() = default;

	inline reference operator*() { return m_currpos; }
    inline pointer operator->() { return &m_currpos; }
	inline void operator++() { m_currpos = (*(++m_pos))(); }

    inline bool operator==( const my_type& rhs ) const { return rhs.m_pos == m_pos; }
    inline bool operator<( const my_type& rhs ) const { return m_currpos < rhs.m_currpos; }

	const range_type *m_pos;
	value_type m_currpos;
	const range_type *m_end;
};

template< typename RangeT >
inline bool operator!=( const IntegerSequence_const_iterator<RangeT>& lhs, const IntegerSequence_const_iterator<RangeT>& rhs ) { return !(lhs == rhs); }
template< typename RangeT >
inline bool operator>( const IntegerSequence_const_iterator<RangeT>& lhs, const IntegerSequence_const_iterator<RangeT>& rhs ) { return (rhs < lhs); }
template< typename RangeT >
inline bool operator<=( const IntegerSequence_const_iterator<RangeT>& lhs, const IntegerSequence_const_iterator<RangeT>& rhs ) { return !(lhs > rhs); }
template< typename RangeT >
inline bool operator>=( const IntegerSequence_const_iterator<RangeT>& lhs, const IntegerSequence_const_iterator<RangeT>& rhs ) { return !(lhs < rhs); }

template< typename ContainerT>
struct back_insert_iterator
{
	using container_t = ContainerT;
	using value_type = typename container_t::value_type;
	using my_type = back_insert_iterator<container_t>;

	back_insert_iterator() = delete;
	back_insert_iterator( container_t& container ) : m_container(container) { }
	back_insert_iterator( const my_type& other ) : m_container( other.m_container ) { }

	constexpr back_insert_iterator( const my_type& other, container_t& new_container )
	: m_container( new_container )
	{
	}

	inline my_type& operator=( const value_type& value ) { m_container.push_back(value); return *this; }

	container_t& m_container;
};

} // namespace iterator

template< typename ContainerT >
apegrunt::iterator::back_insert_iterator<ContainerT> back_inserter( ContainerT& container ) { return apegrunt::iterator::back_insert_iterator<ContainerT>( container ); }

} // namespace apegrunt

#endif // APEGRUNT_INTEGERSEQUENCE_ITERATOR_HPP
