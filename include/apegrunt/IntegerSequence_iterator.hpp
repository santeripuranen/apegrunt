/** @file IntegerSequence_iterator.hpp
 
	Copyright (c) 2016-2019 Santeri Puranen. All rights reserved.
 
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

    inline bool operator==( const my_type& rhs ) const { return rhs.m_pos == m_pos; } // this is wrong if *m_pos is not a simple integer
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
