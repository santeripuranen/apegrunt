/** @file IndexVector_interface.hpp
 
	Copyright (c) 2016 Santeri Puranen. All rights reserved.
 
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

#ifndef APEGRUNT_INDEXVECTOR_INTERFACE_HPP
#define APEGRUNT_INDEXVECTOR_INTERFACE_HPP

#include <string>
#include <iterator>
#include <limits>
#include <ostream>
#include <stdexcept>

#include <boost/serialization/serialization.hpp>
#include <boost/serialization/nvp.hpp>
#include <boost/range/iterator_range_core.hpp>
#include <boost/iterator.hpp>

#include "IndexVector_forward.h"

namespace apegrunt {

template< typename IndexT >
struct RLE_index
{
	RLE_index( IndexT value, IndexT count=1 ) : m_value(value), m_count(count) { }

	RLE_index<IndexT>& operator=( const RLE_index<IndexT>& other ) { m_value = other.m_value; m_count = other.m_count; return *this; }

	IndexT m_value;
	IndexT m_count;
};

namespace iterator {

template< typename IndexT >
struct IndexVector_const_iterator
{
	using value_type = IndexT;
	using reference = const value_type&; // almost always T& or const T&
    using pointer = const value_type*; //almost always T* or const T*
	using iterator_category = std::random_access_iterator_tag;
	//using iterator_category = std::forward_iterator_tag;
	using difference_type = std::ptrdiff_t; //almost always ptrdiff_t

	using my_type = IndexVector_const_iterator<IndexT>;
	IndexVector_const_iterator() = delete;
	IndexVector_const_iterator( RLE_index<IndexT> const *pos, RLE_index<IndexT> const * const end ) : m_pos(pos), m_count(end==pos ? 0 : pos->m_count-1), m_end(end) { }
	~IndexVector_const_iterator() = default;

	IndexVector_const_iterator( const my_type& other ) : m_pos(other.m_pos), m_count(other.m_count), m_end(other.m_end) { }
	IndexVector_const_iterator( my_type&& other ) noexcept : m_pos(std::move(other.m_pos)), m_count(std::move(other.m_count)), m_end(std::move(other.m_end)) { }

    inline my_type& operator++()
    { /*std::cout << "increment" << std::endl;*/
    	m_count = ( 0 < m_count ? m_count-1 : (m_end == m_pos ? 0 : (m_end == ++m_pos ? 0 : m_pos->m_count-1 ) ) );
       	//if( 0 < m_count ) { --m_count; } else { m_count = (m_end == ++m_pos ? 0 : m_pos->m_count-1 ); }
    	return *this;
    }
    inline my_type operator++(int) { my_type tmp(*this); operator++(); return tmp; }
    //inline my_type& operator+=( std::size_t n ) { return this->addition_assignment_operator_impl( n ); }
    inline reference operator*() const { /*std::cout << "deref" << std::endl;*/ return m_pos->m_value; }
    inline pointer operator->() { return &(*this); }
    inline bool operator==( const my_type& rhs ) const { /*std::cout << "operator==(" << m_pos << "," << rhs.m_pos << "): " << (rhs.m_pos == m_pos ? "true" : "false") << std::endl;*/ return rhs.m_pos == m_pos && rhs.m_count == m_count; }
    inline bool operator<( const my_type& rhs ) const { /*std::cout << "operator<(" << m_pos << "," << rhs.m_pos << "): " << (m_pos < rhs.m_pos ? "true" : "false") << std::endl;*/ return m_pos < rhs.m_pos; }

    // random-access iterator requirements

    // operator- is not a constant time op, but otherwise behaves
    // as per the requirements for a random_access_iterator
    inline difference_type operator-( const my_type& rhs ) const
    {
    	std::cout << "IndexVector_const_iterator.operator-()" << std::endl;
    	difference_type diff = 0;
    	if( this->m_pos < rhs.m_pos )
    	{
    		auto pos = this->m_pos;
    		++pos;
    		while( rhs.m_pos != pos )
    		{
    			diff -= pos.m_count;
    			++pos;
    		}
    	}
    	else if( this->m_pos > rhs.m_pos )
    	{
    		auto pos = rhs.m_pos;
    		++pos;
    		while( this->m_pos != pos )
    		{
    			diff += pos.m_count;
    			++pos;
    		}
    	}
    	// this->mpos == rhs.m_pos
    	diff += this->m_count - rhs.m_count;

		return diff;
    }

    RLE_index<IndexT> const *m_pos;
    RLE_index<IndexT> const * const m_end;
    mutable IndexT m_count;
};

template< typename IndexT >
inline bool operator!=( const IndexVector_const_iterator<IndexT>& lhs, const IndexVector_const_iterator<IndexT>& rhs ) { return !(lhs == rhs); }
template< typename IndexT >
inline bool operator>( const IndexVector_const_iterator<IndexT>& lhs, const IndexVector_const_iterator<IndexT>& rhs ) { return (rhs < lhs); }
template< typename IndexT >
inline bool operator<=( const IndexVector_const_iterator<IndexT>& lhs, const IndexVector_const_iterator<IndexT>& rhs ) { return !(lhs > rhs); }
template< typename IndexT >
inline bool operator>=( const IndexVector_const_iterator<IndexT>& lhs, const IndexVector_const_iterator<IndexT>& rhs ) { return !(lhs < rhs); }

} // namespace iterator

template< typename IndexT >
inline typename iterator::IndexVector_const_iterator<IndexT>::difference_type operator-( const iterator::IndexVector_const_iterator<IndexT>& lhs, const iterator::IndexVector_const_iterator<IndexT>& rhs )
{
	return lhs.operator-( rhs );
}

template< typename IndexT >
inline typename iterator::IndexVector_const_iterator<IndexT>::difference_type distance( const iterator::IndexVector_const_iterator<IndexT>& lhs, const iterator::IndexVector_const_iterator<IndexT>& rhs )
{
	return rhs - lhs;
}

template< typename IndexT >
class IndexVector
{

public:
	using value_type = IndexT;
	using const_iterator = apegrunt::iterator::IndexVector_const_iterator<IndexT>;

	IndexVector(): m_size(0), m_cached(std::cend(m_internal_storage)), m_dirty(true) { } // = default;
	~IndexVector() = default;

	//static constexpr std::size_t max() { return s_mask; }

	inline std::size_t size() const { return m_size; }
	inline std::size_t bytesize() const { return m_internal_storage.size()*sizeof(typename internal_storage_t::value_type); }

	inline value_type operator[]( std::size_t index ) const
	{
// /*
		if( index > this->size() )
		{
			throw std::out_of_range("IndexVector out_of_range access");
		}

		if( m_dirty || index < m_cached_pos || m_cached->m_count <= index-m_cached_pos )
		{
			m_cached_pos = 0;
			for( auto token = std::cbegin(m_internal_storage); token != std::cend(m_internal_storage); ++token )
			{
				//m_cached_pos += token.m_count;
				if( index < m_cached_pos+token->m_count ) { m_cached = token; /*std::cout << "operator[" << index << "] token=(" << token.m_value << "," << token.m_count << ")\n";*/ break; }
				m_cached_pos += token->m_count;
			}
			m_dirty = false;
		}
		//std::cout << "operator[" << index << "]=(" << m_cached->m_value << "," << m_cached->m_count << ")\n";  std::cout.flush();
		return m_cached->m_value;

// */
/*
		if( index < this->size() )
		{
			value_type pos = 0;
			for( const auto& token: m_internal_storage )
			{
				pos += token.m_count;
				if( index < pos ) { return token.m_value; }
			}
		}
*/
/*
		if( index < this->size() )
		{
			//std::cout << "index=" << index << std::endl;
			std::size_t pos=0;
			for( auto i: *this )
			{
				//std::cout << "  pos=" << pos << " m_size=" << m_size << std::endl;
				if( pos == index )
				{
					//std::cout << "  break&return=" << i << "\n" << std::endl;
					return i;
				}
				++pos;
			}
			//std::cout << "We went too faaar" << std::endl;
		}
*/
	}

	template< typename T >
	inline void push_back( T value, T count=1 )
	{
		m_dirty = true;
		if( m_size != 0 && value == m_internal_storage.back().m_value )
		{
			if( T(s_maxval) < T(m_internal_storage.back().m_count) + count )
			{
		    	const auto a = std::div( int64_t(m_internal_storage.back().m_count + count), int64_t(s_maxval) );
				m_internal_storage.back().m_count = s_maxval;
				for( std::size_t i=a.quot-1; i != 0; --i ) { m_internal_storage.emplace_back( IndexT(value), IndexT(s_maxval) ); }
				if( a.rem != 0 ) { m_internal_storage.emplace_back( IndexT(value), IndexT(a.rem) ); }
			}
			else
			{
				m_internal_storage.back().m_count += IndexT(count);
			}
/*
			auto free_capacity = s_maxval-m_internal_storage.back().m_count;
			//std::cout << "s_maxval=" << s_maxval << " m_count=" << m_internal_storage.back().m_count << " free_capacity=" << free_capacity << "\n";

			while( free_capacity < IndexT(count) )
			{

				m_internal_storage.back().m_count = s_maxval;
				m_internal_storage.emplace_back( IndexT(value), IndexT( count-free_capacity > s_maxval ? s_maxval : count ) );
				count = ( count-free_capacity > s_maxval ? count-free_capacity : 0 );
				free_capacity = s_maxval-m_internal_storage.back().m_count;
				//std::cout << "  insert loop: m_count=" << m_internal_storage.back().m_count << "\n";
				std::cout << "m_value=" << m_internal_storage.back().m_value << " m_count=" << m_internal_storage.back().m_count << " count=" << count << " free_capacity=" << free_capacity << std::endl;
			}
			m_internal_storage.back().m_count += IndexT(count);
*/
		}
		else
		{
			//m_dirty = true;
			//std::cout << "new count=" << count << " s_maxval=" << s_maxval << std::endl;
			//m_internal_storage.push_back( RLE_index<IndexT>(IndexT(value),IndexT(count)) );
			m_internal_storage.emplace_back( IndexT(value), IndexT(count) );
		}
		m_size += count;
		//std::cout << "  done" << std::endl;
	}

	inline void reserve( std::size_t size ) { } // do nothing

	inline void clear() { m_internal_storage.clear(); m_size = 0; }
/*
	inline std::ostream& print( std::ostream& os ) const
	{
		os << "IndexVector internal:";
		for( const auto& pair: m_internal_storage ) { os << " " << pair.m_value << ":" << pair.m_count; }
		os << "\n            external:";
		for( auto i: *this ) { os << " " << i; }
		os << "\n";
		return os;
	}
*/
	inline const_iterator cbegin() const { return const_iterator( m_internal_storage.data(), this->one_past_end() ); }
	inline const_iterator begin() const { return this->cbegin(); }
	inline const_iterator cend() const { return const_iterator( this->one_past_end(), this->one_past_end() ); }
	inline const_iterator end() const { return this->cend(); }

private:
	using internal_storage_t = std::vector< RLE_index< IndexT > >;

	RLE_index<IndexT> const * const one_past_end() const { return &(m_internal_storage.back())+1; }

	//static constexpr IndexT s_mask = std::numeric_limits<IndexT>::max() >> 1;
	//static constexpr IndexT s_halfway = std::numeric_limits<IndexT>::max() / 2;

	static constexpr IndexT s_maxval = std::numeric_limits<IndexT>::max();
	internal_storage_t m_internal_storage;
	std::size_t m_size;
	mutable typename internal_storage_t::const_iterator m_cached;
	mutable IndexT m_cached_pos;
	mutable bool m_dirty;
};

template< typename IndexT >
std::ostream& operator<< ( std::ostream& os, const IndexVector<IndexT>& index_vector )
{
	//index_vector.print(os);
	for( auto i: index_vector ) { os << " " << i; }
	return os;
}

template< typename IndexT >
typename IndexVector<IndexT>::const_iterator cbegin( const IndexVector<IndexT>& index_vector ) { return index_vector.cbegin(); }

template< typename IndexT >
typename IndexVector<IndexT>::const_iterator cend( const IndexVector<IndexT>& index_vector ) { return index_vector.cend(); }

template< typename IndexT >
typename IndexVector<IndexT>::const_iterator begin( const IndexVector<IndexT>& index_vector ) { return index_vector.cbegin(); }

template< typename IndexT >
typename IndexVector<IndexT>::const_iterator end( const IndexVector<IndexT>& index_vector ) { return index_vector.cend(); }

template< typename T >
std::size_t bytesize( const IndexVector<T>& v ) { return v.bytesize(); }

} // namespace apegrunt

#endif // APEGRUNT_INDEXVECTOR_INTERFACE_HPP
