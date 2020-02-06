/** @file IntegerSequence_interface.hpp
 
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

#ifndef APEGRUNT_INTEGERSEQUENCE_INTERFACE_HPP
#define APEGRUNT_INTEGERSEQUENCE_INTERFACE_HPP

#include <vector>
#include <numeric> // for std::accumulate

#include <boost/container_hash/hash.hpp>
#include <boost/dynamic_bitset.hpp>
//#define BOOST_POOL_NO_MT
//#include <boost/pool/pool_alloc.hpp>

#include "IntegerSequence_forward.h"
#include "IntegerSequence_operations_forward.h"
#include "IntegerSequence_iterator.hpp"
#include "aligned_allocator.hpp"

namespace apegrunt {

template< typename RangeT >
class IntegerSequence
{
public:
	using range_t = RangeT;
	using my_type = IntegerSequence<range_t>;
	using value_type = typename range_t::value_type;
	using aligned = typename range_t::aligned;
	using const_iterator = apegrunt::iterator::IntegerSequence_const_iterator<range_t>;
	using iterator = const_iterator;

	using container_t = std::vector< range_t >;

	IntegerSequence() : m_back_inserter( m_storage ), m_size(0), m_hash(0), m_block_mask(0) { }

	IntegerSequence( const value_type& value ) : m_back_inserter( m_storage ), m_size(0), m_hash(0), m_block_mask(0) { this->push_back(value); }

	IntegerSequence( const my_type& other )
	: m_storage(other.m_storage),
	  m_back_inserter(other.m_back_inserter, m_storage), // need to pass reference to the newly _copied_ m_storage when copy-constructing back_inserter. The regular back_inserter copy constructor would copy a reference to the old m_storage. This is an issue for example when reallocating std::vector<IntegerSequence>.
	  m_size(other.m_size),
	  m_hash(other.m_hash),
	  m_block_mask(other.m_block_mask)
	{ }

	IntegerSequence( my_type&& other )
	: m_storage( std::move(other.m_storage) ),
	  m_back_inserter( other.m_back_inserter, m_storage ), // reference to the newly _moved_ m_storage; see comment above
	  m_size( std::move(other.m_size) ),
	  m_hash( std::move(other.m_hash) ),
	  m_block_mask( std::move(other.m_block_mask) )
	{ }

	template< typename BlockT >
	IntegerSequence( const boost::dynamic_bitset<BlockT>& bitset )
	: m_back_inserter( m_storage ), m_size(0), m_hash(0), m_block_mask(0)
	{
		using dynamic_bitset_t = boost::dynamic_bitset<BlockT>;
		// we could do better here by using block access to the underlying storage of dynamic_bitset
		auto pos = bitset.find_first();
		while( pos != dynamic_bitset_t::npos )
		{
			this->push_back(pos);
			pos = bitset.find_next(pos);
		}
	}

	~IntegerSequence() // = default;
	{
		//boost::singleton_pool<boost::pool_allocator_tag, sizeof(range_t)>::purge_memory();
	}

	inline my_type& operator=( const my_type& other )
	{
		m_storage = other.m_storage; // copy
		m_back_inserter.reset(other.m_back_inserter,m_storage); // reference to the newly _copied_ m_storage; see comment above
		m_size = other.m_size;
		m_hash = other.m_hash;
		m_block_mask = other.m_block_mask;
		return *this;
	}

	inline my_type& operator=( my_type&& other )
	{
		m_storage = std::move(other.m_storage);
		m_back_inserter.reset(other.m_back_inserter,m_storage); // reference to the newly _moved_ m_storage; see comment above
		m_size = std::move(other.m_size);
		m_hash = std::move(other.m_hash);
		m_block_mask = std::move(other.m_block_mask);
		return *this;
	}

	inline void push_back( const value_type& value ) { m_back_inserter = value; ++m_size; this->update_block_mask( my_type::generate_block_mask(value) ); }

	// merge() invalidates iterators, pointers, references to elements of 'this'; 'other' remains intact
	inline my_type& merge( const my_type& other ) { *this = std::move( apegrunt::set_union( *this, other ) ); return *this; }

	// Quick check for potential overlap. May produce false positives, i.e. might not actually overlap when 'true', but is guaranteed _not_ to overlap if 'false'.
	inline bool has_overlap( const my_type& other ) const { return m_block_mask & other.m_block_mask; }

	// An IntegerSequence is empty if its internal storage contains nothing.
	// More consistent performance (O(1)) than "size() != 0", as is_empty() has no side effects.
	inline bool is_empty() const { return m_storage.empty(); }

	// check if the container contains the specified value
	inline bool contains( const value_type& value ) const
	{
		// our storage is 'sparse', so can't use direct access
		// this is a safe but horrible worst case O(N) solution
		for( auto& e: m_storage ) {	if( e.holds(value) ) { return true; } }
		return false;
		// however, our storage is sorted, so we really should use binary search, if only I could get the logic right
		//const auto itr = apegrunt::binary_search( m_storage.begin(), m_storage.end(), range_t(value), [](const range_t& lhs, const range_t& rhs){ return lhs < rhs; } );
		//return itr != m_storage.end() && itr->holds(value);
	}

	// Return number of stored values. Best case performance is O(1), worst case O(N), where N is the size of the internal container.
	inline std::size_t size() const { /*if( !m_storage.empty() && !m_size ) */ { this->update_size(); } return m_size; }
	inline std::size_t storagesize() const { return m_storage.size(); }

	//inline double fill_ratio() const { return double(m_storage.size(); } // not sure what I was planning to report here
	inline std::size_t bytesize() const { return sizeof(range_t) * m_storage.size(); }

	inline std::size_t hash() const { return m_hash == 0 ? this->cache_the_hash() : m_hash; } // return cached hash, creating the hash value as needed

	const_iterator cbegin() const { return const_iterator( m_storage.data(), this->one_past_end() ); }
	const_iterator cend() const { return const_iterator( this->one_past_end(), this->one_past_end() ); }

	iterator begin() const { return this->cbegin(); }
	iterator end() const { return this->cend(); }

	const container_t& storage() const { return m_storage; }

private:
	using block_mask_type = uint64_t;
	/*
		using allocator_t = boost::pool_allocator<range_t,
		    boost::default_user_allocator_new_delete,
		    boost::details::pool::default_mutex,
		    0xFF, 0xFFFF>;
	*/
	// the aligned allocator can be slow -- don't use it unless it's really needed (= only if aligned loads/stores absolutely required)
	//using allocator_t = memory::AlignedAllocator<range_t,alignof(range_t)>;
	//using container_t = std::vector< range_t, allocator_t >;
	container_t m_storage;
	apegrunt::iterator::back_insert_iterator< container_t > m_back_inserter;
	mutable std::size_t m_size;
	mutable std::size_t m_hash;
	block_mask_type m_block_mask;

	inline std::size_t cache_the_hash() const
	{
		m_hash = boost::hash_range( this->cbegin(), this->cend() );
		return m_hash;
	}
// /*
	// single element block mask
	static constexpr const block_mask_type generate_block_mask( value_type value ) noexcept { return block_mask_type(1) << (value >> range_t::BLOCK_NSTEPS); }

	// range element block mask
	static constexpr const block_mask_type generate_block_mask( value_type begin, value_type end ) noexcept
	{
		block_mask_type mask(0);
		for( auto shift(begin >> range_t::BLOCK_NSTEPS); shift < (end >> range_t::BLOCK_NSTEPS)+1; ++shift )
		{
			mask = mask | (block_mask_type(1)<<shift);
		}
		return mask;
	}

	static constexpr const block_mask_type generate_block_mask( const range_t& element )
	{
		return element.is_range() ? my_type::generate_block_mask( element.range_begin(), element.range_end() ) : my_type::generate_block_mask( element() );
	}
// */
	// beware O(N)
	inline void update_size() const
	{
		m_size = std::accumulate( std::cbegin(m_storage), std::cend(m_storage), 0, [](std::size_t sum, const auto& a) { return sum+a.size(); } );
	}

	//inline void update_block_mask( value_type value ) { m_block_mask = m_block_mask | my_type::generate_block_mask(value); }
	inline void update_block_mask( block_mask_type mask ) { m_block_mask = m_block_mask | mask; }

	inline const range_t * const one_past_end() const { return &(m_storage.back())+1; } // this pointer should never be dereferenced

	friend std::ostream& operator<< <value_type> ( std::ostream& os, const my_type& container );
	//friend my_type set_intersection<value_type,aligned::value>( const my_type& a, const my_type& b );
	friend my_type set_intersection<value_type>( const my_type& a, const my_type& b );
	//friend double intersect_and_gather<value_type,double,aligned::value>( const my_type& a, const my_type& b, const std::vector<double>& weights );
	friend double intersect_and_gather<value_type,double>( const my_type& a, const my_type& b, const std::vector<double>& weights );
	friend my_type set_union<value_type>( const my_type& a, const my_type& b );
	friend my_type set_union<value_type>( const std::vector< std::reference_wrapper< const my_type > >& sets );

};

template< typename RangeT >
bool operator==( const IntegerSequence< RangeT >& lhs, const IntegerSequence< RangeT >& rhs )
{
	// _if_ hash values are already cached (and they provide low collision probability), then comparing them first could be faster
	//return /*lhs.hash() == rhs.hash()*/ && std::equal( lhs.cbegin(), lhs.cend(), rhs.cbegin();
	return lhs.storagesize() == rhs.storagesize() && std::equal( lhs.cbegin(), lhs.cend(), rhs.cbegin() );
}

template< typename RangeT >
struct hash< apegrunt::IntegerSequence< RangeT > >
{
	inline std::size_t operator()( const IntegerSequence< RangeT >& s ) const { return s.hash(); }
};

template< typename RangeT >
struct equal_to< apegrunt::IntegerSequence< RangeT > >
{
	inline bool operator()( const IntegerSequence< RangeT >& lhs, const IntegerSequence< RangeT >& rhs ) const
	{
		return lhs == rhs;
	}
};


} // namespace apegrunt

#endif // APEGRUNT_INTEGERSEQUENCE_INTERFACE_HPP
