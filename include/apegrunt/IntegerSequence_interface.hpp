/** @file IntegerSequence_interface.hpp
 
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

#ifndef APEGRUNT_INTEGERSEQUENCE_INTERFACE_HPP
#define APEGRUNT_INTEGERSEQUENCE_INTERFACE_HPP

#include <vector>
#include <numeric> // for std::accumulate

#include <boost/container_hash/hash.hpp>
#include <boost/dynamic_bitset.hpp>

#include "IntegerSequence_forward.h"
#include "IntegerSequence_operations_forward.h"
#include "IntegerSequence_iterator.hpp"

#include "Apegrunt_utility.hpp"

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
	using back_inserter_t = apegrunt::iterator::back_insert_iterator< container_t >;

	IntegerSequence() : m_back_inserter( m_storage ), m_size(0), m_hash(0), m_block_mask(0) { }
	IntegerSequence( value_type value ) : m_back_inserter( m_storage ), m_size(0), m_hash(0), m_block_mask(0) { this->push_back(value); }
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

	IntegerSequence( const std::vector<value_type>& dense )
	: m_back_inserter( m_storage ), m_size(0), m_hash(0), m_block_mask(0) { for( auto pos: dense ) { this->push_back(pos); } }

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

	inline void push_back( value_type value ) { m_back_inserter = value; ++m_size; this->update_block_mask( my_type::generate_block_mask(value) ); }
	inline void push_back( const range_t& element ) { m_back_inserter = element; m_size+=element.size(); this->update_block_mask( my_type::generate_block_mask(element) ); }
	inline void emplace_back( range_t&& element ) { m_size+=element.size(); this->update_block_mask( my_type::generate_block_mask(element) ); m_back_inserter = std::move(element); }

	// merge() invalidates iterators, pointers, references to elements of 'this'; 'other' remains intact
	inline my_type& merge( const my_type& other ) { *this = std::move( apegrunt::set_union( *this, other ) ); return *this; }

	// Quick check for potential overlap. Bitsets are guaranteed _not_ to overlap if 'false' and may overlap if 'true' (but don't necessarily do so)
	inline bool has_overlap( const my_type& other ) const { return m_block_mask & other.m_block_mask; }

	// An IntegerSequence is empty if its internal storage contains nothing.
	// More consistent performance (O(1)) than "size() != 0", as is_empty() has no side effects.
	inline bool empty() const { return m_storage.empty(); }

	// check if the container contains the specified value. O(log N) complexity.
	inline bool contains( const value_type& value ) const
	{
		if( !( m_block_mask & my_type::generate_block_mask(value) ) ) { return false; }

		const auto itr = apegrunt::binary_search( m_storage.begin(), m_storage.end(), value, [](const auto& lhs, const auto& rhs) {	return lhs < rhs; } );
		return itr != m_storage.end() && itr->holds(value);
	}

	// Return number of stored values. Best case performance is O(1), worst case O(N), where N is the size of the internal container.
	inline std::size_t size() const { /*if( !m_storage.empty() && !m_size ) */ { this->update_size(); } return m_size; }
	inline std::size_t storagesize() const { return m_storage.size(); }
	inline std::size_t bytesize() const { return sizeof(range_t) * m_storage.size() /*+ sizeof(container_t)+sizeof(back_inserter_t)+sizeof(std::size_t)*2+sizeof(block_mask_type)*/; }

	inline std::size_t hash() const { return m_hash == 0 ? this->cache_the_hash() : m_hash; } // return cached hash, creating the hash value as needed

	const_iterator cbegin() const { return const_iterator( m_storage.data(), this->one_past_end() ); }
	const_iterator cend() const { return const_iterator( this->one_past_end(), this->one_past_end() ); }

	const_iterator begin() const { return this->cbegin(); }
	const_iterator end() const { return this->cend(); }

	const container_t& storage() const { return m_storage; }
	container_t& storage() { return m_storage; }

	inline void append_block_mask( const range_t& element ) { this->update_block_mask( my_type::generate_block_mask( element ) ); }
	inline void clear_block_mask() { m_block_mask = 0; }

	void shrink_to_fit() { m_storage.shrink_to_fit(); }

	// for testing purposes only
	inline std::size_t compressed_bytesize() const
	{
		static const auto N=std::numeric_limits<uint8_t>::digits;

		std::size_t bs(0);

		for( const auto& element: m_storage )
		{
			++bs; // first byte in compressed element would never be null
			for( std::size_t i=1; i<range_t::ELEMENT_SIZE/N; ++i )
			{
				bs += (0 == element.template block<N>(i));
			}
		}

		return bs;
	}

private:
	using block_mask_type = uint64_t;

	container_t m_storage;
	back_inserter_t m_back_inserter;
	mutable std::size_t m_size;
	mutable std::size_t m_hash;
	block_mask_type m_block_mask;

	inline std::size_t cache_the_hash() const
	{
		m_hash = boost::hash_range( this->cbegin(), this->cend() );
		return m_hash;
	}

	// single element block mask
	static constexpr const block_mask_type generate_block_mask( value_type value ) noexcept { return block_mask_type(1) << (value >> range_t::BLOCK_NSTEPS); }

	// range element block mask
	static constexpr const block_mask_type generate_block_mask( value_type begin, value_type end ) noexcept
	{
		block_mask_type mask(0);
		for( auto shift(begin >> range_t::BLOCK_NSTEPS); shift < (end >> range_t::BLOCK_NSTEPS)+1; ++shift )
		{
			mask |= (block_mask_type(1) << shift);
		}
		return mask;
	}

	static constexpr const block_mask_type generate_block_mask( const range_t& element )
	{
		return element.is_range() ? my_type::generate_block_mask( element.range_begin(), element.range_end() ) : my_type::generate_block_mask( element() );
	}

	inline void update_size() const
	{
		m_size = std::accumulate( std::cbegin(m_storage), std::cend(m_storage), 0, [](const auto& sum, const auto& a) { return sum+a.size(); } );
	}

	inline void update_block_mask( block_mask_type mask ) { m_block_mask |= mask; }

	// For internal use only. The returned pointer should never be dereferenced.
	inline const range_t * const one_past_end() const { return &(m_storage.back())+1; }

	iterator begin() { return iterator( m_storage.data(), this->one_past_end() ); }
	iterator end() { return iterator( this->one_past_end(), this->one_past_end() ); }

	// friend declarations for intrusive set operations -- should be used sparingly
	// These functions either modify set contents in-place or the implementation simply benefits from access to private members
	template< typename UnaryFunction > friend UnaryFunction inplace_set_differences_and_for_each_in_intersection( my_type& a, my_type& b, UnaryFunction f );
	friend my_type set_union<value_type>( const my_type& a, const my_type& b );
	friend my_type set_union<value_type>( const std::vector< std::reference_wrapper< const my_type > >& sets );
};

template< typename RangeT >
bool operator==( const IntegerSequence< RangeT >& a, const IntegerSequence< RangeT >& b ) { return a.hash() == b.hash(); }

template< typename RangeT >
struct hash< apegrunt::IntegerSequence< RangeT > >
{
	std::size_t operator()( const IntegerSequence< RangeT >& s ) const { return s.hash(); }
};

} // namespace apegrunt

#endif // APEGRUNT_INTEGERSEQUENCE_INTERFACE_HPP
