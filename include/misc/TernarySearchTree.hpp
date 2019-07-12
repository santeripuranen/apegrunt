/** @file TernarySearchTree.hpp

	Copyright (c) 2018-2019 Santeri Puranen.

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
#ifndef APEGRUNT_TERNARY_SEARCH_TREE_HPP
#define APEGRUNT_TERNARY_SEARCH_TREE_HPP

#include <cstdint>
#include <limits>
#include <memory>
//#include <vector>

namespace apegrunt {

template< typename KeyT, typename ValueT >
struct alignas(8) Node
{
	using my_type = Node<KeyT,ValueT>;

	using key_t = KeyT;
	using value_t = ValueT;

	using node_ptr_t = std::unique_ptr< my_type >; // let unique_ptr handle the clean-up

	enum { FLAG_MASK = value_t(1) << std::numeric_limits<value_t>::digits-1 };

	Node() = delete; // = default;
	~Node() = default;

	Node( my_type&& other )
	: m_equal( std::move(other.m_equal) ),
	  m_left( std::move(other.m_left) ),
	  m_right( std::move(other.m_right ) ),
	  m_key( std::move(other.m_key) ),
	  m_value( std::move(other.m_value ) )
	{
	}

	Node( key_t key ) : m_key(key), m_value(0) { }

	node_ptr_t m_equal;
	node_ptr_t m_left;
	node_ptr_t m_right;
	const key_t m_key;
	value_t m_value;

	inline bool operator==( const key_t& query ) const { return m_key == query; }
	inline bool operator!=( const key_t& query ) const { return m_key != query; }
	inline bool operator<( const key_t& query ) const { return m_key < query; }
	inline bool operator>( const key_t& query ) const { return m_key > query; }

	inline bool is_leaf() const { return m_value & FLAG_MASK; }
	inline value_t value() const { return m_value & ~FLAG_MASK; }
	inline void set_value( value_t value ) { m_value = value | FLAG_MASK; } // set value and leaf flag as only leaf nodes carry values

	inline void set_leaf( bool state=true ) { m_value = state ? (m_value | FLAG_MASK) : (m_value & ~FLAG_MASK); }
};

template< typename KeyT, typename ValueT >
class TernarySearchTree
{
public:
	using my_type = TernarySearchTree<KeyT,ValueT>;

	using key_fragment_t = typename KeyT::block64_t;
	using node_t = Node<key_fragment_t,ValueT>;
	using node_ptr_t = node_t*;
	using stored_node_t = typename node_t::node_ptr_t;

	TernarySearchTree() //= default;
	: m_root(), // only default initialized here
	  m_n_total_keys(0),
	  m_n_unique_keys(0),
	  m_n_nodes(0)
	{
	}

	~TernarySearchTree() = default;

	TernarySearchTree( my_type&& other )
	: m_root( std::move( other.m_root ) ),
	  m_n_total_keys( std::move( other.m_n_total_keys) ),
	  m_n_unique_keys( std::move( other.m_n_unique_keys) ),
	  m_n_nodes( std::move( other.m_n_nodes ) )
	{
	}

	TernarySearchTree( const KeyT& key, const ValueT& value )
	: m_root(), // only default initialized here
	  m_n_total_keys(0),
	  m_n_unique_keys(0),
	  m_n_nodes(0)
	{
		this->insert( key, value );
	}

	// returns valid node_ptr_t if key is found, else nullptr
	inline node_ptr_t find( KeyT key )
	{
		return nullptr; // dummy
	}

	// always returns a valid node_ptr_t corresponding to 'key', but inserts 'value' only if 'key' is new
	inline node_ptr_t insert( const KeyT& key, const ValueT& value )
	{
		if( !m_root ) { m_root = my_type::create_node( key.block64(0) ); ++m_n_nodes; }

		auto node = this->recursive_insert( m_root.get(), key, 0, KeyT::N64 );

		++m_n_total_keys;
		if( !(node->is_leaf()) ) { node->set_value( value ); ++m_n_unique_keys; }

		return node;
	}

	void statistics( std::ostream *out=nullptr ) const
	{
		if( out )
		{
			*out << "apegrunt::TernarySearchTree:"
					<< " #keys/unique=" << m_n_total_keys << "/" << m_n_unique_keys
					<< " #nodes=" << m_n_nodes
					<< " mem=" << apegrunt::memory_string( this->bytesize() )
					<< "\n"
			;
		}
	}

	inline std::size_t nkeys() const { return m_n_total_keys; }
	inline std::size_t nunique_keys() const { return m_n_unique_keys; }
	inline std::size_t nnodes() const { return m_n_nodes; }

	inline std::size_t bytesize() const { return m_n_nodes*sizeof(node_t)+3*sizeof(std::size_t); }

private:
	stored_node_t m_root;
	std::size_t m_n_total_keys;
	std::size_t m_n_unique_keys;
	std::size_t m_n_nodes;

	static stored_node_t create_node( key_fragment_t key ) { return std::make_unique<node_t>( key ); }

	inline node_ptr_t recursive_insert( node_ptr_t node, const KeyT& key, std::size_t i, const std::size_t N )
	{
		if( *node == key.block64(i) )
		{
			if( !(++i < N) ) { return node; }
			if( !node->m_equal )
			{
				//node->m_equal = node_t::create( key[i] );
				node->m_equal = my_type::create_node( key.block64(i) );
				++m_n_nodes;
			}
			node = node->m_equal.get();
		}
		else if( *node > key.block64(i) )
		{
			if( !node->m_left )
			{
				//node->m_left = node_t::create( key[i] );
				node->m_left = my_type::create_node( key.block64(i) );
				++m_n_nodes;
			}
			node = node->m_left.get();
		}
		else // if( *node < key.block32(i) )
		{
			if( !node->m_right )
			{
				//node->m_right = node_t::create( key[i] );
				node->m_right = my_type::create_node( key.block64(i) );
				++m_n_nodes;
			}
			node = node->m_right.get();
		}

		// ..and bake until done
		return this->recursive_insert( node, key, i, N );
	}
};

template< typename KeyT, typename ValueT >
std::size_t bytesize( const TernarySearchTree<KeyT,ValueT>& tree ) { return tree.bytesize(); }

} // namespace apegrunt

#endif // APEGRUNT_TERNARY_SEARCH_TREE_HPP
