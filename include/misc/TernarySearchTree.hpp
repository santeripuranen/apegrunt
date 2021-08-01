/** @file TernarySearchTree.hpp

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
#ifndef APEGRUNT_TERNARY_SEARCH_TREE_HPP
#define APEGRUNT_TERNARY_SEARCH_TREE_HPP

#include <cstdint>
#include <limits>
#include <memory>

namespace apegrunt {

template< typename KeyT, typename ValueT >
struct key_value
{
	key_value( const KeyT& k ) : key(k), value() { }
	key_value( const KeyT& k, const ValueT& v ) : key(k), value(v) { }

	KeyT key;
	ValueT value;
};

template< typename KeyFragmentT, typename KeyT, typename ValueT >
struct alignas(8) Node
{
	using my_type = Node<KeyFragmentT,KeyT,ValueT>;

	using key_t = KeyFragmentT;
	using value_t = ValueT;
	using key_value_t = key_value<KeyT,value_t>;


	using node_ptr_t = std::unique_ptr< my_type >; // let unique_ptr handle the clean-up
	using key_value_ptr_t = std::unique_ptr< key_value_t >; // let unique_ptr handle the clean-up

	Node() = delete; // = default;
	~Node() = default;

	Node( key_t key ) : m_key(key) { } //, m_equal(nullptr), m_left(nullptr), m_right(nullptr) { }

	Node( my_type&& other )
	: m_key( std::move(other.m_key) ),
	  m_equal( std::move(other.m_equal) ),
	  m_left( std::move(other.m_left) ),
	  m_right( std::move(other.m_right ) ),
	  //m_equal( other.m_equal ),
	  //m_left( other.m_left ),
	  //m_right( other.m_right ),
	  m_key_value( std::move(other.m_key_value ) )
	{
	}

	const key_t m_key;
	node_ptr_t m_equal;
	node_ptr_t m_left;
	node_ptr_t m_right;
	//my_type* m_equal;
	//my_type* m_left;
	//my_type* m_right;
	key_value_ptr_t m_key_value;

	inline bool operator==( const key_t& query ) const { return m_key == query; }
	inline bool operator!=( const key_t& query ) const { return m_key != query; }
	inline bool operator<( const key_t& query ) const { return m_key < query; }
	inline bool operator>( const key_t& query ) const { return m_key > query; }

	inline bool is_leaf() const { return bool(m_key_value); } // only leaf nodes carry values
	inline key_value_t* value() { return m_key_value.get(); }
	inline void set_value( key_value_ptr_t key_value=std::make_unique<key_value_t>() ) noexcept { m_key_value = std::move(key_value); } // set value
};

template< typename KeyT, typename ValueT, std::size_t BlockSize=1 >
class TernarySearchTree
{
public:
	using my_type = TernarySearchTree<KeyT,ValueT,BlockSize>;

	using key_fragment_t = typename KeyT::template sub_block_t<BlockSize>;
	using node_t = Node<key_fragment_t,KeyT,ValueT>;
	using node_ptr_t = node_t*;
	using stored_node_t = typename node_t::node_ptr_t;

	using key_t = KeyT;
	using value_t = ValueT;
	using key_value_t = key_value<key_t,value_t>;

	enum { N=KeyT::N/BlockSize };

	TernarySearchTree() //= default;
	: m_insert_pimpl(&my_type::insert_impl_init),
	  m_latest_entry(),
	  m_root(),//m_root( my_type::create_node( key_fragment_t() ) ), // only default initialized here
	  m_n_total_keys(0),
	  m_n_unique_keys(0),
	  m_n_nodes(0)//, m_buffer( new unsigned char[768*48] )
	{
	}

	~TernarySearchTree() = default;

	TernarySearchTree( my_type&& other )
	: m_insert_pimpl(&my_type::insert_impl_init), //m_insert_pimpl( other.m_insert_pimpl ),
	  m_latest_entry( std::move( other.m_latest_entry) ),
	  m_root( std::move( other.m_root ) ),
	  m_n_total_keys( std::move( other.m_n_total_keys) ),
	  m_n_unique_keys( std::move( other.m_n_unique_keys) ),
	  m_n_nodes( std::move( other.m_n_nodes ) )//, m_buffer( std::move( other.m_buffer) )
	{
	}

	TernarySearchTree( const key_t& key )
	: m_insert_pimpl(&my_type::insert_impl_init),
	  m_latest_entry(),
	  m_root(),//m_root( my_type::create_node( key.template block<BlockSize>(0) ) ),
	  m_n_total_keys(0),
	  m_n_unique_keys(0),
	  m_n_nodes(0)//, m_buffer( new unsigned char[768*48] )
	{
		this->insert( key );
	}

	// Returns a valid key_value_t* if key is found, else nullptr
	inline key_value_t* find( const key_t& key )
	{
		return nullptr; // dummy
	}

	// Always returns a valid key_value_t*, either to a previously stored key_value_t or a freshly created
	// Not thread-safe
	inline key_value_t* insert( const key_t& key )
	{
		return (this->*m_insert_pimpl)(key);
	}

	template< typename Ftor >
	void for_each( Ftor& f ) //const
	{
		this->for_each_impl( m_root.get(), f );
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
	friend inline void swap( my_type& lhs, my_type& rhs )
	{
		using std::swap;
		swap( lhs.m_root, rhs.m_root );
		swap( lhs.m_latest_entry, rhs.m_latest_entry );
		swap( lhs.m_n_total_keys, rhs.m_n_total_keys );
		swap( lhs.m_n_unique_keys, rhs.m_n_unique_keys );
		swap( lhs.m_n_nodes, rhs.m_n_nodes );
		//swap( lhs.m_buffer, rhs.m_buffer );
	}

private:
	key_value_t* (my_type::*m_insert_pimpl)(const key_t& key);
	key_value_t* m_latest_entry;
	stored_node_t m_root;
	std::size_t m_n_total_keys;
	std::size_t m_n_unique_keys;
	std::size_t m_n_nodes;

	// init function, used only on the first insert() call; initialize root node and last-entry cache
	inline key_value_t* insert_impl_init(const key_t& key)
	{
		my_type::m_insert_pimpl = &my_type::insert_impl_cached;
		m_root = my_type::create_node( key.template blockview<BlockSize>(0) );
		++m_n_total_keys;
		++m_n_unique_keys;
		++m_n_nodes;
		return m_latest_entry = this->recursive_insert( m_root.get(), key, 0 );
	}

	inline key_value_t* insert_impl_cached(const key_t& key)
	{
		++m_n_total_keys;
		return key == m_latest_entry->key ? m_latest_entry : m_latest_entry = this->recursive_insert( m_root.get(), key, 0 );
	}

	//alignas(sizeof(stored_node_t)) unsigned char m_buffer[sizeof(stored_node_t)*512];
	//std::unique_ptr<unsigned char []> m_buffer;
	//alignas(8) unsigned char m_buffer[40*1024];

	template< typename Key >
	static inline stored_node_t create_node( const Key& key ) { return std::make_unique<node_t>( key ); }
	//static inline stored_node_t create_node( const Key& key ) { auto node = std::make_unique<node_t>( key ); std::cout << ((std::size_t)node.get() & 0x3F) << " "; return node; }

	//static inline stored_node_t create_node( key_fragment_t key ) { return std::make_unique<node_t>( key ); }

	template< typename Key >
	static inline node_ptr_t create_node( const Key& key, unsigned char *where ) { return new (where) node_t(key); }

	inline key_value_t* recursive_insert( node_ptr_t node, const key_t& key, std::size_t i )
	{
		const auto c = key.template blockview<BlockSize>(i).cmp(node->m_key);

		if( c < 0 )
		{
			if( !node->m_left && ++m_n_nodes )
			{
				node->m_left = my_type::create_node( key.template blockview<BlockSize>(i) );
			}
			node = node->m_left.get();
		}
		else
		{
			if( c > 0 )
			{
				if( !node->m_right && ++m_n_nodes )
				{
					node->m_right = my_type::create_node( key.template blockview<BlockSize>(i) );
				}
				node = node->m_right.get();
			}
			else
			{
				if( ++i < N )
				{
					if( !node->m_equal && ++m_n_nodes )
					{
						node->m_equal = my_type::create_node( key.template blockview<BlockSize>(i) );
					}
					node = node->m_equal.get();
				}
				else
				{
					if( !node->is_leaf() &&  ++m_n_unique_keys ) { node->set_value( std::make_unique<key_value_t>( key ) ); }
					return node->value(); // get out
				}
			}
		}

		// ..and bake until done
		return this->recursive_insert( node, key, i );
	}

	template< typename Ftor >
	inline void for_each_impl( node_ptr_t const node, Ftor& f ) const
	{
		if( node )
		{
			if( node->is_leaf() ) { f(node->value()); }
			this->for_each_impl( node->m_left.get(), f );
			this->for_each_impl( node->m_equal.get(), f );
			this->for_each_impl( node->m_right.get(), f );
		}
	}
};

template< typename KeyT, typename ValueT >
std::size_t bytesize( const TernarySearchTree<KeyT,ValueT>& tree ) { return tree.bytesize(); }

} // namespace apegrunt

#endif // APEGRUNT_TERNARY_SEARCH_TREE_HPP
