/** @file Graph.hpp

	Copyright (c) 2018 Santeri Puranen.

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
#ifndef APEGRUNT_GRAPH_HPP
#define APEGRUNT_GRAPH_HPP

#include <memory> // for std::shared_ptr and std::make_shared
#include <vector>
#include <set>
#include <cmath>
#include <mutex> // for std::mutex
#include <limits>

#include "Graph_forward.h"
#include "apegrunt/Apegrunt_utility.hpp"

namespace apegrunt {

template< typename IntegerT >
constexpr IntegerT create_flag_mask() { return IntegerT(1) << (std::numeric_limits<IntegerT>::digits-1); }

class EdgeID : public apegrunt::extend_comparison_operators
{
public:
	using node_id_t = uint32_t;
	using id_t = uint64_t;
	using my_type = EdgeID;

	enum { ID_FLAG_MASK = create_flag_mask<id_t>() };
	enum { ID_UNFLAG_MASK = ~create_flag_mask<id_t>() };
	enum { NODE_ID_UNFLAG_MASK = ~create_flag_mask<node_id_t>() };

	EdgeID() = delete;
	inline EdgeID( node_id_t i, node_id_t j, bool set=false )
	: m_i( std::min(i,j) & NODE_ID_UNFLAG_MASK ),
	  m_j( std::max(i,j) & NODE_ID_UNFLAG_MASK )
	{
		if(set) { this->set(); }
	}
	inline EdgeID( id_t id ) { m_id.store(id); }
	~EdgeID() = default;

	inline EdgeID( const my_type& other ) { m_id.store( other.m_id.load() ); }
	inline EdgeID( my_type&& other ) { m_id.store( other.m_id.load() ); }

	inline my_type& operator=( const my_type& other ) { m_id.store(other.m_id.load()); return *this; }

	inline node_id_t first() const { return m_i & NODE_ID_UNFLAG_MASK; }
	inline node_id_t second() const { return m_j & NODE_ID_UNFLAG_MASK; }
	inline id_t id() const { return m_id.load() & ID_UNFLAG_MASK; }

	inline operator bool() const { return m_id.load() & ID_FLAG_MASK; }
	inline bool is_set() const { return bool(*this); }
	inline void set() { m_id.fetch_or(ID_FLAG_MASK); }
	inline void unset() { m_id.fetch_and(ID_UNFLAG_MASK); }

	inline bool operator==( const my_type& rhs ) const { return this->id() == rhs.id(); }
	inline bool operator<( const my_type& rhs ) const { return this->id() < rhs.id(); }
	inline bool operator!=( const my_type& rhs ) const { return this->id() != rhs.id(); }

private:
	using node_id_intl_t = node_id_t;
	using id_intl_t = std::atomic<id_t>;

	union
	{
		struct {
			node_id_intl_t m_i;
			node_id_intl_t m_j; // we use sign as flag bit
		};
		id_intl_t m_id;
	};
};

static std::ostream& operator<< ( std::ostream& os, const EdgeID& id )
{
	os << id.first() << " " << id.second();
	return os;
}

// an undirected edge class
class Edge : public apegrunt::extend_comparison_operators
{
public:
	using id_t = EdgeID;
	using node_id_t = typename id_t::node_id_t;
	using real_t = double;
	using my_type = Edge;

	Edge() = delete;
	inline Edge( node_id_t i, node_id_t j, real_t w ) : m_id(i,j), m_w(w) { }
	~Edge() = default;

	inline Edge( const my_type& other ) : m_id( other.m_id ), m_w( other.m_w ) { }
	inline Edge( my_type&& other ) : m_id( std::move(other.m_id) ), m_w( std::move(other.m_w) ) { }

	inline my_type& operator=( const my_type& rhs ) { m_id = rhs.m_id; m_w = rhs.m_w; return *this; }

	inline id_t id() const { return m_id; }
	inline typename id_t::node_id_t node1() const { return m_id.first(); }
	inline typename id_t::node_id_t node2() const { return m_id.second(); }
	inline real_t weight() const { return m_w; }
	inline operator bool() const { return bool(m_id); }
	inline void set() { m_id.set(); }
	inline void unset() { m_id.unset(); }

	inline bool operator==( const my_type& rhs ) const { return m_id == rhs.m_id; }
	inline bool operator<( const my_type& rhs ) const { return m_id < rhs.m_id; }

private:
	id_t m_id;
	real_t m_w;
};

static std::ostream& operator<< ( std::ostream& os, const Edge& edge )
{
	os << edge.id() << " " << edge.weight();
	return os;
}

class Graph
{
public:
	using my_type = Graph;

	using edge_t = Edge;
	using node_id_t = typename edge_t::node_id_t;
	using real_t = typename edge_t::real_t;

	using node_storage_t = std::set<node_id_t>;
	using node_itr_t = typename std::set<node_id_t>::iterator;
	using node_const_itr_t = typename std::set<node_id_t>::const_iterator;

	using edge_storage_t = std::vector< edge_t >;
	using edge_itr_t = typename edge_storage_t::iterator;
	using edge_const_itr_t = typename edge_storage_t::const_iterator;

	using lock_t = std::unique_lock<std::mutex>;

	Graph() = default;
	~Graph() = default;

	Graph( my_type&& other )
	: //m_nodes( std::move(other.m_nodes) ),
	  m_edges( std::move(other.m_edges) )
	{
	}

	inline edge_t& operator[]( std::size_t i ) { return m_edges[i]; }
	inline const edge_t& operator[]( std::size_t i ) const { return m_edges[i]; }

	// Add a now node. Caller should lock before calling add(), if used in a parallel setting.
	inline void add( node_id_t i, node_id_t j, real_t weight )
	{
		m_edges.emplace_back( i, j, weight );
	}

	inline void join( my_type& other )
	{
		//auto&& scope_lock_this = this->acquire_lock();
		//auto&& scope_lock_other = other.acquire_lock();
		this->lock(); other.lock();
		//std::cout << "apegrunt: Graph.join() - this->size()=" << this->size() << " other.size()=" << other.size() << std::endl;
		this->reserve( this->size()+other.size() );
		for( auto& edge: other ) { m_edges.push_back( edge ); }
		other.clear();
		other.unlock(); this->unlock();
	}

	inline lock_t acquire_lock() { lock_t lock(m_modify); lock.lock(); return std::move( lock ); }
	inline void lock() { m_modify.lock(); }
	inline void unlock() { m_modify.unlock(); }
	inline void clear() { m_edges.clear(); }

	inline node_id_t capacity() const { return m_edges.capacity(); }
	inline void reserve( std::size_t n ) { m_edges.reserve(n); }
	inline std::size_t size() const { return m_edges.size(); }

	inline edge_itr_t begin() { using std::begin; return begin(m_edges); }
	inline edge_itr_t end() { using std::end; return end(m_edges); }

	inline edge_const_itr_t begin() const { using std::cbegin; return cbegin(m_edges); }
	inline edge_const_itr_t end() const { using std::cend; return cend(m_edges); }

	// Convenience functions

	// Find edge in graph using binary search; caller should make sure that edges are sorted in ascending order of id() first,
	// i.e. this->sort( []( const auto& a, const auto& b ) { return a.id() < b.id(); } )
	inline edge_itr_t find( const edge_t& edge ) { return apegrunt::binary_search( this->begin(), this->end(), edge ); }
	inline edge_const_itr_t find( const edge_t& edge ) const { return apegrunt::binary_search( this->begin(), this->end(), edge ); }

	// sort edges in descending order of weight, or whatever you'd prefer
	template< typename Compare=bool(const edge_t&, const edge_t&) >
	inline void sort( Compare comp=[](const edge_t& a, const edge_t& b) { return (a.weight() == b.weight()) ? (a.id() < b.id()) : (a.weight() > b.weight()); } )
	{
		std::sort( this->begin(), this->end(), comp );
	}

private:
	//node_storage_t m_nodes;
	edge_storage_t m_edges;
	std::mutex m_modify;
};

static std::ostream& operator<< ( std::ostream& os, const Graph& graph )
{
	for( const auto& edge: graph )
	{
		os << edge << " " << bool(edge) << "\n";
	}
	return os;
}

static std::ostream& operator<< ( std::ostream& os, const Graph_ptr& graph_ptr )
{
	os << *(graph_ptr.get());
	return os;
}

} // namespace apegrunt

#endif // APEGRUNT_GRAPH_HPP
