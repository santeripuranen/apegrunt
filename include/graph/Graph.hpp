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

#include "Graph_forward.h"
#include "apegrunt/Apegrunt_utility.hpp"

namespace apegrunt {

struct EdgeID : public apegrunt::extend_comparison_operators
{
	using node_id_t = int32_t;
	using id_t = int64_t;
	using my_type = EdgeID;

	enum { FLAG_MASK=-id_t(0) };

	//EdgeID() : m_id(0) { }

	EdgeID( node_id_t i, node_id_t j ) : m_i( std::min(i,j) ), m_j( std::max(i,j) )
	{
		//std::cout << "EdgeID(): m_i=" << m_i << " m_j=" << m_j << " m_id=" << m_id << std::endl;
	}
	EdgeID( id_t id ) : m_id(id) { }
	~EdgeID() = default;

	EdgeID( const my_type& other ) : m_id( other.m_id ) { }
	EdgeID( my_type&& other ) : m_id( std::move(other.m_id) ) { }

	constexpr my_type& operator=( const my_type& other ) { m_id = other.m_id; }

	union
	{
		struct {
			node_id_t m_i; // we use sign as flag bit
			node_id_t m_j;
		};
		id_t m_id;
	};

	inline node_id_t first() const { return std::abs(m_i); }
	inline node_id_t second() const { return std::abs(m_j); }
	inline id_t id() const { return std::abs(m_id); }

	inline constexpr operator bool() const { return std::signbit(m_id); }
	inline constexpr void set() { m_id = m_id | FLAG_MASK; }
	inline constexpr void unset() { m_id = m_id ^ FLAG_MASK; }

	inline bool operator==( const my_type& rhs ) const { return this->id() == rhs.id(); }
	inline bool operator<( const my_type& rhs ) const { return this->id() < rhs.id(); }
};

static std::ostream& operator<< ( std::ostream& os, const EdgeID& id )
{
	os << id.first() << " " << id.second();
	return os;
}

// an undirected edge class
struct Edge : public apegrunt::extend_comparison_operators
{
	using id_t = EdgeID;
	using node_id_t = typename id_t::node_id_t;
	using real_t = double;
	using my_type = Edge;

	Edge() = delete;
	Edge( node_id_t i, node_id_t j, real_t w ) : m_id(i,j), m_w(w) { }
	~Edge() = default;

	Edge( const my_type& other ) : m_id( other.m_id ), m_w( other.m_w ) { }
	Edge( my_type&& other ) : m_id( std::move(other.m_id) ), m_w( std::move(other.m_w) ) { }

	constexpr my_type& operator=( const my_type& rhs ) { m_id = rhs.m_id; m_w = rhs.m_w; }

	id_t m_id;
	real_t m_w;

	inline id_t id() const { return m_id; }
	inline typename id_t::node_id_t node1() const { return m_id.first(); }
	inline typename id_t::node_id_t node2() const { return m_id.second(); }
	inline constexpr real_t weight() const { return m_w; }
	inline constexpr operator bool() const { return bool(m_id); }
	inline constexpr void set() { m_id.set(); }
	inline constexpr void unset() { m_id.unset(); }

	inline bool operator==( const my_type& rhs ) const { return m_id == rhs.m_id; }
	inline bool operator<( const my_type& rhs ) const { return m_id < rhs.m_id; }
};

static std::ostream& operator<< ( std::ostream& os, const Edge& edge )
{
	os << edge.m_id << " " << edge.m_w;
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
	: m_nodes( std::move(other.m_nodes) ),
	  m_edges( std::move(other.m_edges) )
	{
	}

	inline void add( node_id_t i, node_id_t j, real_t weight )
	{
		//std::cout << "apegrunt::Graph.add(" << i << "," << j << "," << weight << ")" << std::endl;
		//m_modify.lock();
		//std::cout << "apegrunt::Graph: insert node" << std::endl;
		//m_nodes.insert(i); m_nodes.insert(j);
		//std::cout << "apegrunt::Graph: insert edge" << std::endl;
		m_edges.emplace_back( i, j, weight );
		//m_modify.unlock();
	}

	void join( my_type& other )
	{
		auto&& scope_lock = this->lock();
		other.lock();
		this->reserve( this->size()+other.size() );
		for( auto& edge: other ) { m_edges.push_back( edge ); }
	}

	inline lock_t lock() { return std::move( lock_t(m_modify) ); }

	inline node_id_t capacity() const { return m_edges.capacity(); }
	inline void reserve( std::size_t n ) { m_edges.reserve(n); }
	inline std::size_t size() const { return m_edges.size(); }

	// sort in descending order of weight
	inline void sort() { std::sort( this->begin(), this->end(), [](edge_t a, edge_t b) { return a.m_w > b.m_w; } ); }

	inline edge_itr_t begin() { using std::begin; return begin(m_edges); }
	inline edge_itr_t end() { using std::end; return end(m_edges); }

	inline edge_const_itr_t begin() const { using std::cbegin; return cbegin(m_edges); }
	inline edge_const_itr_t end() const { using std::cend; return cend(m_edges); }

private:
	node_storage_t m_nodes;
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
