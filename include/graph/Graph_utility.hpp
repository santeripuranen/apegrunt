/** @file Graph_utility.hpp

	Copyright (c) 2019-2020 Santeri Puranen.

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
#ifndef APEGRUNT_GRAPH_UTILITY_HPP
#define APEGRUNT_GRAPH_UTILITY_HPP

#include <vector>
#include <unordered_set>
#include <algorithm>

#include "graph/Graph.hpp"
#include "apegrunt/Loci.h"
#include "apegrunt/Loci_parsers.hpp"

namespace apegrunt {

// UnaryPredicate should accept a const Edge& as argument and return
// a bool indicating whether to extract the node indices of that edge
template< typename UnaryPredicate >
Loci_ptr extract_node_indices( const Graph_ptr& graph, UnaryPredicate p )
{
	std::vector<std::size_t> indices;
	for( const auto& edge: *graph )
	{
		if( p(edge) )
		{
			indices.push_back( edge.node1() );
			indices.push_back( edge.node2() );
		}
	}
	std::sort( indices.begin(), indices.end(), std::less<std::size_t>() );
	auto end = std::unique( indices.begin(), indices.end() );
	indices.erase( end, indices.end() );
	indices.shrink_to_fit();

	return apegrunt::make_Loci_list(std::move(indices),0);
}

} // namespace apegrunt

#endif // APEGRUNT_GRAPH_UTILITY_HPP
