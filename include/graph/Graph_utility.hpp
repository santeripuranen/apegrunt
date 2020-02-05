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
#include <algorithm>

#include "graph/Graph.hpp"
#include "apegrunt/Loci.h"
#include "apegrunt/Loci_parsers.hpp"

namespace apegrunt {

Loci_ptr extract_set_node_indices( const Graph_ptr& graph )
{
	std::vector<std::size_t> indices;
	for( const auto& edge: *graph )
	{
		if( bool(edge) )
		{
			indices.push_back( edge.node1() );
			indices.push_back( edge.node2() );
		}
		std::sort( indices.begin(), indices.end(), std::less<std::size_t>() );
	}

	return apegrunt::make_Loci_list(std::move(indices),0);
}

} // namespace apegrunt

#endif // APEGRUNT_GRAPH_UTILITY_HPP
