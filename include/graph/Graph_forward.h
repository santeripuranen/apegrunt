/** @file Graph_forward.h

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

#ifndef APEGRUNT_GRAPH_FORWARD_H
#define APEGRUNT_GRAPH_FORWARD_H

#include <memory> // std::shared_ptr, std::make_shared

namespace apegrunt {

// forward declarations
class Graph;

using Graph_ptr = std::shared_ptr< Graph >;

template< typename GraphT, typename... Args >
Graph_ptr make_Graph_ptr( Args&&... args )
{
	return std::make_shared< GraphT >( std::forward<Args>(args)... );
}

template< typename GraphT >
Graph_ptr make_Graph_ptr( const GraphT& graph )
{
	return std::make_shared< GraphT >( graph );
}

template< typename GraphT >
Graph_ptr make_Graph_ptr()
{
	return std::make_shared< GraphT >();
}

} // namespace apegrunt

#endif // APEGRUNT_LOCI_FORWARD_H

