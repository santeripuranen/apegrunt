/** @file Graph_output_format.hpp

	Copyright (c) 2019 Santeri Puranen.

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
#ifndef APEGRUNT_GRAPH_OUTPUT_FORMAT_HPP
#define APEGRUNT_GRAPH_OUTPUT_FORMAT_HPP

#include "graph/Graph.hpp"
#include "apegrunt/Alignment.h"

namespace apegrunt {

template< typename StateT >
struct Graph_output_formatter
{
	Graph_output_formatter( const Graph_ptr& graph, const Alignment_ptr<StateT>& alignment ) //, std::size_t (*	const distance)(std::size_t,std::size_t) )
	: m_graph(graph),
	  m_alignment(alignment)
	  //m_distance(distance)
	{
	}

	const Graph_ptr m_graph;
	const Alignment_ptr<StateT> m_alignment;
	//std::size_t (*const m_distance)(std::size_t,std::size_t);
};

template< typename StateT >
static std::ostream& operator<< ( std::ostream& os, const Graph_output_formatter<StateT>& gof )
{
	const auto& index_translation = *(gof.m_alignment->get_loci_translation());

	const std::size_t base_index = apegrunt::Apegrunt_options::get_output_indexing_base();

	for( const auto& edge: *(gof.m_graph) )
	{
		const auto index1 = index_translation[edge.id().first()]+base_index;
		const auto index2 = index_translation[edge.id().second()]+base_index;

		os << index1 << " "
				<< index2 << " "
				//<< m_distance(index1,index2) << " "
				//<< edge << " "
				<< edge.weight()
				<< bool(edge)
				<< "\n";
	}
	return os;
}

} // namespace apegrunt

#endif // APEGRUNT_GRAPH_OUTPUT_FORMAT_HPP
