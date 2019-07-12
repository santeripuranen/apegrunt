/** @file Alignment_forward.h

	Copyright (c) 2016-2017 Santeri Puranen.

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

#ifndef APEGRUNT_ALIGNMENT_INTERFACE_FORWARD_H
#define APEGRUNT_ALIGNMENT_INTERFACE_FORWARD_H

#include <memory> // for std::shared_ptr

namespace apegrunt {

enum { StateBlock_size = 16 };

static constexpr std::size_t get_last_block_size( std::size_t n_columns ) { return n_columns == 0 ? 0 : ((n_columns-1) % StateBlock_size )+1; }
static constexpr std::size_t get_last_block_index( std::size_t n_columns ) { return n_columns == 0 ? 0 : (n_columns / StateBlock_size) - ( n_columns % StateBlock_size == 0 ? 1 : 0 ); }
static constexpr std::size_t get_number_of_blocks( std::size_t n_columns ) { return n_columns == 0 ? 0 : get_last_block_index(n_columns) + 1; }
static constexpr std::size_t get_block_index( std::size_t pos ) { return pos == 0 ? 0 : pos/StateBlock_size; }
static constexpr std::size_t get_pos_in_block( std::size_t pos ) { return pos%StateBlock_size; }

} // namespace apegrunt

namespace apegrunt {

// forward declarations
template< typename StateT >
class Alignment;

template< typename StateT >
using Alignment_ptr = std::shared_ptr< Alignment<StateT> >;
//using Alignment_ptr = std::shared_ptr< Alignment >;

template< typename AlignmentT, typename... Args >
Alignment_ptr< typename AlignmentT::state_t > make_Alignment_ptr( Args&&... args )
{
	return std::make_shared< AlignmentT >( std::forward<Args>(args)... );
}

template< typename AlignmentT >
Alignment_ptr< typename AlignmentT::state_t > make_Alignment_ptr( const AlignmentT& alignment )
{
	return std::make_shared< AlignmentT >( alignment );
}

} // namespace apegrunt

#endif // APEGRUNT_ALIGNMENT_INTERFACE_FORWARD_H

