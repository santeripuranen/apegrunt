/** @file Loci_generators.hpp

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

#ifndef APEGRUNT_LOCI_GENERATORS_HPP
#define APEGRUNT_LOCI_GENERATORS_HPP

#include <ostream>

#include "Loci_interface.hpp"

namespace apegrunt {

template< typename LociListT >
bool generate_Loci_list( Loci_ptr loci, std::ostream *outstream )
{
// this is temporary file output code, used until we have a proper generator implementation.
	if( outstream && outstream->good() )
	{
		const auto base_index = Apegrunt_options::get_output_indexing_base();
		for( const auto locus: loci )
		{
			*outstream << locus+base_index << "\n";
		}
		return true;
	}
	return false;
}

} // namespace apegrunt

#endif // APEGRUNT_LOCI_GENERATORS_HPP

