/** @file Loci_operations.hpp
 
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
#ifndef APEGRUNT_LOCI_OPERATIONS_HPP
#define APEGRUNT_LOCI_OPERATIONS_HPP

#include <string>
#include <algorithm> // for std::find

#include <boost/serialization/serialization.hpp>
#include <boost/serialization/nvp.hpp>

#include "Loci_forward.h"
#include "Loci_iterator.hpp"
#include "Loci_impl_default_storage.hpp"

namespace apegrunt {

template<typename LociT=Loci_impl_default_storage<std::size_t> >
Loci_ptr operator-( Loci_ptr lhs, Loci_ptr rhs )
{
	std::vector<std::size_t> difference; difference.reserve( lhs->size() - rhs->size() ); // speculative
	for( auto index: lhs )
	{
		if( rhs->cend() == std::find( rhs->cbegin(), rhs->cend(), index ) ) { difference.push_back( index ); }
	}
	return make_Loci_ptr<LociT>( std::move( difference ) );
}

std::ostream& operator<< ( std::ostream& os, const Loci* locilist )
{
	for( auto i: *locilist ) { os << " " << i; }
	return os;
}

std::ostream& operator<< ( std::ostream& os, const Loci_ptr& locilist )
{
	os << locilist.get();
	return os;
}


} // namespace apegrunt

#endif // APEGRUNT_LOCI_INTERFACE_HPP
