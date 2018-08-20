/** @file distribution_generator_base.hpp

	Copyright (c) 2017 Santeri Puranen.

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
#ifndef APEGRUNT_ACCMULATORS_DISTRIBUTION_GENERATOR_BASE_HPP
#define APEGRUNT_ACCMULATORS_DISTRIBUTION_GENERATOR_BASE_HPP

#include "accumulators/distribution.hpp"

namespace apegrunt {

namespace accumulators {

template< typename Impl, typename Iterable, typename Bin >
struct distribution_generator_base
{
	using const_cast_t = const Impl* const;

	distribution_generator_base( Iterable&& arg_range, const Bin&(*arg_access_bin)(const typename Iterable::value_type&)=&boost::accumulators::access_bin<typename Iterable::value_type> ) :
		range(std::move(arg_range)),
		access_bin(arg_access_bin)
	{
	}

	inline void operator()( std::ostream& os ) const
	{
		return static_cast<const_cast_t>(this)->operator()(os); // forward to
	}

	Iterable&& range;
	const Bin&(*access_bin)(const typename Iterable::value_type&);
};

template< typename Impl, typename Iterator, typename Bin >
std::ostream& operator<< ( std::ostream& os, const distribution_generator_base<Impl, Iterator,Bin>& generator )
{
	generator( os );
	return os;
}

} // namespace accumulators

} // namespace apegrunt

#endif // APEGRUNT_ACCMULATORS_DISTRIBUTION_GENERATOR_BASE_HPP
