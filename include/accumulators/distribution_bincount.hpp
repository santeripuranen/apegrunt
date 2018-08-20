/** @file distribution_bincount.hpp

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
#ifndef APEGRUNT_ACCUMULATORS_DISTRIBUTION_BINCOUNT_HPP
#define APEGRUNT_ACCUMULATORS_DISTRIBUTION_BINCOUNT_HPP

#include <numeric> // for std::accumulate

#include <boost/mpl/placeholders.hpp>
#include <boost/accumulators/framework/accumulator_base.hpp>
#include <boost/accumulators/framework/extractor.hpp>
#include <boost/accumulators/framework/parameters/sample.hpp>
#include <boost/accumulators/framework/depends_on.hpp>

#include <boost/accumulators/accumulators.hpp>

#include "accumulators/distribution.hpp"

namespace boost {
namespace accumulators {

namespace impl {

template< typename Sample >
struct distribution_bincount_accumulator : accumulator_base
{
	using result_type = std::size_t;

	template< typename Args >
	distribution_bincount_accumulator( Args const & args ) {}

	template< typename Args >
	result_type result( Args const & args ) const
	{
		using std::cbegin; using std::cend;

		const auto& d = distribution( args[accumulator] );
		return std::accumulate( cbegin(d), cend(d), result_type(0), [=]( auto bincount ) { return ++bincount; } );
	}
};

} // namespace impl

namespace tag {

// sample count calculated from a binned set of samples
struct distribution_bincount : depends_on< distribution >
{
	using impl = accumulators::impl::distribution_bincount_accumulator< mpl::_1 >;
};

} // namespace tag

namespace extract {

extractor< tag::distribution_bincount > const distribution_bincount = {};

} // namespace extract

using extract::distribution_bincount;

} // namspace accumulators
} // namspace boost

#endif // APEGRUNT_ACCUMULATORS_DISTRIBUTION_BINCOUNT_HPP
