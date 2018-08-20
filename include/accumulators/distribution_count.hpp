/** @file distribution_count.hpp

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
#ifndef APEGRUNT_ACCUMULATORS_DISTRIBUTION_COUNT_HPP
#define APEGRUNT_ACCUMULATORS_DISTRIBUTION_COUNT_HPP

#include <numeric> // for std::accumulate

#include <boost/mpl/placeholders.hpp>
#include <boost/accumulators/framework/accumulator_base.hpp>
#include <boost/accumulators/framework/extractor.hpp>
#include <boost/accumulators/framework/parameters/sample.hpp>
#include <boost/accumulators/framework/depends_on.hpp>

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/count.hpp>

#include "accumulators/distribution.hpp"

namespace boost {
namespace accumulators {

namespace impl {

template< typename Sample >
struct distribution_count_accumulator : accumulator_base
{
	using result_type = std::size_t;

	template< typename Args >
	distribution_count_accumulator( Args const & args ) {}

	template< typename Args >
	result_type result( Args const & args ) const
	{
		using std::cbegin; using std::cend;

		const auto& d = distribution( args[accumulator] );
		return std::accumulate( cbegin(d), cend(d), result_type(0), [=]( auto sum, const auto it ) { return sum += count(access_bin(it)); } );
	}
};

} // namespace impl

namespace tag {

// sample count calculated from a binned set of samples
struct distribution_count : depends_on< distribution >
{
	using impl = accumulators::impl::distribution_count_accumulator< mpl::_1 >;
};

} // namespace tag

namespace extract {

extractor< tag::distribution_count > const distribution_count = {};

} // namespace extract

using extract::distribution_count;

// For the purposes of feature-based dependency resolution,
// distribution_count provides the same feature as count
template<> struct feature_of<tag::distribution_count> : feature_of<tag::count> {};

// enable use of tag::count(distribution) as an alias of tag::sample_count
template<> struct as_feature<tag::count(from_distribution)> { using type = tag::distribution_count; };

} // namspace accumulators
} // namspace boost

namespace apegrunt {
namespace accumulators {

template< typename Iterable >
std::size_t distribution_count( Iterable range )
{
	using result_type = std::size_t;
	using boost::accumulators::access_bin;
	using boost::accumulators::count;

	return std::accumulate( cbegin(range), cend(range), result_type(0), [=]( auto sum, const auto it ) { return sum += count(access_bin(it)); } );
}

} // accumulators
} // apegrunt

#endif // APEGRUNT_ACCUMULATORS_DISTRIBUTION_COUNT_HPP
