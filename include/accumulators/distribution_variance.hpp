/** @file distribution_variance.hpp

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
#ifndef APEGRUNT_ACCUMULATORS_DISTRIBUTION_VARIANCE_HPP
#define APEGRUNT_ACCUMULATORS_DISTRIBUTION_VARIANCE_HPP

#include <cmath>
#include <numeric> // for std::accumulate

#include <boost/mpl/placeholders.hpp>
#include <boost/accumulators/framework/accumulator_base.hpp>
#include <boost/accumulators/framework/extractor.hpp>
#include <boost/accumulators/framework/depends_on.hpp>

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/count.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/variance.hpp>

#include "accumulators/distribution.hpp"
#include "accumulators/distribution_count.hpp"
#include "accumulators/distribution_mean.hpp"

namespace boost {
namespace accumulators {

namespace impl {

template< typename Sample >
struct distribution_variance_accumulator : accumulator_base
{
	using result_type = Sample;

	template< typename Args >
	distribution_variance_accumulator( Args const & args ) {}

	template< typename Args >
	result_type result( Args const & args ) const
	{
		using std::pow;
		using std::max;
		using std::cbegin; using std::cend;

		const auto& d = distribution( args[accumulator] );
		const result_type varsum = std::accumulate( cbegin(d), cend(d), result_type(0), [=]( auto sum, const auto it ) { return sum += ( variance(access_bin(it)) + pow(mean(access_bin(it)),2) ) * count(access_bin(it)); } );

		return ( varsum / result_type(max(std::size_t(1),distribution_count(args[accumulator])-1)) ) - pow(distribution_mean(args[accumulator]),2);
	}
};

} // namespace impl

namespace tag {

// distribution variance calculated from a binned set of samples
struct distribution_variance : depends_on< distribution_mean >
{
    using impl = accumulators::impl::distribution_variance_accumulator< mpl::_1 >;
};

} // namespace tag

namespace extract {

extractor< tag::distribution_variance > const distribution_variance = {};

} // namespace extract

using extract::distribution_variance;

// For the purposes of feature-based dependency resolution,
// distribution_variance provides the same feature as variance
template<> struct feature_of<tag::distribution_variance> : feature_of<tag::variance> {};

// enable use of tag::variance(distribution) as an alias of tag::distribution_variance
template<> struct as_feature<tag::variance(from_distribution)> { using type = tag::distribution_variance; };

// enable use of tag::variance(exact) as an alias of tag::variance
template<> struct as_feature<tag::variance(exact)> { using type = tag::variance; };

} // namspace accumulators
} // namspace boost

namespace apegrunt {
namespace accumulators {

template< typename Iterable, typename RealT=double >
RealT distribution_variance( Iterable range )
{
	using std::pow;
	using std::max;
	using std::cbegin; using std::cend;
	using result_type = RealT;
	using boost::accumulators::access_bin;
	using boost::accumulators::mean;
	using boost::accumulators::variance;
	using boost::accumulators::count;

	const result_type varsum = std::accumulate( cbegin(range), cend(range), result_type(0), [=]( auto sum, const auto it ) { return sum += ( variance(access_bin(it)) + pow(mean(access_bin(it)),2) ) * count(access_bin(it)); } );

	return ( varsum / result_type(max(std::size_t(1),distribution_count(range)-1)) ) - pow(distribution_mean(range),2);
}

} // accumulators
} // apegrunt

#endif // APEGRUNT_ACCUMULATORS_DISTRIBUTION_VARIANCE_HPP
