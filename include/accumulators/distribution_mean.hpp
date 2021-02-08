/** @file distribution_mean.hpp

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
#ifndef APEGRUNT_ACCUMULATORS_DISTRIBUTION_MEAN_HPP
#define APEGRUNT_ACCUMULATORS_DISTRIBUTION_MEAN_HPP

#include <cmath>

#include <boost/mpl/placeholders.hpp>
#include <boost/accumulators/framework/accumulator_base.hpp>
#include <boost/accumulators/framework/extractor.hpp>
#include <boost/accumulators/framework/depends_on.hpp>

#include <boost/accumulators/statistics/count.hpp>
#include <boost/accumulators/statistics/mean.hpp>

#include "accumulators/distribution.hpp"
#include "accumulators/distribution_count.hpp"

namespace boost {
namespace accumulators {

namespace impl {

template< typename Sample >
struct distribution_mean_accumulator : accumulator_base
{
	using result_type = Sample;

	template< typename Args >
	distribution_mean_accumulator( Args const & args ) {}

	template< typename Args >
	result_type result( Args const & args ) const
	{
		using std::max;
		using std::cbegin; using std::cend;

		const auto& d = distribution( args[accumulator] );
		return std::accumulate( cbegin(d), cend(d), result_type(0),
				[=]( auto sum, const auto it ) { return sum += mean(access_bin(it)) * count(access_bin(it)); }
		) / result_type( distribution_count(args[accumulator])-1 );
	}
};

} // namespace impl

namespace tag {

// distribution mean calculated from a binned set of samples
struct distribution_mean : depends_on< distribution_count >
{
	using impl = accumulators::impl::distribution_mean_accumulator< mpl::_1 >;
};

} // namespace tag

namespace extract {

extractor< tag::distribution_mean > const distribution_mean = {};

} // namespace extract

using extract::distribution_mean;

// For the purposes of feature-based dependency resolution,
// distribution_mean provides the same feature as mean
template<> struct feature_of<tag::distribution_mean> : feature_of<tag::mean> {};

// enable use of tag::mean(distribution) as an alias of tag::distribution_mean
template<> struct as_feature<tag::mean(from_distribution)> { using type = tag::distribution_mean; };

// enable use of tag::mean(exact) as an alias of tag::mean
template<> struct as_feature<tag::mean(exact)> { using type = tag::mean; };

} // namspace accumulators
} // namspace boost

namespace apegrunt {
namespace accumulators {

template< typename Iterable, typename RealT=double >
std::size_t distribution_mean( Iterable range )
{
	using result_type = RealT;
	using boost::accumulators::access_bin;
	using std::max;
	using std::cbegin; using std::cend;
	using boost::accumulators::mean;
	using boost::accumulators::count;

	return std::accumulate( cbegin(range), cend(range), result_type(0),
			[=]( auto sum, const auto it ) { return sum += mean(access_bin(it)) * count(access_bin(it)); }
	) / result_type( distribution_count(range)-1 );

}

} // accumulators
} // apegrunt

#endif // APEGRUNT_ACCUMULATORS_DISTRIBUTION_MEAN_HPP
