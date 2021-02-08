/** @file distribution_cumulative.hpp

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
#ifndef APEGRUNT_ACCUMULATORS_DISTRIBUTION_CUMULATIVE_HPP
#define APEGRUNT_ACCUMULATORS_DISTRIBUTION_CUMULATIVE_HPP

#include <algorithm> // for std::transform

#include <boost/mpl/placeholders.hpp>
#include <boost/accumulators/framework/accumulator_base.hpp>
#include <boost/accumulators/framework/extractor.hpp>
#include <boost/accumulators/framework/parameters/sample.hpp>
#include <boost/accumulators/framework/depends_on.hpp>

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/count.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/variance.hpp>

#include "accumulators/distribution.hpp"
#include "accumulators/distribution_bin.hpp"
#include "accumulators/distribution_bincount.hpp"

namespace apegrunt {
namespace accumulators {

template< typename Sample, typename ... Args >
std::vector< distribution_bin<Sample> > distribution_cumulative( boost::accumulators::accumulator_set< Sample, Args... >& acc )
{
	using value_type = distribution_bin<Sample>;
	using result_type = std::vector< value_type >;
	using std::cbegin; using std::cend; using std::begin;
	using boost::accumulators::count;
	using boost::accumulators::mean;
	using boost::accumulators::variance;
	using boost::accumulators::distribution;
	using boost::accumulators::distribution_bincount;

	const auto& d = distribution( acc );
	result_type bins( distribution_bincount( acc ) );

	std::transform( cbegin(d), cend(d), rbegin(bins),
		[=]( value_type local, const auto node ) {
			const auto& dbin = access_bin(node);
			local.m_count += count(dbin);
			local.m_variance = variance(dbin);
			local.m_mean = mean(dbin);
			return local;
		}
	);
	return bins;
};

} // namspace accumulators
} // namspace apegrunt


#endif // APEGRUNT_ACCUMULATORS_DISTRIBUTION_CUMULATIVE_HPP
