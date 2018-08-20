/** @file distribution_std.hpp

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
#ifndef APEGRUNT_ACCUMULATORS_DISTRIBUTION_STD_HPP
#define APEGRUNT_ACCUMULATORS_DISTRIBUTION_STD_HPP

#include <cmath>

#include <boost/mpl/placeholders.hpp>
#include <boost/accumulators/framework/accumulator_base.hpp>
#include <boost/accumulators/framework/extractor.hpp>
#include <boost/accumulators/framework/depends_on.hpp>

#include "accumulators/distribution_variance.hpp"

namespace boost {
namespace accumulators {

namespace impl {

template< typename Sample >
struct distribution_std_accumulator : accumulator_base
{
	using result_type = Sample;

	template< typename Args >
	distribution_std_accumulator( Args const & args ) {}

	template< typename Args >
	result_type result( Args const & args ) const
	{
		using std::sqrt;

		return sqrt( distribution_variance(args[accumulator]) );
	}
};

} // namespace impl

namespace tag {

// for convenience: distribution standard deviation calculated from a binned set of samples
struct distribution_std : depends_on< distribution_variance >
{
    using impl = accumulators::impl::distribution_std_accumulator< mpl::_1 >;
};

} // namespace tag

namespace extract {

extractor< tag::distribution_std > const distribution_std = {};

} // namespace extract

using extract::distribution_std;

// For the purposes of feature-based dependency resolution,
// distribution_std provides the same feature as std
template<> struct feature_of<tag::distribution_std> : feature_of<tag::std> {};

// enable use of tag::std(from_distribution) as an alias of tag::distribution_std
template<> struct as_feature<tag::std(from_distribution)> { typedef tag::distribution_std type; };

// enable use of tag::std(exact) as an alias of tag::std
template<> struct as_feature<tag::std(exact)> { typedef tag::std type; };

} // namspace accumulators
} // namspace boost

namespace apegrunt {
namespace accumulators {

template< typename Iterable, typename RealT=double >
RealT distribution_std( Iterable range ) { using std::sqrt; return sqrt( distribution_variance(range) ); }

} // accumulators
} // apegrunt

#endif // APEGRUNT_ACCUMULATORS_DISTRIBUTION_STD_HPP
