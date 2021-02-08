/** @file distribution_generator_csv.hpp

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
#ifndef APEGRUNT_ACCMULATORS_DISTRIBUTION_GENERATOR_CSV_HPP
#define APEGRUNT_ACCMULATORS_DISTRIBUTION_GENERATOR_CSV_HPP

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/count.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/variance.hpp>

#include "accumulators/distribution_generator_base.hpp"

namespace apegrunt {

namespace accumulators {

template< typename Iterable, typename Bin >
struct distribution_generator_csv : distribution_generator_base< distribution_generator_csv<Iterable,Bin>, Iterable, Bin >
{
	using base_type = distribution_generator_base< distribution_generator_csv<Iterable,Bin>, Iterable, Bin >;
	using base_type::base_type; // pull in base class constructor(s)

	inline void operator()( std::ostream& os ) const
	{
		using boost::accumulators::count;
		using boost::accumulators::mean;
		using boost::accumulators::variance;

		const auto dcount = distribution_count( base_type::range );
		const auto std = distribution_std( base_type::range );
		const double binwidth = 0.001;

		os.precision(3);
		for( const auto& node : base_type::range )
		{
			const auto& bin = base_type::access_bin(node);
			const auto cnt = count(bin);
			if( 0 != cnt )
			{
				// produce output for each populated bin
				os << mean(bin) << " " << cnt << " " << variance(bin) << " " << mean(bin)/std << " " << double(cnt)/double(dcount)*binwidth << "\n";
			}
		}
	}
};

// generic factory function matches std containers and boost::iterator_range
template< typename Iterable >
distribution_generator_csv<Iterable,typename Iterable::value_type> csv( Iterable&& range )
{
	return distribution_generator_csv<Iterable,typename Iterable::value_type>(std::move(range));
}

// specialized factory function for distribution_accumulator
template< typename Sample=double >
distribution_generator_csv< boost::iterator_range< typename boost::accumulators::impl::distribution_accumulator<Sample>::const_iterator >, typename boost::accumulators::impl::distribution_accumulator<Sample>::bin_type >
csv( boost::iterator_range< typename boost::accumulators::impl::distribution_accumulator<Sample>::const_iterator >&& range )
{
	namespace bai = boost::accumulators::impl;
	using boost::accumulators::access_bin;
	using range_type = boost::iterator_range< typename boost::accumulators::impl::distribution_accumulator<Sample>::const_iterator >;
	return distribution_generator_csv<range_type,typename bai::distribution_accumulator<Sample>::bin_type>(std::move(range),&access_bin<Sample>);
}

} // namespace accumulators

} // namespace apegrunt

#endif // APEGRUNT_ACCMULATORS_DISTRIBUTION_GENERATOR_CSV_HPP
