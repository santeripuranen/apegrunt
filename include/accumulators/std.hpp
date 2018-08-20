/** @file std.hpp

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
#ifndef APEGRUNT_ACCUMULATORS_STD_HPP
#define APEGRUNT_ACCUMULATORS_STD_HPP

#include <cmath>

#include <boost/mpl/placeholders.hpp>
#include <boost/accumulators/framework/accumulator_base.hpp>
#include <boost/accumulators/framework/extractor.hpp>
#include <boost/accumulators/framework/depends_on.hpp>

#include <boost/accumulators/statistics/variance.hpp>

namespace boost {
namespace accumulators {

namespace impl {

template< typename Sample >
struct std_accumulator : accumulator_base
{
	using result_type = Sample;

	template< typename Args >
	std_accumulator( Args const & args ) {}
/*
	template< typename Args >
	void operator()( dont_care ) {}
*/
	template< typename Args >
	result_type result( Args const & args ) const
	{
		using std::sqrt;
		return sqrt( variance(args[accumulator]) );
	}
};

} // namespace impl

namespace tag {

// for convenience: standard deviation accumulator
struct std : depends_on< variance(immediate) >
{
    using impl = accumulators::impl::std_accumulator< mpl::_1 >;
};

} // namespace tag

namespace extract {

extractor< tag::std > const std = {};

} // namespace extract

using extract::std;

} // namespace accumulators
} // namespace boost

#endif // APEGRUNT_ACCUMULATORS_STD_HPP
