/** @file distribution_extractor.hpp

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
#ifndef APEGRUNT_ACCMULATORS_DISTRIBUTION_EXTRACTOR_HPP
#define APEGRUNT_ACCMULATORS_DISTRIBUTION_EXTRACTOR_HPP

namespace apegrunt {

namespace accumulators {
template< typename Feature, typename AccumulatorSet >
std::ostream& extract_result( std::ostream& os, AccumulatorSet const &acc )
{

}

namespace tag {

// for convenience: distribution standard deviation calculated from a binned set of samples
struct distribution_std : depends_on< distribution_variance >
{
    using impl = accumulators::impl::distribution_std_accumulator< mpl::_1 >;
};

} // namespace tag

namespace extractor {

template<typename Feature>
struct extractor
{
    typedef extractor<Feature> this_type;

    inline void operator()( std::ostream& os ) const
	{
		return static_cast<const_cast_t>(this)->operator()(os); // forward to
	}
};

} // namespace extractor

} // namespace accumulators

} // namespace apegrunt

#endif // APEGRUNT_ACCMULATORS_DISTRIBUTION_EXTRACTOR_HPP
