/** @file distribution_ordered.hpp

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
#ifndef APEGRUNT_ACCUMULATORS_DISTRIBUTION_ORDERED_HPP
#define APEGRUNT_ACCUMULATORS_DISTRIBUTION_ORDERED_HPP

#include <algorithm> // for std::sort

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

namespace apegrunt {
namespace accumulators {

template< typename Sample >
struct bin_wrapper
{
	using wrapped_bin_type = typename boost::accumulators::impl::distribution_accumulator<Sample>::container_type::value_type;
	using bin_type = typename boost::accumulators::impl::distribution_accumulator<Sample>::bin_type;
	using my_type = bin_wrapper<Sample>;

	bin_wrapper( const wrapped_bin_type& bin ) : bin_ptr(&bin) { }
	const wrapped_bin_type *bin_ptr;

	// bin_wrappers are equal if the addresses of the wrapped bins are equal
	inline bool operator==( const my_type& rhs ) const { return bin_ptr == rhs.bin_ptr; }

	// *this is less than rhs if the bin range of the wrapped bin is less than that of wrapped bin in rhs
	inline bool operator<( const my_type& rhs ) const { return bin_ptr->first < rhs.bin_ptr->first; }

	static const wrapped_bin_type& access_bin( const my_type& vt ) { return boost::accumulators::access_bin(*(vt.bin_ptr)); }
};

template< typename Sample >
const typename bin_wrapper<Sample>::bin_type& access_bin( const bin_wrapper<Sample>& vt ) { bin_wrapper<Sample>::access_bin(vt); }

template< typename Sample, typename ... Args >
struct distribution_ordered_impl
{
	using wrapped_bin_type = bin_wrapper<Sample>;
	using bin_type = typename wrapped_bin_type::bin_type;
	using container_type = std::vector< wrapped_bin_type >;
	using value_type = typename  container_type::value_type;
	using const_iterator = typename container_type::const_iterator;

	distribution_ordered_impl() = delete;
	distribution_ordered_impl( boost::accumulators::accumulator_set< Sample, Args... >& acc ) :
		//m_distribution( distribution(acc) ),
		m_ordered( cbegin(distribution(acc)), cend(distribution(acc)) )
	{
		using std::begin; using std::end;
		std::sort( begin(m_ordered), end(m_ordered) );
	}

	const_iterator cbegin() const { using std::cbegin; return cbegin(m_ordered); }
	const_iterator cend() const { using std::cend; return cend(m_ordered); }

    static const bin_type& access_bin( const wrapped_bin_type& vt ) { return wrapped_bin_type::access_bin(vt); }

private:
	//const typename distribution_accumulator<Sample>::result_type m_distribution;
	container_type m_ordered;
};

template< typename Sample, typename ... Args >
const typename boost::accumulators::impl::distribution_accumulator<Sample>::bin_type& access_bin( const typename distribution_ordered_impl<Sample,Args...>::wrapped_bin_type& vt ) { return distribution_ordered_impl<Sample,Args...>::access_bin(vt); }

template< typename Sample, typename ... Args >
distribution_ordered_impl<Sample,Args...> distribution_ordered( boost::accumulators::accumulator_set< Sample, Args... >& acc ) { return distribution_ordered_impl<Sample,Args...>(acc); }

} // namspace accumulators
} // namspace apegrunt


#endif // APEGRUNT_ACCUMULATORS_DISTRIBUTION_ORDERED_HPP
