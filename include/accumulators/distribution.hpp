/** @file distribution.hpp

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
#ifndef APEGRUNT_ACCUMULATORS_DISTRIBUTION_HPP
#define APEGRUNT_ACCUMULATORS_DISTRIBUTION_HPP

#include <cmath>
#include <numeric> // for std::accumulate
//#include <map>
#include <unordered_map>

#include <boost/mpl/placeholders.hpp>
#include <boost/parameter/keyword.hpp>
#include <boost/accumulators/framework/accumulator_base.hpp>
#include <boost/accumulators/framework/extractor.hpp>
#include <boost/accumulators/framework/parameters/sample.hpp>
#include <boost/accumulators/framework/depends_on.hpp>

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>

#include <boost/range/iterator_range.hpp>

#include "accumulators/std.hpp"

namespace boost {
namespace accumulators {

BOOST_PARAMETER_NESTED_KEYWORD(tag, distribution_binwidth, binwidth)
BOOST_ACCUMULATORS_IGNORE_GLOBAL(distribution_binwidth)

namespace impl {

/** distribution_accumulator -- a histogram distribution estimator similar
 * to boost::accumulators::density, but with extended functionalty.
 */

template< typename Sample >
struct distribution_accumulator : accumulator_base
{
	using my_type = distribution_accumulator<Sample>;
    using bin_type = accumulator_set<Sample, stats<tag::std> >;
    using key_type = int64_t;
    //using container_type = std::unordered_map<key_type,value_type>;
    using container_type = std::map<key_type,bin_type>;
    using const_iterator = typename container_type::const_iterator;
    using result_type = boost::iterator_range< typename my_type::const_iterator >;

    template< typename Args >
    distribution_accumulator( Args const & args )
    	: //m_bins( cbegin(args[ distribution ]), cend(args[ distribution ]) ) ),
		  m_binwidth( args[distribution_binwidth|Sample(0.001)] )
    {
    	// 1.0/m_binwidth <- assume input values are in interval [0,1]
    	//m_bins.reserve( 1.0/m_bins.max_load_factor()*1.0/m_binwidth+std::numeric_limits<Sample>::epsilon() );
    }

    distribution_accumulator( const my_type& other )
    	: m_bins(other.m_bins),
		  m_binwidth(other.m_binwidth)
    {
    }

    template< typename Args >
    void operator ()( Args const & args )
    {
    	const Sample& sample_ = args[sample];
        const auto key = key_type(std::round(sample_ / m_binwidth)); // round makes binning behave nice on both sides of zero
        m_bins[key](sample_);
    }

    result_type result(dont_care) const
    {
        return boost::make_iterator_range( this->cbegin(), this->cend() );
    }

    static const bin_type& access_bin( const typename container_type::value_type& vt ) { return vt.second; }

    Sample binwidth() const { return m_binwidth; }

private:
    container_type m_bins;
    const Sample m_binwidth;

    const_iterator cbegin() const { using std::cbegin; return cbegin(m_bins); }
    const_iterator cend() const { using std::cend; return cend(m_bins); }
};

template< typename Accumulator >
Accumulator join( const Accumulator& lhs, const Accumulator& rhs ) {}

} // namespace impl

template< typename T >
const T& access_bin( const T& t ) { return t; }

template< typename Sample=double >
const typename impl::distribution_accumulator<Sample>::bin_type& access_bin( const typename impl::distribution_accumulator<Sample>::container_type::value_type& vt ) { return impl::distribution_accumulator<Sample>::access_bin(vt); }
/*
// gcc has trouble deducing template parameter Sample, so here are the two most typical overloads of the above template function
const typename impl::distribution_accumulator<double>::bin_type& get_bin( const typename impl::distribution_accumulator<double>::container_type::value_type& vt ) { return impl::distribution_accumulator<double>::get_bin(vt); }
const typename impl::distribution_accumulator<float>::bin_type& get_bin( const typename impl::distribution_accumulator<float>::container_type::value_type& vt ) { return impl::distribution_accumulator<float>::get_bin(vt); }

template< typename Sample >
typename impl::distribution_accumulator<Sample>::bin_type& get_bin( typename impl::distribution_accumulator<Sample>::container_type::value_type& vt ) { return impl::distribution_accumulator<Sample>::get_bin(vt); }

// gcc has trouble deducing template parameter Sample, so here are the two most typical overloads of the above template function
typename impl::distribution_accumulator<double>::bin_type& get_bin( typename impl::distribution_accumulator<double>::container_type::value_type& vt ) { return impl::distribution_accumulator<double>::get_bin(vt); }
typename impl::distribution_accumulator<float>::bin_type& get_bin( typename impl::distribution_accumulator<float>::container_type::value_type& vt ) { return impl::distribution_accumulator<float>::get_bin(vt); }
*/

namespace tag {

// accumulate the shape of a distribution by sample distribution
struct distribution : depends_on< >, distribution_binwidth
{
	using impl = accumulators::impl::distribution_accumulator< mpl::_1 >;
    static boost::parameter::keyword< distribution_binwidth > const binwidth;
};

} // namespace tag

namespace extract {

extractor< tag::distribution > const distribution = {};
BOOST_ACCUMULATORS_IGNORE_GLOBAL(distribution)

} // namespace extract

using extract::distribution;

// tag for feature(distribution) aliasing of features depending on binning
struct from_distribution {};
struct exact {};

} // namespace accumulators
} // namespace boost

#endif // APEGRUNT_ACCUMULATORS_DISTRIBUTION_HPP
