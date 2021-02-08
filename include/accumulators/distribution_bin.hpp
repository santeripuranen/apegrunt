/** @file distribution_bin.hpp

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
#ifndef APEGRUNT_ACCUMULATORS_DISTRIBUTION_BIN_HPP
#define APEGRUNT_ACCUMULATORS_DISTRIBUTION_BIN_HPP

namespace apegrunt {
namespace accumulators {

template< typename Sample >
struct distribution_bin { std::size_t m_count; Sample m_variance; Sample m_mean; };

template< typename Sample >
std::size_t count( const distribution_bin<Sample>& bin ) { return bin.m_count; }

template< typename Sample >
Sample mean( const distribution_bin<Sample>& bin ) { return bin.m_mean; }

template< typename Sample >
Sample variance( const distribution_bin<Sample>& bin ) { return bin.m_variance; }

template< typename Sample >
Sample std( const distribution_bin<Sample>& bin ) { using std::sqrt; return sqrt(bin.m_variance); }

} // namspace accumulators
} // namspace apegrunt


#endif // APEGRUNT_ACCUMULATORS_DISTRIBUTION_BIN_HPP
