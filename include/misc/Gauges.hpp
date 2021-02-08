/** @file Gauges.hpp

	Copyright (c) 2016-2017 Santeri Puranen.

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
#ifndef APEGRUNT_GAUGES_HPP
#define APEGRUNT_GAUGES_HPP

#include <array>

#include "misc/Matrix_math.hpp"
#include "misc/Coupling_matrix_view.hpp"

namespace apegrunt {

template< typename MatrixViewT >
auto ising_gauge( MatrixViewT&& J_ij, bool transpose_input=false )
{
	// J_ij: J_ij as estimated from from g_i.

	// Shift the coupling estimates into the Ising gauge.
	const auto Jt_ij = transpose(J_ij);

	if( transpose_input )
	{
	    const auto J_ij_mean = mean(Jt_ij);
	    const auto J_ij_mean_mean = mean(J_ij_mean);
	    const auto J_ij_mean2 = mean( J_ij );

	    return Jt_ij - repmat(J_ij_mean) - repmat(J_ij_mean2,true) + J_ij_mean_mean;
	}
	else
	{
		const auto J_ij_mean = mean(J_ij);
		const auto J_ij_mean_mean = mean(J_ij_mean);
		const auto J_ij_mean2 = mean( Jt_ij );

		return J_ij - repmat(J_ij_mean) - repmat(J_ij_mean2,true) + J_ij_mean_mean;
	}
}

template< typename AccessOrder, std::size_t BlockSize, typename RealT >
auto ising_gauge( apegrunt::Coupling_matrix_view<AccessOrder,BlockSize,RealT>& Js, std::size_t n, bool transpose_input=false )
{
	// J_ij: J_ij as estimated from from g_i.

	auto J_ij = convert( Js, n );

	// Shift the coupling estimates into the Ising gauge.
	const auto Jt_ij = transpose(J_ij);

	if( transpose_input )
	{
	    const auto J_ij_mean = mean(Jt_ij);
	    const auto J_ij_mean_mean = mean(J_ij_mean);
	    const auto J_ij_mean2 = mean( J_ij );

		return Jt_ij - repmat(J_ij_mean) - repmat(J_ij_mean2,true) + J_ij_mean_mean;
	}
	else
	{
		const auto J_ij_mean = mean(J_ij);
		const auto J_ij_mean_mean = mean(J_ij_mean);
		const auto J_ij_mean2 = mean( Jt_ij );

		return J_ij - repmat(J_ij_mean) - repmat(J_ij_mean2,true) + J_ij_mean_mean;
	}
}

} // namespace apegrunt

#endif // APEGRUNT_GAUGES_HPP

