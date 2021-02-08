/** @file Distance.hpp

	Copyright (c) 2019 Santeri Puranen.

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
#ifndef APEGRUNT_DISTANCE_HPP
#define APEGRUNT_DISTANCE_HPP

#include <cmath>

namespace apegrunt {

inline std::size_t linear_distance( std::size_t pos1, std::size_t pos2 )
{
    return std::abs( int64_t(pos1) - int64_t(pos2) );
}

inline std::size_t circular_distance( std::size_t pos1, std::size_t pos2, std::size_t circle_size )
{
	const auto d = linear_distance(pos1,pos2);
    return std::min( circle_size-d, d );
}

struct CircularDistance
{
	CircularDistance() = delete;
	CircularDistance( std::size_t circle_size ) : m_circle_size( circle_size ) { }

	inline std::size_t operator()( std::size_t pos1, std::size_t pos2 ) const
	{
		return circular_distance( pos1, pos2, m_circle_size );
	}

	const std::size_t m_circle_size;
};

struct LinearDistance
{
	LinearDistance() = delete;
	LinearDistance( std::size_t circle_size=0 ) { }

	inline std::size_t operator()( std::size_t pos1, std::size_t pos2 ) const
	{
		return linear_distance( pos1, pos2 );
	}
};

template< typename DistanceT >
struct GenomeDistance
{
	GenomeDistance() = delete;
	template< typename StateT >
	GenomeDistance( const Alignment_ptr<StateT> alignment )
	: m_distance( alignment->n_loci() ),
	  m_mapping_ptr( alignment->get_loci_translation() ),
	  m_mapping( *m_mapping_ptr )
	{
	}

	inline std::size_t operator()( std::size_t pos1, std::size_t pos2 ) const
	{
		return m_distance( m_mapping[pos1], m_mapping[pos2] );
	}

	const DistanceT m_distance;
	const apegrunt::Loci_ptr m_mapping_ptr;
	const apegrunt::Loci& m_mapping;
};

} // namespace apegrunt

#endif // APEGRUNT_DISTANCE_HPP
