/** @file CouplingStorage.hpp

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
#ifndef APEGRUNT_TOOLS_COUPLING_STORAGE_HPP
#define APEGRUNT_TOOLS_COUPLING_STORAGE_HPP

#include "apegrunt/Apegrunt_options.h"
#include "apegrunt/aligned_allocator.hpp"
#include "misc/type_traits.hpp"
#include "misc/Array_view.hpp"
#include "misc/Math.hpp" // for ipow

namespace apegrunt {

//namespace tools {

uint64_t get_pool_size( uint64_t n_loci ) { return ipow(n_loci,2)-n_loci; }

template< typename RealT, uint States >
class CouplingStorage
{
public:

	enum { Q=States };

	using real_t = RealT;
	using internal_real_t = float;
	using allocator_t = typename apegrunt::memory::AlignedAllocator<internal_real_t>;

	using my_type = CouplingStorage<real_t,Q>;

	using raw_matrix_t = std::array< internal_real_t, Q*Q >;

	using array_view_t = Array_view< internal_real_t, Q >;
	using matrix_view_t = Array_view< array_view_t, extent<array_view_t>::value >;

	using matrix_view_array_t = Array_view< matrix_view_t >;

	CouplingStorage() : m_dim1(0), m_dim2(0), m_has_dim1_mapping(false), m_has_dim2_mapping(false) { }

	CouplingStorage( apegrunt::Loci_ptr dim1_loci, apegrunt::Loci_ptr dim2_loci )
	: m_dim1(dim1_loci->size()),
	  m_has_dim1_mapping(true),
	  m_dim2(dim2_loci->size()),
	  m_has_dim2_mapping(true)
	{
	    // transfer dim1 mapping
		std::size_t i=0;
		for( const auto locus: dim1_loci ) { m_dim1_mapping[locus] = i; ++i; }

		// transfer dim2 mapping
		std::size_t j=0;
		for( const auto locus: dim2_loci ) { m_dim2_mapping[locus] = j; ++j; }
	}

	bool allocate_matrix_pool()
	{
	    const uint64_t pool_size = m_dim1*m_dim2;

	    try { m_matrix_storage.resize( pool_size ); }
		catch(...)
		{
			*Apegrunt_options::get_err_stream() << "apegrunt error: unable to allocate " << my_type::memory_estimate(m_dim1,m_dim2) << " of memory\n\n";
			return false;
		}

		return true;
	}

	bool allocate_scalar_pool()
	{
	    const uint64_t pool_size = m_dim1*m_dim2;

	    try { m_coupling_storage.resize( pool_size ); }
		catch(...)
		{
			*Apegrunt_options::get_err_stream() << "apegrunt error: unable to allocate " << my_type::memory_estimate(m_dim1,m_dim2) << " of memory\n\n";
			return false;
		}

		return true;
	}

	static std::size_t memory_estimate( std::size_t dim1, std::size_t dim2, bool store_matrices=false )
	{
		return dim1*dim2 * ( store_matrices ? sizeof(raw_matrix_t) : sizeof(internal_real_t) );
	}

	static std::string memory_estimate_string( std::size_t dim1, std::size_t dim2, bool store_matrices=false )
	{
		std::ostringstream mem_string; mem_string << apegrunt::memory_string( my_type::memory_estimate(dim1,dim2,store_matrices) );
		return std::string( "coupling storage will require approximately " + mem_string.str() + " of memory" );
	}

	inline std::size_t get_matrix_storage_size() const { return m_matrix_storage.size(); }
	inline std::size_t get_coupling_storage_size() const { return m_coupling_storage.size(); }

	inline matrix_view_array_t get_Ji_matrices( std::size_t i ) // zero-based locus index; returns a view to a row vector
	{
		if( m_has_dim1_mapping )
		{
			i = m_dim1_mapping[i];
		}
		assert( i < m_dim1 );
		//return matrix_view_array_t( m_matrix_storage[i*(m_dim1-1)].data(), (m_dim2-1) );
		return matrix_view_array_t( m_matrix_storage[i*m_dim1].data(), m_dim2 );
	}

	inline matrix_view_t get_Jij_matrix( std::size_t i, std::size_t j ) // zero-based locus index
	{
		if( m_has_dim1_mapping && m_has_dim2_mapping )
		{
			i = m_dim1_mapping[i];
			j = m_dim2_mapping[j];
		}
		assert(  i < m_dim1 && j < m_dim2 );
		return this->get_Ji_matrices(i)[j]; //coupling_storage[i*(m_dim1-1)+j]; // we store diagonal elements, too, even though they might never be used.
	}

	inline internal_real_t& get_Jij_score( std::size_t i, std::size_t j ) // zero-based locus index
	{
		if( m_has_dim1_mapping )
		{
			i = m_dim1_mapping[i];
			j = m_dim2_mapping[j];
		}
		assert( i < m_dim1 && j < m_dim2 );
		//return m_coupling_storage[i*(m_dim1-1)+j]; // we store diagonal elements, too, even though they might never be used.
		return m_coupling_storage[i*m_dim1+j]; // we store diagonal elements, too, even though they might never be used.
	}

private:
    std::vector< raw_matrix_t > m_matrix_storage;
    std::vector< internal_real_t, allocator_t > m_coupling_storage;
    std::size_t m_dim1;
    std::size_t m_dim2;

    using loci_mapping_t = std::map<std::size_t,std::size_t>;

    bool m_has_dim1_mapping;
    loci_mapping_t m_dim1_mapping;
    bool m_has_dim2_mapping;
    loci_mapping_t m_dim2_mapping;

};

//} // namespace tools

} // namespace apegrunt

#endif // APEGRUNT_TOOLS_COUPLING_STORAGE_HPP
