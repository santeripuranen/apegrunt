/** @file Coincidence_matrix.hpp

	Copyright (c) 2016-2018 Santeri Puranen.

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

#ifndef APEGRUNT_COINCIDENCE_MATRIX_HPP
#define APEGRUNT_COINCIDENCE_MATRIX_HPP

#include <vector>

#include "Vector.h"
#include "Matrix_kernel_access_order.hpp"

namespace apegrunt {

// Construct a coincidence matrix (aka a crosstable)
template< std::size_t VectorSize, typename RealT >
class Weighted_coincidence_single_kernel
{
public:
	using real_t = RealT;
	enum { N=VectorSize };

	using vector_t = Vector<real_t,N>;
	using vector_view_t = Vector<real_t,N,true>;
	using my_type = Weighted_coincidence_single_kernel<VectorSize,real_t>;

	Weighted_coincidence_single_kernel( real_t* const data )
	: m_data(data)
	{ }

	~Weighted_coincidence_single_kernel() = default;

	Weighted_coincidence_single_kernel( const my_type& other )
	: m_data(other.m_data)
	{ }

	my_type& operator=( const my_type& other )
	{
		m_data = other.m_data;
	}

	template< typename StateT >
	inline void accumulate( StateT istate, StateT jstate, real_t weight )
	{
		const auto i = std::size_t( istate );
		const auto j = std::size_t( jstate );
		*( m_data + j*N + i ) += weight;
	}

	inline void set( real_t value )
	{
		vector_t val(value);
		for( std::size_t state=0; state < N; ++state )
		{
			vector_view_t( m_data + state*N ) = val;
		}
	}

private:

	real_t* const m_data;
};

// Construct a coincidence matrix (aka a crosstable)
template< typename AccessOrder, std::size_t Size, typename RealT >
class Weighted_coincidence_1Dblock_kernel
{
public:

	using real_t = RealT;
	enum { N=AccessOrder::N };
	enum { BlockSize=Size };

	using vector_t = Vector<real_t,N>;
	using vector_view_t = Vector<real_t,N,true>;

	using my_type = Weighted_coincidence_1Dblock_kernel<AccessOrder,BlockSize,real_t>;

	Weighted_coincidence_1Dblock_kernel( real_t* const data, std::size_t extent )
	: m_data(data), m_extent(extent)
	{ }

	~Weighted_coincidence_1Dblock_kernel() = default;

	Weighted_coincidence_1Dblock_kernel( const my_type& other )
	: m_data(other.m_data), m_extent(other.m_extent)
	{ }

	my_type& operator=( const my_type& other )
	{
		m_data = other.m_data;
		m_extent = other.m_extent;
	}

	template< typename StateT, typename VectorT >
	inline void accumulate( apegrunt::State_block<StateT,BlockSize> stateblock, const VectorT& vector )
	{
		for( std::size_t i=0; i < m_extent; ++i )
		{
			const auto state = std::size_t( stateblock[i] );
			//std::cout << "i=" << i << " state=" << state << " incr=" << AccessOrder::ptr_increment(state,i,m_extent) << " " << vector << "\n";
			vector_view_t( m_data + AccessOrder::ptr_increment(state,i,m_extent) ) += vector;
		}
	}

	template< typename StateT, typename VectorT >
	inline void accumulate( apegrunt::State_block<StateT,BlockSize> stateblock, const VectorT& vector, std::size_t exclude )
	{
		for( std::size_t i=0; i < m_extent; ++i )
		{
			if( i != exclude )
			{
				const auto state = std::size_t( stateblock[i] );
				vector_view_t( m_data + AccessOrder::ptr_increment(state,i,m_extent) ) += vector;
			}
		}
	}

	inline void clear()
	{
		vector_t zero; // default-constructed Vector is all zeros
		for( std::size_t i=0; i < m_extent; ++i )
		{
			for( std::size_t state=0; state < N; ++state )
			{
				vector_view_t( m_data + MATRICES_AccessOrder_tag<N>::ptr_increment(state,i,m_extent) ) = zero;
			}
		}
	}

	inline void set( real_t value )
	{
		vector_t val(value);
		for( std::size_t i=0; i < m_extent; ++i )
		{
			for( std::size_t state=0; state < N; ++state )
			{
				vector_view_t( m_data + MATRICES_AccessOrder_tag<N>::ptr_increment(state,i,m_extent) ) = val;
			}
		}
	}

private:

	real_t* const m_data;
	const std::size_t m_extent;
};

// Construct a coincidence matrix (aka a crosstable)
template< typename AccessOrder, std::size_t Size, typename RealT >
class Weighted_coincidence_2Dblock_kernel
{
public:
	using real_t = RealT;
	enum { N=AccessOrder::N };
	enum { BlockSize=Size };

	using vector_t = Vector<real_t,N>;
	using vector_view_t = Vector<real_t,N,true>;

	using my_type = Weighted_coincidence_2Dblock_kernel<AccessOrder,BlockSize,real_t>;

	Weighted_coincidence_2Dblock_kernel( real_t* const data, std::size_t iextent, std::size_t jextent )
	: m_data(data), m_iextent(iextent), m_jextent(jextent)
	{ }

	~Weighted_coincidence_2Dblock_kernel() = default;

	Weighted_coincidence_2Dblock_kernel( const my_type& other )
	: m_data(other.m_data), m_iextent(other.m_iextent), m_jextent(other.m_jextent)
	{ }

	my_type& operator=( const my_type& other )
	{
		m_data = other.m_data;
		m_iextent = other.m_iextent;
		m_jextent = other.m_jextent;
	}

	template< typename StateT >
	inline void accumulate( apegrunt::State_block<StateT,BlockSize> istateblock, apegrunt::State_block<StateT,BlockSize> jstateblock, real_t weight )
	{
		for( std::size_t i=0; i < m_iextent; ++i )
		{
			const auto data = m_data + N*N*m_jextent*i + std::size_t( istateblock[i] );
			for( std::size_t j=0; j < m_jextent; ++j )
			{
				const auto jstate = std::size_t( jstateblock[j] );
				*( data + MATRICES_AccessOrder_tag<N>::ptr_increment(jstate,j,m_jextent) ) += weight;
			}
		}
	}

	inline void set( real_t value )
	{
		vector_t val(value);
		for( std::size_t i=0; i < m_iextent; ++i )
		{
			const auto data = m_data + N*N*m_jextent*i;
			for( std::size_t j=0; j < m_jextent; ++j )
			{
				for( std::size_t state=0; state < N; ++state )
				{
					vector_view_t( data + MATRICES_AccessOrder_tag<N>::ptr_increment(state,j,m_jextent) ) = val;
				}
			}
		}
	}

private:

	real_t* const m_data;
	const std::size_t m_iextent;
	const std::size_t m_jextent;
};


template< typename AccessOrder, std::size_t StateBlockSize, typename RealT >
class Matrix_storage_view
{
public:
	using real_t = RealT;
	enum { N=AccessOrder::N };
	enum { BlockSize=StateBlockSize };

	using my_type = Matrix_storage_view<AccessOrder,BlockSize,real_t>;

	Matrix_storage_view( real_t* const data, std::size_t extent ) : m_data(data), m_extent(extent)
	{
	}
	~Matrix_storage_view() { }

	inline real_t* data() { return m_data; }
	inline const real_t* data() const { return m_data; }
	inline std::size_t size() const { return m_extent; }

	inline auto get_weighted_coincidence_single_kernel()
	{
		using wcoincidence_kernel_t = Weighted_coincidence_single_kernel<N,real_t>;

		return wcoincidence_kernel_t( m_data );
	}
/*
	inline auto get_weighted_coincidence_1Dblock_kernel( std::size_t bn, std::size_t block_size )
	{
		using wcoincidence_kernel_t = Weighted_coincidence_1Dblock_kernel<AccessOrder,BlockSize,real_t>;

		//real_t *const data = m_data+bn*BlockSize*N*N; // need the BlockSize here, not block_size
		real_t *const data = m_data+bn*block_size*N*N; // ..well, it depends on how the calling/called code maps the data, and..
		return wcoincidence_kernel_t( data, block_size ); // ..should be reflected here, too
	}
*/
	inline auto get_weighted_coincidence_2Dblock_kernel( std::size_t iblock_size, std::size_t jblock_size )
	{
		using wcoincidence_kernel_t = Weighted_coincidence_2Dblock_kernel<AccessOrder,BlockSize,real_t>;

		return wcoincidence_kernel_t( m_data, iblock_size, jblock_size );
	}

private:
	real_t* const m_data;
	const std::size_t m_extent;
};
/*
template< typename RealT >
struct stateweight { std::size_t state; RealT weight; };

template< typename StateT, typename RealT >
class weighted_coincidence_1Dblock
{
public:
	using state_t = StateT;
	using real_t = RealT;
	using my_type = weighted_coincidence_1Dblock<state_t,real_t>;
	enum { BlockSize=apegrunt::StateBlock_size };
	enum { N=apegrunt::number_of_states<state_t>::value };
	using AccessOrder = apegrunt::MATRICES_AccessOrder_tag<N>;
	using vector_t = Vector<real_t,N>;

	weighted_coincidence_1Dblock( const Alignment_ptr<state_t> alignment, const std::shared_ptr< std::vector<real_t> > weights, real_t initialize_with=0.0 )
	: m_alignment(alignment),
	  m_weights(weights),
	  m_r_states(alignment->size()),
	  m_initialize_with(initialize_with),
	  m_n_loci(alignment->n_loci()),
	  m_last_block_size(apegrunt::get_last_block_size(m_n_loci)),
	  m_last_block_index(apegrunt::get_last_block_index(m_n_loci)),
	  m_cached_block_weights(),
	  m_cached_block_col(0)
	{
	}

	weighted_coincidence_1Dblock( const my_type& other )
	: m_alignment( other.m_alignment ),
	  m_weights( other.m_weights ),
	  m_r_states( other.m_alignment->size() ),
	  m_initialize_with( other.m_initialize_with ),
	  m_n_loci( other.m_n_loci ),
	  m_last_block_size( other.m_last_block_size ),
	  m_last_block_index( other.m_last_block_index ),
	  m_cached_block_weights( other.m_cached_block_weights ),
	  m_cached_block_col( other.m_cached_block_col )
	{
	}

	void operator()( real_t* dest, std::size_t ref_locus, std::size_t block_col ) const
	{
		const auto& alignment = *m_alignment.get();

		const auto isize = ( apegrunt::get_block_index(ref_locus) == m_last_block_index ? m_last_block_size : BlockSize );
		const auto jsize = ( block_col == m_last_block_index ? m_last_block_size : BlockSize );

		const auto& blocks = (*m_alignment->get_block_storage())[block_col];
		const auto& sequence_indices = (*m_alignment->get_block_accounting())[block_col];

		for( std::size_t i=0; i < alignment.size(); ++i )
		{
			m_r_states[i] = std::size_t( (*alignment[i])[ref_locus] );
		}

		const auto& block_weights = this->get_block_weights( block_col );

		Matrix_storage_view<MATRICES_AccessOrder_tag<N>,BlockSize,real_t> wcoincidence_matrix_storage( dest, BlockSize*N*N );
		auto&& wcoincidence_block_kernel = wcoincidence_matrix_storage.get_weighted_coincidence_1Dblock_kernel( 0, jsize );
		//wcoincidence_block_kernel.clear();
		wcoincidence_block_kernel.set(m_initialize_with);

		for( std::size_t i=0; i < blocks.size(); ++i )
		{
			wcoincidence_block_kernel.accumulate( blocks[i], block_weights[i] );
		}
	}

private:
	const Alignment_ptr<state_t> m_alignment;
	const std::shared_ptr< std::vector< real_t > > m_weights;
	mutable std::vector<std::size_t> m_r_states;
	const real_t m_initialize_with;
	const std::size_t m_n_loci;
	const std::size_t m_last_block_size;
	const std::size_t m_last_block_index;

	mutable std::vector< vector_t > m_cached_block_weights;
	mutable std::size_t m_cached_block_col;

	const std::vector< vector_t >& get_block_weights( std::size_t block_col ) const
	{
		if( block_col != m_cached_block_col || m_cached_block_weights.empty() )
		{
			m_cached_block_col = block_col;
			m_cached_block_weights.clear();
			m_cached_block_weights.reserve( (*m_alignment->get_block_storage())[block_col].size() );

			const auto& sequence_indices = (*m_alignment->get_block_accounting())[block_col];
			auto& weights = *m_weights.get();

			for( const auto& sequence_indices_list: sequence_indices )
			{
				m_cached_block_weights.emplace_back(); // add new empty
				auto& wcache = m_cached_block_weights.back();
				for( const auto si: sequence_indices_list ) // loop over sequence indices for this block
				{
					wcache[ m_r_states[si] ] += weights[si];
				}
			}
		}
		return m_cached_block_weights;
	}

};
*/
template< typename StateT, typename RealT >
class weighted_coincidence_2Dblock
{
public:
	using state_t = StateT;
	using real_t = RealT;
	using my_type = weighted_coincidence_2Dblock<state_t,real_t>;
	enum { BlockSize=apegrunt::StateBlock_size };
	enum { N=apegrunt::number_of_states<state_t>::value };
	using AccessOrder = apegrunt::MATRICES_AccessOrder_tag<N>;
	using vector_t = Vector<real_t,N>;

	weighted_coincidence_2Dblock( const Alignment_ptr<state_t> alignment, const std::shared_ptr< std::vector<real_t> > weights, real_t initialize_with=0.0 )
	: m_alignment(alignment),
	  m_weights(weights),
	  m_initialize_with(initialize_with),
	  m_last_block_size(apegrunt::get_last_block_size(alignment->n_loci())),
	  m_last_block_index(apegrunt::get_last_block_index(alignment->n_loci()))
	{
	}

	weighted_coincidence_2Dblock( const my_type& other )
	: m_alignment( other.m_alignment ),
	  m_weights( other.m_weights ),
	  m_initialize_with( other.m_initialize_with ),
	  m_last_block_size( other.m_last_block_size ),
	  m_last_block_index( other.m_last_block_index )
	{
	}

	~weighted_coincidence_2Dblock() = default;

	inline void operator()( real_t* dest, std::size_t iblock, std::size_t jblock ) const
	{
		const auto& alignment = *m_alignment.get();
		const auto& blocks = *(alignment.get_block_storage());
		const auto& block_accounting = *(alignment.get_block_accounting());

		const auto& iblocks = blocks[iblock];
		const auto& jblocks = blocks[jblock];

		// block column intersections code
		//const auto intersection = apegrunt::block_list_intersection( block_accounting[iblock], block_accounting[jblock] );
		//const auto intersection_weights = apegrunt::block_list_intersection_weights( intersection, *m_weights );
		const auto intersection = apegrunt::block_list_intersection_weights( block_accounting[iblock], block_accounting[jblock], *m_weights );

		const auto isize = ( iblock == m_last_block_index ? m_last_block_size : BlockSize );
		const auto jsize = ( jblock == m_last_block_index ? m_last_block_size : BlockSize );

		Matrix_storage_view<MATRICES_AccessOrder_tag<N>,BlockSize,real_t> wcoincidence_matrix_storage( dest, BlockSize*BlockSize*N*N );
		auto&& wcoincidence_block_kernel = wcoincidence_matrix_storage.get_weighted_coincidence_2Dblock_kernel( isize, jsize );
		wcoincidence_block_kernel.set(m_initialize_with);

		for( std::size_t i=0; i < intersection.size(); ++i )
		{
			const auto& block_indices = intersection.block_pairs[i];
			//wcoincidence_block_kernel.accumulate( iblocks[block_indices.first], jblocks[block_indices.second], intersection_weights[i] );
			wcoincidence_block_kernel.accumulate( iblocks[block_indices.first], jblocks[block_indices.second], intersection.block_pair_weights[i] );
		}
	}

	inline void single( real_t* dest, std::size_t ipos, std::size_t jpos ) const
	{
		const auto& alignment = *m_alignment.get();
		const auto& blocks = *(alignment.get_block_storage());
		const auto& block_accounting = *(alignment.get_block_accounting());

		const auto iblock = apegrunt::get_block_index(ipos);
		const auto jblock = apegrunt::get_block_index(jpos);

		const auto& iblocks = blocks[iblock];
		const auto& jblocks = blocks[jblock];

		// block column intersections code
		//const auto intersection = apegrunt::block_list_intersection( block_accounting[iblock], block_accounting[jblock] );
		//const auto intersection_weights = apegrunt::block_list_intersection_weights( intersection, *m_weights );
		const auto intersection = apegrunt::block_list_intersection_weights( block_accounting[iblock], block_accounting[jblock], *m_weights );

		Matrix_storage_view<MATRICES_AccessOrder_tag<N>,BlockSize,real_t> wcoincidence_matrix_storage( dest, N*N );
		auto&& wcoincidence_block_kernel = wcoincidence_matrix_storage.get_weighted_coincidence_single_kernel();
		wcoincidence_block_kernel.set(m_initialize_with);

		const auto ipos_local = ipos%BlockSize;
		const auto jpos_local = jpos%BlockSize;

		for( std::size_t i=0; i < intersection.size(); ++i )
		{
			const auto& block_indices = intersection.block_pairs[i];
			//wcoincidence_block_kernel.accumulate( iblocks[block_indices.first][ipos_local], jblocks[block_indices.second][jpos_local], intersection_weights[i] );
			wcoincidence_block_kernel.accumulate( iblocks[block_indices.first][ipos_local], jblocks[block_indices.second][jpos_local], intersection.block_pair_weights[i] );
		}
	}

private:
	const Alignment_ptr<state_t> m_alignment;
	const std::shared_ptr< std::vector< real_t > > m_weights;
	const real_t m_initialize_with;
	const std::size_t m_last_block_size;
	const std::size_t m_last_block_index;
};

} // namespace apegrunt

#endif // APEGRUNT_COINCIDENCE_MATRIX_HPP

