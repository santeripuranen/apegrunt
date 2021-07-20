/** @file Crosstable.hpp

	Copyright (c) 2016-2021 Santeri Puranen.

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

// Construct a contingency table (aka a crosstable)
template< std::size_t VectorSize, typename RealT >
class Weighted_crosstable_single_kernel
{
public:
	using real_t = RealT;
	enum { N=VectorSize };

	using vector_t = Vector<real_t,N>;
	using vector_view_t = Vector<real_t,N,true>;
	using my_type = Weighted_crosstable_single_kernel<VectorSize,real_t>;

	Weighted_crosstable_single_kernel() = delete;

	Weighted_crosstable_single_kernel( real_t* const data )
	: m_data(data)
	{ }

	~Weighted_crosstable_single_kernel() = default;
/*
	Weighted_crosstable_single_kernel( const my_type& other )
	: m_data(other.m_data)
	{ }

	my_type& operator=( const my_type& other )
	{
		m_data = other.m_data;
	}
*/
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
/*
// Construct a contingency table (aka a crosstable)
template< typename AccessOrder, std::size_t Size, typename RealT >
class Weighted_crosstable_1Dblock_kernel
{
public:
 	using real_t = RealT;
 	enum { N=AccessOrder::N };
 	enum { BlockSize=Size };

	using access_order_tag = AccessOrder;
 	using vector_t = Vector<real_t,N>;
 	using vector_view_t = Vector<real_t,N,true>;
	using my_type = Weighted_crosstable_1Dblock_kernel<AccessOrder,BlockSize,real_t>;

	Weighted_crosstable_1Dblock_kernel( real_t* const data, std::size_t extent )
	: m_data(data), m_extent(extent)
	{ }

	~Weighted_crosstable_1Dblock_kernel() = default;

	Weighted_crosstable_1Dblock_kernel( const my_type& other )
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
*/
// Construct a contingency table (aka a crosstable)
template< typename AccessOrder, std::size_t Size, typename RealT >
class Weighted_crosstable_2Dblock_kernel
{
public:
	using real_t = RealT;
	enum { N=AccessOrder::N };
	enum { BlockSize=Size };

	using access_order_tag = AccessOrder;
	using vector_t = Vector<real_t,N>;
	using vector_view_t = Vector<real_t,N,true>;

	using my_type = Weighted_crosstable_2Dblock_kernel<AccessOrder,BlockSize,real_t>;

	Weighted_crosstable_2Dblock_kernel() = delete;

	Weighted_crosstable_2Dblock_kernel( real_t* const data, std::size_t iextent, std::size_t jextent )
	: m_data(data), m_iextent(iextent), m_jextent(jextent)
	{ }

	~Weighted_crosstable_2Dblock_kernel() = default;
/*
	Weighted_crosstable_2Dblock_kernel( const my_type& other )
	: m_data(other.m_data), m_iextent(other.m_iextent), m_jextent(other.m_jextent)
	{ }

	my_type& operator=( const my_type& other )
	{
		m_data = other.m_data;
		m_iextent = other.m_iextent;
		m_jextent = other.m_jextent;
	}
*/
	template< typename StateT >
	inline void accumulate( const apegrunt::State_block<StateT,BlockSize>& istateblock, const apegrunt::State_block<StateT,BlockSize>& jstateblock, real_t weight )
	{
		const auto jcrement = N*N*m_jextent;
		for( std::size_t i=0; i < m_iextent; ++i )
		{
			const auto data = m_data + jcrement*i + std::size_t( istateblock[i] );
			for( std::size_t j=0; j < m_jextent; ++j )
			{
				const auto jstate = std::size_t( jstateblock[j] );
				//*( data + MATRICES_AccessOrder_tag<N>::ptr_increment(jstate,j,m_jextent) ) += weight;
				*( data + AccessOrder::ptr_increment(jstate,j,m_jextent) ) += weight;
			}
		}
	}

	template< typename StateT >
	inline void accumulate_iblocked( const apegrunt::State_block<StateT,BlockSize>& istateblock, const apegrunt::State_block<StateT,BlockSize>& jstateblock, std::size_t i0, std::size_t iN, real_t weight )
	{
		const auto jcrement = N*N*m_jextent;
		//for( std::size_t i=0; i < m_iextent; ++i )
		for( std::size_t i=i0; i < iN; ++i )
		{
			//const auto istate = std::size_t( istateblock[i] );
			//const auto data = m_data + N*N*m_jextent*i;
			const auto data = m_data + jcrement*i + std::size_t( istateblock[i] ); // access row-by-row
			for( std::size_t j=0; j < m_jextent; ++j )
			//for( std::size_t j=j0; j < jN; ++j )
			{
				const auto jstate = std::size_t( jstateblock[j] );
				//*( data + MATRICES_AccessOrder_tag<N>::ptr_increment(jstate,j,m_jextent) ) += weight;
				*( data + AccessOrder::ptr_increment(jstate,j,m_jextent) ) += weight;
			}
		}
	}

	template< typename StateT >
	//inline void accumulate( const apegrunt::State_block<StateT,BlockSize>& istateblock, const apegrunt::State_block<StateT,BlockSize>& jstateblock, std::size_t i0, std::size_t iN, real_t weight )
	inline void accumulate_jblocked( const apegrunt::State_block<StateT,BlockSize>& istateblock, const apegrunt::State_block<StateT,BlockSize>& jstateblock, std::size_t j0, std::size_t jN, real_t weight )
	{
		//std::cout << "[\n";
		const auto jcrement = N*N*m_jextent;
		for( std::size_t i=0; i < m_iextent; ++i )
		//for( std::size_t i=i0; i < iN; ++i )
		{
			//const auto istate = std::size_t( istateblock[i] );
			//const auto data = m_data + N*N*m_jextent*i;
			const auto data = m_data + jcrement*i + std::size_t( istateblock[i] ); // access row-by-row
			//for( std::size_t j=0; j < m_jextent; ++j )
			//std::cout << "  [\n";
			for( std::size_t j=j0; j < jN; ++j )
			{
				const auto jstate = std::size_t( jstateblock[j] );
				//*( data + MATRICES_AccessOrder_tag<N>::ptr_increment(jstate,j,m_jextent) ) += weight;
				*( data + AccessOrder::ptr_increment(jstate,j,m_jextent) ) += weight;
				//std::cout	<< "    i=" << i << " +" << jcrement*i + std::size_t( istateblock[i] )
				//			<< " j=" << j << " +" << AccessOrder::ptr_increment(jstate,j,m_jextent)
				//			<< " data=" << jcrement*i + std::size_t( istateblock[i] ) + AccessOrder::ptr_increment(jstate,j,m_jextent) << "\n";
			}
			//std::cout << "  ]\n";
		}
		//std::cout << "]\n";
	}

	inline void set( real_t value )
	{
		const auto jcrement = N*N*m_jextent;
		vector_t val(value);
		for( std::size_t i=0; i < m_iextent; ++i )
		{
			const auto data = m_data + jcrement*i;
			for( std::size_t j=0; j < m_jextent; ++j )
			{
				for( std::size_t state=0; state < N; ++state )
				{
					//vector_view_t( data + MATRICES_AccessOrder_tag<N>::ptr_increment(state,j,m_jextent) ) = val;
					vector_view_t( data + AccessOrder::ptr_increment(state,j,m_jextent) ) = val;
				}
			}
		}
	}

	inline void set_iblocked( real_t value, std::size_t i0, std::size_t iN )
	{
		const auto jcrement = N*N*m_jextent;
		vector_t val(value);
		for( std::size_t i=i0; i < iN; ++i )
		{
			const auto data = m_data + jcrement*i;
			for( std::size_t j=0; j < m_jextent; ++j )
			{
				for( std::size_t state=0; state < N; ++state )
				{
					//vector_view_t( data + MATRICES_AccessOrder_tag<N>::ptr_increment(state,j,m_jextent) ) = val;
					vector_view_t( data + AccessOrder::ptr_increment(state,j,m_jextent) ) = val;
				}
			}
		}
	}

	inline void set_jblocked( real_t value, std::size_t j0, std::size_t jN )
	{
		const auto jcrement = N*N*m_jextent;
		vector_t val(value);
		for( std::size_t i=0; i < m_iextent; ++i )
		//for( std::size_t i=i0; i < iN; ++i )
		{
			const auto data = m_data + jcrement*i;
			//for( std::size_t j=0; j < m_jextent; ++j )
			for( std::size_t j=j0; j < jN; ++j )
			{
				for( std::size_t state=0; state < N; ++state )
				{
					//vector_view_t( data + MATRICES_AccessOrder_tag<N>::ptr_increment(state,j,m_jextent) ) = val;
					vector_view_t( data + AccessOrder::ptr_increment(state,j,m_jextent) ) = val;
				}
			}
		}
	}

private:

	real_t* const m_data;
	const std::size_t m_iextent;
	const std::size_t m_jextent;
};

template< typename StateT, typename RealT >
class Weighted_crosstable_2Dblock
{
public:
	using state_t = StateT;
	using real_t = RealT;
	using my_type = Weighted_crosstable_2Dblock<state_t,real_t>;
	enum { BlockSize=apegrunt::StateBlock_size };
	enum { N=apegrunt::number_of_states<state_t>::value };
	using AccessOrder = apegrunt::MATRICES_AccessOrder_tag<N>;
	using vector_t = Vector<real_t,N>;

	Weighted_crosstable_2Dblock() = delete;

	Weighted_crosstable_2Dblock( const Alignment_ptr<state_t> alignment, const std::shared_ptr< std::vector<real_t> > weights, real_t initialize_with=0.0 )
	: m_alignment(alignment),
	  m_weights(weights),
	  m_initialize_with(initialize_with),
	  m_last_block_size(apegrunt::get_last_block_size(alignment->n_loci())),
	  m_last_block_index(apegrunt::get_last_block_index(alignment->n_loci()))
	{
	}

	Weighted_crosstable_2Dblock( const my_type& other )
	: m_alignment( other.m_alignment ),
	  m_weights( other.m_weights ),
	  m_initialize_with( other.m_initialize_with ),
	  m_last_block_size( other.m_last_block_size ),
	  m_last_block_index( other.m_last_block_index )
	{
	}

	~Weighted_crosstable_2Dblock() = default;

	inline void operator()( real_t* dest, std::size_t iblock, std::size_t jblock ) const
	{
		const auto& alignment = *m_alignment.get();
		const auto& blocks = *(alignment.get_block_storage());
		const auto& block_accounting = *(alignment.get_block_accounting());

		const auto& iblocks = blocks[iblock];
		const auto& jblocks = blocks[jblock];

		// block column intersections code
/*
		{ // all in one go
			//const auto intersection = apegrunt::block_list_intersection_weights( block_accounting[iblock], block_accounting[jblock], *m_weights );
			const auto intersection = apegrunt::block_list_intersection_weights_dynamic( block_accounting[iblock], block_accounting[jblock], *m_weights );

			const auto isize = ( iblock == m_last_block_index ? m_last_block_size : BlockSize );
			const auto jsize = ( jblock == m_last_block_index ? m_last_block_size : BlockSize );

			//using wcrosstable_kernel_t = Weighted_crosstable_2Dblock_kernel<MATRICES_AccessOrder_tag<N>,BlockSize,real_t>;
			using wcrosstable_kernel_t = Weighted_crosstable_2Dblock_kernel<AccessOrder,BlockSize,real_t>;
			wcrosstable_kernel_t wcrosstable_block_kernel( dest, isize, jsize );
			wcrosstable_block_kernel.set(m_initialize_with);

			for( const auto& edge: intersection )
			{
				wcrosstable_block_kernel.accumulate( iblocks[edge.node1()], jblocks[edge.node2()], edge.weight() );
			}
		}

		else
*/		{ // split into two pieces

			//using wcrosstable_kernel_t = Weighted_crosstable_2Dblock_kernel<MATRICES_AccessOrder_tag<N>,BlockSize,real_t>;
			using wcrosstable_kernel_t = Weighted_crosstable_2Dblock_kernel<AccessOrder,BlockSize,real_t>;

			{
				const auto intersection = apegrunt::block_list_intersection_weights_dynamic( block_accounting[iblock], block_accounting[jblock], *m_weights );

				const auto isize = ( iblock == m_last_block_index ? m_last_block_size : BlockSize );
				const auto jsize = ( jblock == m_last_block_index ? m_last_block_size : BlockSize );

				wcrosstable_kernel_t wcrosstable_block_kernel( dest, isize, jsize );

				const std::size_t last = jsize;
				const std::size_t midpoint = std::min(std::size_t(BlockSize/2),last);

				{ // first half
					wcrosstable_block_kernel.set_jblocked(m_initialize_with, 0, midpoint );

					for( const auto& edge: intersection )
					{
						wcrosstable_block_kernel.accumulate_jblocked( iblocks[edge.node1()], jblocks[edge.node2()], 0, midpoint, edge.weight() );
					}
				}
				{ // second half
					wcrosstable_block_kernel.set_jblocked(m_initialize_with, midpoint, last );

					for( const auto& edge: intersection )
					{
						wcrosstable_block_kernel.accumulate_jblocked( iblocks[edge.node1()], jblocks[edge.node2()], midpoint, last, edge.weight() );
					}
				}
			}
		}
// */
	}

	inline void single( real_t* dest, std::size_t ipos, std::size_t jpos ) const
	{
		const auto& alignment = *m_alignment.get();
		const auto& m_blocks = *(alignment.get_block_storage());
		const auto& m_accounting = *(alignment.get_block_accounting());

		const auto iblock = apegrunt::get_block_index(ipos);
		const auto jblock = apegrunt::get_block_index(jpos);

		using block_type = typename apegrunt::template Alignment<state_t>::block_type::template sub_block_t<1>;
		using accounting_t = typename apegrunt::template Alignment<state_t>::block_accounting_container_t;

		std::vector< block_type > a_bs, b_bs;
		std::vector< accounting_t > a_as, b_as;

		{
			// It turns out that it's on average much faster to isolate each of the interacting data columns, than
			// it is to just straight intersect the indexes of the entire column blocks.
			// The time saved by performing fewer intersect ops, due to the reduced number of data blocks, more than
			// makes up for the extraction overhead.
			using union_t = apegrunt::lazy_set_union<accounting_t>;
			using tst_type = apegrunt::TernarySearchTree< block_type, union_t /*, TST_block_size<state_t>::value*/ >;
			using boost::get;

			{
				tst_type a;
				const auto ipos_local = apegrunt::get_pos_in_block( ipos );

				for( auto&& b_a: apegrunt::zip_range(m_blocks[iblock],m_accounting[iblock]) ) { a.insert( get<0>(b_a).template block<1>(ipos_local) )->value |= get<1>(b_a); }
				auto fa( [&a_as, &a_bs](auto node) mutable { a_as.emplace_back(std::move(node->value)); a_bs.emplace_back(std::move(node->key)); } );
				a.for_each( fa );
			}
			{
				tst_type b;
				const auto jpos_local = apegrunt::get_pos_in_block( jpos );

				for( auto&& b_a: apegrunt::zip_range(m_blocks[jblock],m_accounting[jblock]) ) { b.insert( get<0>(b_a).template block<1>(jpos_local) )->value |= get<1>(b_a); }
				auto fb( [&b_as, &b_bs](auto node) mutable { b_as.emplace_back(std::move(node->value)); b_bs.emplace_back(std::move(node->key)); } );
				b.for_each( fb );
			}
		}

		//const auto intersection = apegrunt::inplace_block_list_intersection_weights( a_as, b_as, *m_weights );
		const auto intersection = apegrunt::inplace_block_list_intersection_weights_dynamic( a_as, b_as, *m_weights ); // a_as and b_as are modified in-place: they will be empty when done

		Weighted_crosstable_single_kernel<N,real_t> wcrosstable_block_kernel(dest);
		wcrosstable_block_kernel.set(m_initialize_with);

		for( const auto& edge: intersection )
		{
			wcrosstable_block_kernel.accumulate( a_bs[edge.node1()][0], b_bs[edge.node2()][0], edge.weight() ); // access the first, and only, elements of blocks
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

