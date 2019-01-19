/** @file Alignment_impl_block_compressed_storage.hpp

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

#ifndef APEGRUNT_ALIGNMENT_IMPL_BLOCK_COMPRESSED_STORAGE_HPP
#define APEGRUNT_ALIGNMENT_IMPL_BLOCK_COMPRESSED_STORAGE_HPP

//#include <iosfwd>
#include <cmath>
#include <cfloat>
#include <memory> // for std::enable_shared_from_this
#include <algorithm> // for std::find
#include <numeric> // for std::accumulate
#include <iterator> // for std::distance
#include <mutex> // for std::mutex

#include <boost/range/adaptor/indexed.hpp>

#include "Alignment_forward.h"
#include "Alignment_impl_base.hpp"
#include "StateVector_forward.h"

#include "Alignment_parser_forward.h"
#include "Alignment_generator_forward.h"

#include "Alignment_iterator.h"
//#include "Alignment_iterator_impl_block_compressed_storage_forward.h"
#include "Alignment_iterator_impl_block_compressed_storage.hpp"

#include "State_block.hpp"

#include "accumulators/distribution_std.hpp"
#include "accumulators/distribution_bincount.hpp"
#include "accumulators/distribution_generator_csv.hpp"

#include "Apegrunt_utility.hpp"
#include "misc/Math.hpp"

namespace apegrunt {

template< typename T >
std::size_t bytesize( const std::vector<T>& v ) { return v.size()*sizeof(T); }

template< typename T, std::size_t VectorSize, std::size_t StateBlockSize >
class Column_frequency_accumulator
{
public:
	using frequency_t = T;
	enum { N=VectorSize };
	enum { BlockSize=StateBlockSize };

	using my_type = Column_frequency_accumulator<frequency_t,N,BlockSize>;

	Column_frequency_accumulator( frequency_t* const data, std::size_t extent )
	: m_data(data), m_extent(extent)
	{ }

	~Column_frequency_accumulator() { }

	Column_frequency_accumulator( const my_type& other )
	: m_data(other.m_data), m_extent(other.m_extent)
	{ }

	my_type& operator=( const my_type& other )
	{
		m_data = other.m_data;
		m_extent = other.m_extent;
	}

	template< typename StateT, typename IndexContainerT >
	inline void accumulate( apegrunt::State_block<StateT,BlockSize> stateblock, const IndexContainerT& sequence_indices, frequency_t* const weights=nullptr )
	{
		//const auto n_seq = frequency_t( sequence_indices.size() );
		const frequency_t weight = weights
			? std::accumulate( cbegin(sequence_indices), cend(sequence_indices), frequency_t(0), [=]( frequency_t sum, const auto seqindex ) { return sum += ( *(weights+seqindex) ); } )
			: frequency_t( sequence_indices.size() );
		//std::cout << "weight=" << weight << std::endl;
		for( std::size_t i=0; i < m_extent; ++i )
		{
			const auto state = std::size_t( stateblock[i] );
			*(m_data+i*N+state) += weight;
		}
	}

	inline void setZero()
	{
		//vector_t zero;
		for( std::size_t i=0; i < m_extent; ++i )
		{
			for( std::size_t j=0; j < N; ++j )
			{
				*(m_data+i*N+j) = frequency_t(0);
			}
			//vector_view_t( m_data+i*N ) = zero();
		}
	}

private:

	frequency_t* const m_data;
	const std::size_t m_extent;
};

} // namespace apegrunt

namespace apegrunt {

// forward declarations

template< typename AlignmentT > class Alignment_factory;

template< typename StateVectorT > //=StateVector_impl_block_compressed_alignment_storage<nucleic_acid_state_t> >
class Alignment_impl_block_compressed_storage : public Alignment_impl_base< Alignment_impl_block_compressed_storage<StateVectorT>, typename StateVectorT::state_t >
{
public:
	using statevector_t = StateVectorT;
	using state_t = typename statevector_t::state_t;

	using my_type = Alignment_impl_block_compressed_storage<statevector_t>;
	using base_type = Alignment_impl_base<my_type,state_t>;
	using const_iterator = typename base_type::const_iterator;
	using iterator = typename base_type::iterator;
	using value_type = typename base_type::value_type;

	using block_type = typename StateVector<state_t>::block_type;
	enum { N=block_type::N };

	using frequency_t = typename base_type::frequency_t;
	using frequencies_t = typename base_type::frequencies_t;
	using frequencies_ptr = typename base_type::frequencies_ptr;

	using w_frequency_t = typename base_type::w_frequency_t;
	using w_frequencies_t = typename base_type::w_frequencies_t;
	using w_frequencies_ptr = typename base_type::w_frequencies_ptr;

	using distance_matrix_t = typename base_type::distance_matrix_t;
	using distance_matrix_ptr = typename base_type::distance_matrix_ptr;

	using block_index_t = typename base_type::block_index_t;
	using block_accounting_t = typename base_type::block_accounting_t;
	using block_accounting_ptr = typename base_type::block_accounting_ptr;

	using block_weight_t = typename base_type::block_weight_t;
	using block_weights_t = typename base_type::block_weights_t;
	using block_weights_ptr = typename base_type::block_weights_ptr;

	using block_storage_t = typename base_type::block_storage_t;
	using block_storage_ptr = typename base_type::block_storage_ptr;

	using block_indices_t = typename base_type::block_indices_t;
	using block_indices_ptr = typename base_type::block_indices_ptr;

	using statecount_t = typename base_type::statecount_t;
	using statecount_block_t = typename base_type::statecount_block_t;
	using statecount_block_storage_t = typename base_type::statecount_block_storage_t;
	using statecount_block_storage_ptr = typename base_type::statecount_block_storage_ptr;

	using statepresence_t = typename base_type::statepresence_t;
	using statepresence_block_t = typename base_type::statepresence_block_t;
	using statepresence_block_storage_t = typename base_type::statepresence_block_storage_t;
	using statepresence_block_storage_ptr = typename base_type::statepresence_block_storage_ptr;

	Alignment_impl_block_compressed_storage() : m_rows(), m_block_storage( std::make_shared<block_storage_t>() ) { }

	~Alignment_impl_block_compressed_storage() = default;

	Alignment_impl_block_compressed_storage( const my_type& other )
		: base_type( other.id_string() ), m_rows( other.m_rows ), m_block_storage( other.m_block_storage )
	{
	}

	Alignment_impl_block_compressed_storage( my_type&& other ) noexcept
		: base_type( other.id_string() ), m_rows( std::move(other.m_rows) ), m_block_storage(other.m_block_storage) //m_block_storage( std::move(other.m_block_storage) )
	{
	}

	my_type& operator=( const my_type& other )
	{
		this->set_id_string( other.id_string() );
		m_rows = other.m_rows;
		m_block_storage = other.m_block_storage;
		return *this;
	}

	my_type& operator=( my_type&& other ) noexcept
	{
		this->set_id_string( other.id_string() );
		m_rows = std::move( other.m_rows );
		m_block_storage = other.m_block_storage; // std::move( other.m_block_storage );
		return *this;
	}

	Alignment_ptr<state_t> clone() const
	{
		return make_Alignment_ptr( *this );
	}

    inline const_iterator cbegin() const { return const_iterator( std::make_shared<const_iterator_impl>( m_rows.cbegin() ) ); }
    inline const_iterator cend() const { return const_iterator( std::make_shared<const_iterator_impl>( m_rows.cend() ) ); }

    inline const_iterator begin() const { return this->cbegin(); }
    inline const_iterator end() const { return this->cend(); }

    inline iterator begin() { return iterator( std::make_shared<iterator_impl>( m_rows.begin() ) ); }
    inline iterator end() { return iterator( std::make_shared<iterator_impl>( m_rows.end() ) ); }

    inline value_type operator[]( std::size_t index ) const { return m_rows[index]; }

    iterator erase( iterator first, iterator last )
    {
    	using std::cbegin; using std::cend;
    	const auto beg = std::find( cbegin(m_rows), cend(m_rows), *first );
    	const auto end = std::find( beg, cend(m_rows), *last );
    	return iterator( std::make_shared<iterator_impl>( m_rows.erase(beg,end) ) );
    }

    inline std::size_t size() const { return m_rows.size(); }
    inline std::size_t n_loci() const { return m_rows.front()->size(); } // this may not always be safe/produce the correct result

    inline const std::type_info& type() const { return typeid(my_type); }

    inline frequencies_ptr frequencies() const
    {
		//std::cout << " {get frequencies"; std::cout.flush();
    	m_cache_frequencies_mutex.lock();
    	if( !m_frequencies ) { this->cache_column_frequencies(); }
    	m_cache_frequencies_mutex.unlock();
		//std::cout << "}"; std::cout.flush();
    	return m_frequencies;
    }

    inline w_frequencies_ptr w_frequencies() const
    {
		//std::cout << " {get w_frequencies"; std::cout.flush();
    	m_cache_w_frequencies_mutex.lock();
    	if( !m_w_frequencies ) { this->cache_column_w_frequencies(); }
    	m_cache_w_frequencies_mutex.unlock();
		//std::cout << "} " << std::endl;
    	return m_w_frequencies;
    }

    inline distance_matrix_ptr distance_matrix() const
    {
    	m_cache_distance_matrix_mutex.lock();
    	if( !m_distance_matrix ) { this->cache_distance_matrix(); }
    	m_cache_distance_matrix_mutex.unlock();
    	return m_distance_matrix;
    }

    inline Alignment_subscript_proxy< StateVector_ptr<state_t> > subscript_proxy() const { return Alignment_subscript_proxy< StateVector_ptr<state_t> >( &m_rows ); }

	block_accounting_ptr get_block_accounting() const
	{
		m_cache_block_accounting_mutex.lock();
		if( !m_block_accounting ) { this->cache_block_accounting(); }
		m_cache_block_accounting_mutex.unlock();
		return m_block_accounting;
	}

	block_weights_ptr get_block_weights() const
	{
		m_cache_block_weights_mutex.lock();
		if( !m_block_weights ) { this->cache_block_weights(); }
		m_cache_block_weights_mutex.unlock();
		return m_block_weights;
	}

	void get_column_block_distance_matrix() const
	{
		this->cache_column_block_distance();
	}

	block_storage_ptr get_block_storage() const
	{
		return m_block_storage;
	}

	statecount_block_storage_ptr get_statecount_blocks() const
	{
		//std::cout << "apegrunt: get_statecount_blocks()" << std::endl;
		m_cache_statecount_and_statepresence_blocks_mutex.lock();
		if( !m_statecount_blocks ) { this->cache_statecount_and_statepresence_blocks(); }
		m_cache_statecount_and_statepresence_blocks_mutex.unlock();
		//std::cout << "apegrunt: m_statecount_blocks->size()=" << m_statecount_blocks->size() << std::endl;
		//std::cout << "apegrunt: return m_statecount_blocks" << std::endl;
		return m_statecount_blocks;
	}

	statepresence_block_storage_ptr get_statepresence_blocks() const
	{
		//std::cout << "apegrunt: get_statepresence_blocks()" << std::endl;
		m_cache_statecount_and_statepresence_blocks_mutex.lock();
		if( !m_statepresence_blocks ) { this->cache_statecount_and_statepresence_blocks(); }
		m_cache_statecount_and_statepresence_blocks_mutex.unlock();
		//std::cout << "apegrunt: m_statepresence_blocks->size()=" << m_statepresence_blocks->size() << std::endl;
		//std::cout << "apegrunt: return m_statepresence_blocks" << std::endl;
		return m_statepresence_blocks;
	}

	block_indices_ptr get_block_indices() const
    {
		// this is a silly function.. should use something like boost::irange instead
		auto indices_ptr = std::make_shared<block_indices_t>();
		auto& indices = *indices_ptr;
		indices.reserve( m_block_storage->size() );
		for( std::size_t i=0; i < m_block_storage->size(); ++i ) { indices.push_back(i); }
		return indices_ptr;
    }

	void statistics( std::ostream *out=nullptr ) const
	{
		using std::cbegin; using std::cend;
		if( !out ) { return; }

		const std::size_t n_loci = this->n_loci();
		const std::size_t n_seqs = m_rows.size();
		const std::size_t n_block_column_groups = apegrunt::get_number_of_blocks(n_loci); // n_loci / N + ( n_loci % N == 0 ? 0 : 1 );
		const std::size_t indexing_overhead_mem = std::accumulate( cbegin(m_rows), cend(m_rows), std::size_t(0), [=]( std::size_t sum, const auto& seq ) { return sum += seq->bytesize(); } );
		const std::size_t dense_indexing_overhead_mem = sizeof(block_index_t)*n_block_column_groups*n_seqs;
/*
		std::size_t duals = 0;
		std::size_t full_dual_blocks = 0;
		std::size_t full_dual_cols = 0;
*/
		std::size_t n_cblocks = 0;
		std::size_t min_blocks_per_column = m_block_storage->front().size();
		std::size_t max_blocks_per_column = 0;
		std::vector<std::size_t> col_size; col_size.reserve(m_block_storage->size());
		//std::map<std::size_t,std::size_t> nb_bins;

		//std::set<block_type> global_blocks;
		std::vector<block_type> global_blocks;

		for( const auto& block_col: *m_block_storage )
		{
			const std::size_t nb = block_col.size();
			col_size.push_back( nb );

			min_blocks_per_column = std::min(min_blocks_per_column,nb);
			max_blocks_per_column = std::max(max_blocks_per_column,nb);
			n_cblocks += nb;
/*
			using std::cbegin; using std::cend;
			for( const auto& block: block_col )
			{
				auto pos = std::find( cbegin(global_blocks), cend(global_blocks), block );
				if( pos == cend(global_blocks) )
				{
					global_blocks.emplace_back( block );
				}
			}
*/
			// if using std::set
			//for( const auto& block: block_col ) { global_blocks.insert( block ); }
/*
			std::size_t fd_count = 0;
			for( const auto block: block_col )
			{
				std::size_t dcount = 0;
				for( std::size_t i = 0; i < N; ++i ) { dcount += ( i%2 != 0 ? (block[i-1] == block[i]) : 0 ); }
				duals += dcount;
				fd_count += (N/2 == dcount);
			}
			full_dual_blocks += fd_count;
*/
		}

		// Get state statistics
		using std::begin; using std::end;
		std::array<std::size_t,number_of_states<state_t>::value+1> state_stats; std::fill( begin(state_stats), end(state_stats), 0 );

		for( const auto& block: *(this->get_statecount_blocks()) )
		{
			//for( auto count: block ) { state_stats[ count ] += 1; }
			for( std::size_t i=0; i < N; ++i ) { state_stats[ block[i] ] += 1; }
		}

		std::size_t SNPs = 0;
		for( std::size_t i=2; i < state_stats.size(); ++i ) { SNPs += state_stats[i]; }
		// end Get state statistics

		const std::size_t n_gblocks = global_blocks.size();

		std::sort( col_size.begin(), col_size.end() );
		const std::size_t n_col_median = col_size[ col_size.size()/2 ];
		const std::size_t n_col_mean = n_cblocks / n_block_column_groups;

		const std::size_t onucs = n_loci*n_seqs;
		const std::size_t cnucs = n_cblocks*N;
		const std::size_t onucs_mem = onucs*sizeof(state_t);
		const std::size_t cnucs_mem = cnucs*sizeof(state_t);
		const std::size_t gnucs_mem = n_gblocks*N*sizeof(state_t);

		*out << "apegrunt: alignment has " << n_seqs << " samples and " << n_loci << " positions (" << onucs << " outcomes in total)\n";
		*out << "apegrunt: alignment has " << SNPs << " SNPs in total\n";
		*out << "apegrunt: state distribution [#states:#positions] =";
		for( std::size_t i=1; i < state_stats.size(); ++i )
		{
			if( state_stats[i] != 0 )
			{
				*out << " " << i << ":" << state_stats[i];
			}
		}
		*out << "\n";

//		*out << "apegrunt: stateblock compression uses " << n_block_column_groups << " column groups of " << N << "-nucleotide blocks (" << sizeof(block_type) << "B per block)\n";
		*out << "apegrunt: uncompressed size = " << memory_string(onucs_mem) << "\n";
// Compressed blocks -- block-wise lists
		*out << "apegrunt: compressed size = " << memory_string(cnucs_mem+indexing_overhead_mem) << " (" << memory_string(cnucs_mem) << " for " << n_cblocks << " blocks and " << memory_string(indexing_overhead_mem) << " for indexing)\n";
		*out << "apegrunt: compression ratio = " << double(onucs_mem)/double(cnucs_mem+indexing_overhead_mem) << "\n";
		*out << "apegrunt: min/mean/median/max # of stateblocks per column = " << min_blocks_per_column << "/" << n_col_mean << "/" << n_col_median << "/" << max_blocks_per_column << "\n";
//		*out << "apegrunt: globally compressed size = " << memory_string(gnucs_mem+indexing_overhead_mem) << " (" << memory_string(gnucs_mem) << " for " << n_gblocks << " globally unique blocks and " << memory_string(indexing_overhead_mem) << " for indexing)\n";
//		*out << "apegrunt: compression ratio = " << double(onucs_mem)/double(gnucs_mem+indexing_overhead_mem) << " (uncompressed/global) | compression ratio = " << double(n_cblocks)/double(n_gblocks) << " (compressed blocks/globally compressed blocks)" << "\n";
/*
		*out << "apegrunt: duals = " << duals << " (" << double(duals)/double(onucs/2) << ")\n";
		*out << "apegrunt: full dual blocks = " << full_dual_blocks << " (" << double(full_dual_blocks)/double(n_cblocks) << ")\n";
		*out << "apegrunt: full dual block columns = " << full_dual_cols << " (" << double(full_dual_cols)/double(n_block_column_groups) << ")\n";
*/
	}

private:
	using iterator_impl = apegrunt::iterator::Alignment_iterator_impl_block_compressed_storage< StateVector_ptr<state_t> >;
	using const_iterator_impl = apegrunt::iterator::Alignment_const_iterator_impl_block_compressed_storage< StateVector_ptr<state_t> >;

	// The data row proxies and the actual data matrix should in general not
	// be modified during the lifetime of a const Alignment..
	std::vector< StateVector_ptr<state_t> > m_rows;
	block_storage_ptr m_block_storage;

	// ..but all caches may be modified and are therefore declared mutable
	mutable block_accounting_ptr m_block_accounting;
	mutable block_weights_ptr m_block_weights;
	mutable frequencies_ptr m_frequencies;
	mutable w_frequencies_ptr m_w_frequencies;
	mutable distance_matrix_ptr m_distance_matrix;
	mutable statecount_block_storage_ptr m_statecount_blocks;
	mutable statecount_block_storage_ptr m_statepresence_blocks;

	mutable std::mutex m_cache_block_accounting_mutex;
	mutable std::mutex m_cache_block_weights_mutex;
	mutable std::mutex m_cache_frequencies_mutex;
	mutable std::mutex m_cache_w_frequencies_mutex;
	mutable std::mutex m_cache_distance_matrix_mutex;
	mutable std::mutex m_cache_statecount_and_statepresence_blocks_mutex;

	// pull in all
	void cache_block_accounting() const
	{
		namespace acc = boost::accumulators;

		//std::cout << "apegrunt: cache block accounting\n"; std::cout.flush();
		//using boost::get;
		using std::cbegin;
		using std::cend;
		m_block_accounting = std::make_shared<block_accounting_t>( m_block_storage->size() );

		std::size_t flat_mem = 0;
		std::size_t compr_mem = 0;

		// reserve memory
		for( std::size_t col=0; col < m_block_storage->size(); ++col )
		{
			(*m_block_accounting)[col].resize( (*m_block_storage)[col].size() );
		}

		// add block indices
		for( std::size_t seq=0; seq < this->size(); ++seq )
		{
			std::size_t col=0;
			for( const auto& block_index: m_rows[seq]->get_block_indices() )
			{
				(*m_block_accounting)[col][block_index].push_back(seq);
				++col;
			}
		}
/*
		acc::accumulator_set<double, acc::stats<acc::tag::std(acc::from_distribution),acc::tag::distribution_bincount> >
		blockcount_distribution( acc::tag::distribution::binwidth=1 );
*/
		// collect memory use statistics
		for( std::size_t col=0; col < m_block_storage->size(); ++col )
		{
			for( const auto& index_container: (*m_block_accounting)[col] )
			{
				//blockcount_distribution( index_container.size() );
				flat_mem += index_container.size()*sizeof(block_index_t);
				compr_mem += apegrunt::bytesize(index_container);
			}
		}
/*
		{ // block pair statistics -- beware: O(N^2) in number of block cols = very costly
			std::size_t matching_pairs = 0;
			std::size_t naive_pairs = 0;
			const std::size_t total_col_pairs = ( apegrunt::ipow( m_block_storage->size(), 2 ) - m_block_storage->size() ) / 2;
			std::size_t done_col_pairs = 0;

			std::cout << "apegrunt: find column pair matches" << std::endl;
			for( std::size_t col_i=0; col_i < m_block_storage->size(); ++col_i )
			{
				for( std::size_t col_j=0; col_j < col_i; ++col_j )
				{
					std::cout << "\r  processing.. " << double(done_col_pairs)/double(total_col_pairs)*100 << " % -- found " << matching_pairs << " / " << naive_pairs << " matching pairs (" << double(matching_pairs)/double(naive_pairs) << ")"; std::cout.flush();
					naive_pairs += (*m_block_accounting)[col_i].size() * (*m_block_accounting)[col_j].size();
					for( const auto& index_i_container: (*m_block_accounting)[col_i] )
					{
						for( const auto& index_j_container: (*m_block_accounting)[col_j] )
						{
							// reference implementation
							//const auto match = std::find_first_of( std::cbegin(index_i_container), std::cend(index_i_container), std::cbegin(index_j_container), std::cend(index_j_container) );
							//if( match != std::cend(index_i_container) ) { ++matching_pairs; }

							if( index_j_container.front() > index_i_container.back() || index_j_container.back() < index_i_container.front() ) { continue; }

							// our index lists are sorted, so a faster binary search is applicable
							if( index_i_container.size() < index_j_container.size() )
							{
								for( const auto& index: index_j_container )
								{
									if( index > index_i_container.back() ) { break; }
									if( index < index_i_container.front() ) { continue; }
									if( std::binary_search( std::cbegin(index_i_container), std::cend(index_i_container), index ) )
									{
										++matching_pairs;
										break;
									}
								}
							}
							else
							{
								for( const auto& index: index_i_container )
								{
									if( index > index_j_container.back() ) { break; }
									if( index < index_j_container.front() ) { continue; }
									if( std::binary_search( std::cbegin(index_j_container), std::cend(index_j_container), index ) )
									{
										++matching_pairs;
										break;
									}
								}
							}
						}
					}
					++done_col_pairs;
				}
			}

			std::cout << "\rapegrunt: found " << matching_pairs << " / " << naive_pairs << " matching pairs (" << double(matching_pairs)/double(naive_pairs)*100 << " %)" << std::endl;
		}
*/
/*
		auto distribution_file = apegrunt::get_unique_ofstream( "block_index_distribution.csv" );
		*distribution_file->stream() << std::fixed;
		*distribution_file->stream() << apegrunt::accumulators::csv(acc::distribution(blockcount_distribution));
*/
/*
		std::cout << " {level1"; std::cout.flush();
		for( std::size_t col=0; col < m_block_storage->size(); ++col )
		{
			(*m_block_accounting)[col].resize( (*m_block_storage)[col].size() );

			//std::cout << " {level2"; std::cout.flush();
			for( std::size_t seq=0; seq < this->size(); ++seq )
			{
				auto query_block = m_rows[seq]->get_block(col);
				auto block_itr = std::find( cbegin( (*m_block_storage)[col] ), cend( (*m_block_storage)[col] ), query_block );
				if( block_itr != cend( (*m_block_storage)[col] ) )
				{
					using std::distance;
					//std::cout << " {distance"; std::cout.flush();
					(*m_block_accounting)[col][ distance( cbegin( (*m_block_storage)[col] ), block_itr ) ].push_back( seq );
					//std::cout << "} ";
				}
			}
			//std::cout << "} ";
			for( const auto& index_container: (*m_block_accounting)[col] )
			{
				flat_mem += index_container.size()*sizeof(block_index_t);
				compr_mem += apegrunt::bytesize(index_container);
			}
		}
		std::cout << "} ";
*/
		//std::cout << "}"; std::cout.flush();
	}


	void cache_block_weights() const
	{
/*
		const auto& block_accounting = *(this->get_block_accounting());
		const auto& blocks = *(alignment->get_block_storage());

		m_block_weights = std::make_shared<block_weights_t>( block_accounting.size() );

		// reserve memory
		for( std::size_t col=0; col < m_block_weights->size(); ++col )
		{
			(*m_block_weights)[col].resize( block_accounting[col].size() );
		}

	    const std::size_t n_loci_per_block = apegrunt::StateBlock_size;
	    const std::size_t last_block_size = apegrunt::get_last_block_size(this->n_loci());
	    const std::size_t last_block = apegrunt::get_last_block_index(this->n_loci());
	    const std::size_t number_of_blocks = apegrunt::get_number_of_blocks(this->n_loci());

		for( std::size_t n_block=0; n_block < number_of_blocks; ++n_block )
		{
			const auto n_end = ( n_block == last_block ? last_block_size : n_loci_per_block );
			//auto&& fcovariance_block_kernel = fcovariance_matrix_storage.get_fcovariance_kernel_for_block( n_block, n_end );
			//fcovariance_block_kernel.clear();

			const auto& sequence_blocks = blocks[n_block];
			const auto& sequence_indices_for_blocks = block_accounting[n_block];

			for( std::size_t block_index=0; block_index < block_accounting[n_block].size(); ++block_index )
			{
				auto block_weights& = *(m_block_weights)[n_block][block_index];

				for( const auto si: sequence_indices_for_blocks[block_index] ) // loop over sequence indices for this block
				{
					wcache[std::size_t( r_states[si] )] += weights[si];
				}

				fcovariance_block_kernel.accumulate( sequence_blocks[block_index], wcache );
			}
		}
*/
	}

	void cache_column_frequencies() const
	{
		m_frequencies = std::make_shared<frequencies_t>( this->n_loci() );

	    const std::size_t n_loci = this->n_loci(); // number of columns in the alignment
	    const std::size_t n_loci_per_block = apegrunt::StateBlock_size;
	    const std::size_t last_block_size = get_last_block_size(n_loci);
	    const std::size_t last_block_index = get_last_block_index(n_loci);
	    const std::size_t number_of_blocks = apegrunt::get_number_of_blocks(this->n_loci());

	    const auto& block_accounting = *(this->get_block_accounting());
		const auto& blocks = *(this->get_block_storage());

		auto& frequencies = *m_frequencies;

		for( std::size_t n_block=0; n_block < number_of_blocks; ++n_block )
		{
			const auto n_end = ( n_block == last_block_index ? last_block_size : n_loci_per_block );
			auto&& column_frequency_accumulator = apegrunt::Column_frequency_accumulator<typename frequency_t::value_type, number_of_states<state_t>::value, apegrunt::StateBlock_size>( frequencies[n_block*n_loci_per_block].data(), n_end );
			column_frequency_accumulator.setZero();
			const auto& sequence_blocks = blocks[n_block];
			const auto& indices = block_accounting[n_block];

			for( std::size_t block_index=0; block_index < indices.size(); ++block_index )
			{
				column_frequency_accumulator.accumulate( sequence_blocks[block_index], indices[block_index] );
			}
		}
		//std::cout << "}"; std::cout.flush();
	}

	void cache_column_w_frequencies() const
	{
		m_w_frequencies = std::make_shared<w_frequencies_t>( this->n_loci() );

	    const std::size_t n_loci = this->n_loci(); // number of columns in the alignment
	    const std::size_t n_loci_per_block = apegrunt::StateBlock_size;
	    const std::size_t last_block_size = get_last_block_size(n_loci);
	    const std::size_t last_block_index = get_last_block_index(n_loci);
	    const std::size_t number_of_blocks = apegrunt::get_number_of_blocks(this->n_loci());

	    const auto& block_accounting = *(this->get_block_accounting());
		const auto& blocks = *(this->get_block_storage());

		std::vector<typename StateVector<state_t>::weight_type > weights; weights.reserve(this->size());
		for( const auto& sequence: m_rows ) { weights.push_back(sequence->weight()); }

		auto& frequencies = *m_w_frequencies;

		for( std::size_t n_block=0; n_block < number_of_blocks; ++n_block )
		{
			const auto n_end = ( n_block == last_block_index ? last_block_size : n_loci_per_block );
			auto&& column_frequency_accumulator = apegrunt::Column_frequency_accumulator<typename w_frequency_t::value_type, number_of_states<state_t>::value, apegrunt::StateBlock_size>( frequencies[n_block*n_loci_per_block].data(), n_end );
			column_frequency_accumulator.setZero();
			const auto& sequence_blocks = blocks[n_block];
			const auto& indices = block_accounting[n_block];

			for( std::size_t block_index=0; block_index < indices.size(); ++block_index )
			{
				column_frequency_accumulator.accumulate( sequence_blocks[block_index], indices[block_index], weights.data() );
			}
		}
	}
/*
	// the triangular sample-sample Hamming distance matrix
	void cache_distance_matrix()
	{
		const std::size_t n_seqs = this->size();
		m_distance_matrix = std::make_shared<distance_matrix_t>( n_seqs*(n_seqs-1)/2 );

		auto& dmat = *m_distance_matrix;

		// compute matrix of pairwise sequence identities
		for( std::size_t i = 0; i < n_seqs; ++i )
		{
			const auto seq_i = (*this)[i];
			for( std::size_t j = 0; j < i; ++j )
			{
				const std::size_t n_elem = std::min( seq_i->size(), (*this)[j]->size() );
				dmat[i*(i-1)/2+j] = n_elem - ( *seq_i && *(*this)[j] ); // bitwise-AND counts the number of identical positions
			}
		}
	}
*/
// /*
	// the triangular sample-sample Hamming distance matrix
	void cache_distance_matrix() const
	{
		const std::size_t n_seqs = this->size();
		m_distance_matrix = std::make_shared<distance_matrix_t>( n_seqs*(n_seqs-1)/2, 0 );

		auto& dmat = *m_distance_matrix;

		const std::size_t n_loci = this->n_loci(); // number of columns in the alignment
		const std::size_t n_loci_per_block = apegrunt::StateBlock_size;
		const std::size_t last_colblock_size = get_last_block_size(n_loci);
		const std::size_t last_colblockidx = get_last_block_index(n_loci);
	    const std::size_t number_of_blocks = apegrunt::get_number_of_blocks(this->n_loci());

		const auto& block_accounting = *(this->get_block_accounting());
		const auto& blocks = *(this->get_block_storage());

		for( std::size_t colblockidx=0; colblockidx < number_of_blocks; ++colblockidx )
		{
			const auto& sequence_blocks = blocks[colblockidx];
			const auto& indices = block_accounting[colblockidx];

			for( std::size_t bidx_i = 0; bidx_i < sequence_blocks.size(); ++bidx_i )
			{
				const auto block_i = sequence_blocks[bidx_i];
				const auto& block_i_seqs = indices[bidx_i];
				for( std::size_t bidx_j = 0; bidx_j < bidx_i; ++bidx_j )
				{
					const auto ndiff = n_loci_per_block - count_identical( block_i, sequence_blocks[bidx_j] );
					const auto& block_j_seqs = indices[bidx_j];
					for( const auto i : block_i_seqs )
					{
						for( const auto j : block_j_seqs )
						{
							dmat[ i < j ? j*(j-1)/2+i : i*(i-1)/2+j ] += ndiff; // bad inner-loop branched memory access
						}
					}
				}
			}
		}
	}

	void cache_statecount_and_statepresence_blocks() const
	{
		//std::cout << "apegrunt: cache statecount and statepresence accounting" << std::endl;
		const auto& blocks = *(this->get_block_storage());

		m_statecount_blocks = std::make_shared< statecount_block_storage_t >(); m_statecount_blocks->reserve( blocks.size() );
		auto& statecount_blocks = *m_statecount_blocks;

		m_statepresence_blocks = std::make_shared< statepresence_block_storage_t >(); m_statepresence_blocks->reserve( blocks.size() );
		auto& statepresence_blocks = *m_statepresence_blocks;

		// some temp definitions
		statecount_block_t ones( statecount_t(1) );

		// create statepresence masks
		for( const auto& block_column: blocks )
		{
			statecount_block_t statepresence_mask( statecount_t(0) );
			//std::cout << " at onset: \"" << statepresence_mask << "\"" << std::endl;
			for( const auto block: block_column )
			{
				statepresence_mask = statepresence_mask | ( ones << block );
			}
			statepresence_blocks.push_back( statepresence_mask );
			statecount_blocks.push_back( apegrunt::popcnt_per_element(statepresence_mask) );
			//std::cout << " at end:   \"" << statecount_blocks.back() << "\"" << std::endl;
		}

		// zero out padding elements of the last column block
		for( std::size_t i = apegrunt::get_last_block_size(this->n_loci()); i < N; ++i ) { statepresence_blocks.back()[i] = 0; statecount_blocks.back()[i] = 0; }
		//std::cout << "m_statecount_blocks->size()=" << m_statecount_blocks->size() << " m_statepresence_blocks->size()=" << m_statepresence_blocks->size() << std::endl;
	}

// */
	// development code
	void cache_column_block_distance() const
	{
		std::vector< bool > ident; ident.resize( this->size() );
		std::vector< block_type > refblocks; refblocks.resize( this->size() );

		const auto& block_accounting = *(this->get_block_accounting());
		const auto& blocks = *(this->get_block_storage());

		for( std::size_t coli = 0; coli < blocks.size(); ++coli )
		{
			for( std::size_t i = 0; i < blocks[coli].size(); ++i )
			{
				for( auto seq : block_accounting[coli][i] )
				{
					refblocks[seq] = blocks[coli][i];
				}
			}
			for( std::size_t colj = 0; colj < coli; ++ colj )
			{
				for( std::size_t i = 0; i < blocks[colj].size(); ++i )
				{
					for( auto seq : block_accounting[colj][i] )
					{
						ident[seq] = ( refblocks[seq] == blocks[colj][i] );
					}
				}
				std::size_t sum = 0; for( auto val: ident ) { sum += ( val ? 1 : 0 ); }
				std::cout << "coli=" << coli << " colj=" << colj << " ident=" << sum << "\n";
			}
		}
	}
/*
	void loci_distance_matrix() const
	{
		const std::size_t n_loci = this->n_loci();

		if( n_loci > 10000 ) { return; }

		auto m_loci_distance_matrix = std::make_shared<distance_matrix_t>( n_loci*(n_loci-1)/2 );

		auto& dmat = *m_loci_distance_matrix;

		for( auto seqptr : *this )
		{
			const auto *seq = seqptr.get(); // get raw pointer
			for( std::size_t i = 0; i < n_loci; ++i )
			{
				for( std::size_t j = 0; j < i; ++j )
				{
					dmat[i*(i-1)/2+j] += ( (*seq)[i] != (*seq)[j] );
				}
			}
		}

		std::vector<std::size_t> distance_profile(this->size(), 0);

		for( std::size_t i = 0; i < n_loci; ++i )
		{
			for( std::size_t j = 0; j < i; ++j )
			{
				distance_profile[ dmat[i*(i-1)/2+j] ] += 1;
			}
		}

		for( std::size_t i = 0; i < distance_profile.size(); ++i )
		{
			std::cout << "rank=" << i << " count=" << distance_profile[i] << std::endl;
		}

	}
*/
	/// Parser interface

	statevector_t* get_new_sequence( const std::string& name )
	{
		//std::cout << "Add new sequence named \"" << name << "\"" << std::endl;

		std::size_t reserve = 0;
		if( m_rows.size() > 0 ) {
			reserve = m_rows.back()->size() / N + ( m_rows.back()->size() % N == 0 ? 0 : 1 );
		}

		m_rows.emplace_back( make_StateVector_ptr<statevector_t>( m_block_storage, name, reserve ) );
		return static_cast<statevector_t*>(m_rows.back().get());
	}

	// Allow parser access to private members
	ALIGNMENT_PARSER_GRAMMAR_FRIENDS(my_type)

	friend class Alignment_factory<my_type>;

	/// Generator interface

	// Allow generator access to private members
	ALIGNMENT_GENERATOR_GRAMMAR_FRIENDS(my_type)

	/// boost.serialization interface.
	friend class boost::serialization::access;
	template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_NVP(boost::serialization::base_object< base_type >(*this));
        ar & BOOST_SERIALIZATION_NVP(m_rows);
    }

	/// Helper class for memory management thru std::shared_ptr
	class deleter
	{
	public:
		void operator()( my_type* p )
		{
			delete p;
		}
	};
	friend class deleter;

};

} // namespace apegrunt

#endif // APEGRUNT_ALIGNMENT_IMPL_BLOCK_COMPRESSED_STORAGE_HPP

