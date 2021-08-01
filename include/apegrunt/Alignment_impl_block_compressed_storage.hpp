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
#include <vector>
#include <utility> // for std::pair

#include <boost/range/adaptor/indexed.hpp>
#include <boost/range/irange.hpp>

#include "Apegrunt_options.h"

#include "Alignment_forward.h"
#include "Alignment_impl_base.hpp"
#include "StateVector_forward.h"

#include "Alignment_parser_forward.h"
#include "Alignment_generator_forward.h"

#include "Alignment_iterator.h"
//#include "Alignment_iterator_impl_block_compressed_storage_forward.h"
#include "Alignment_iterator_impl_block_compressed_storage.hpp"

#include "State_block.hpp"

#include "accumulators/distribution.hpp"
#include "accumulators/distribution_std.hpp"
#include "accumulators/distribution_bincount.hpp"
#include "accumulators/distribution_cumulative.hpp"
#include "accumulators/distribution_ordered.hpp"
#include "accumulators/distribution_generator_csv.hpp"
#include "Apegrunt_IO_misc.hpp"

#include "Apegrunt_utility.hpp"
#include "misc/TernarySearchTree.hpp"
#include "misc/Math.hpp"
#include "misc/Stopwatch.hpp"

namespace apegrunt {

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

	~Column_frequency_accumulator() = default;

	Column_frequency_accumulator( const my_type& other )
	: m_data(other.m_data), m_extent(other.m_extent)
	{ }

	my_type& operator=( const my_type& other )
	{
		m_data = other.m_data;
		m_extent = other.m_extent;
	}

	template< typename StateT, typename IndexContainerT >
	inline void accumulate( const apegrunt::State_block<StateT,BlockSize>& stateblock, const IndexContainerT& sequence_indices, frequency_t* const weights=nullptr )
	{
		const frequency_t weight = weights
			? apegrunt::gather( sequence_indices, weights )
			: frequency_t( sequence_indices.size() );
		for( std::size_t i=0; i < m_extent; ++i )
		{
			const auto state = std::size_t( stateblock[i] );
			*(m_data+i*N+state) += weight;
		}
	}

	inline void setZero()
	{
		std::memset(m_data,0,m_extent*sizeof(frequency_t));
	}

private:

	frequency_t* const m_data;
	const std::size_t m_extent;
};

namespace acc = boost::accumulators;

namespace apegrunt {

// forward declarations


// Helper for selecting TST-key blocking factor based on state type.
template< typename StateT > struct TST_block_size { static const std::size_t value=std::max(1,std::min(apegrunt::StateBlock_size/2,8)); };
template<> struct TST_block_size<apegrunt::nucleic_acid_state_t> { static const std::size_t value=std::max(1,std::min(apegrunt::StateBlock_size/2,16)); };
template<> struct TST_block_size<char> { static const std::size_t value=std::max(1,std::min(apegrunt::StateBlock_size/2,16)); };

template< typename StateT, typename InternalState=StateT, std::size_t Cache=0 >
class Block_adder
{
public:
	using state_t = StateT;
	using internal_state_t = InternalState;
	using my_type = Block_adder<state_t,internal_state_t>;

	using block_index_t = typename Alignment<state_t>::block_index_t;
	using block_accounting_t = typename Alignment<state_t>::block_accounting_t;
	using block_accounting_ptr = typename Alignment<state_t>::block_accounting_ptr;
	using block_accounting_container_t = typename Alignment<state_t>::block_accounting_container_t;

	using block_storage_t = typename Alignment<state_t>::block_storage_t;
	using block_storage_ptr = typename Alignment<state_t>::block_storage_ptr;

	using block_type = typename Alignment<state_t>::block_type;
	using internal_block_type = typename Alignment<internal_state_t>::block_type;

	using block_finder_t = apegrunt::TernarySearchTree<internal_block_type,block_accounting_container_t,TST_block_size<internal_state_t>::value>;
	//using block_finder_t = apegrunt::TernarySearchTree<block_type,block_index_t>;
	using block_finder_storage_t = std::vector< block_finder_t >;

	Block_adder() = delete;

	Block_adder( block_storage_ptr block_storage, block_accounting_ptr block_accounting )
	: m_block_storage_ptr( block_storage ),
	  m_block_accounting_ptr( block_accounting ),
	  m_extracted(false)
	{
	}

	Block_adder( my_type&& other )
	: m_block_storage_ptr( other.m_block_storage_ptr ),
	  m_block_accounting_ptr( other.m_block_accounting_ptr ),
	  m_block_finders( std::move( other.m_block_finders ) ),
	  m_extracted(false)
	{
	}

	inline std::size_t size() const { return m_block_storage_ptr->size(); }

	// cached version
	template< std::size_t W = Cache, std::enable_if_t<(W!=0),bool> = true >
	inline void reserve( std::size_t size )
	{
		if( size > m_block_buffer.size() )
		{
			const std::lock_guard<std::mutex> lock(m_buffer_resize_mtx);
			m_block_buffer.resize(size);
			//m_block_finders.resize(size);
		}
	}


	// write-thru version
	template< std::size_t W = Cache, std::enable_if_t<(W==0),bool> = false >
	inline void reserve( std::size_t size )
	{
		if( size > m_block_finders.size() )
		{
			const std::lock_guard<std::mutex> lock(m_finder_resize_mtx);
			m_block_finders.resize(size);
		}
	}

	// cached version -- this should be the preferred way most of the time; to keep the search tree
	// hot in cpu cache during insert, we'll perform insertion bursts. This requires that we temporarily
	// store blocks in a block buffer that we then flush in one go into the actual search tree.
	template< typename Block, std::size_t W = Cache, std::enable_if_t<(W!=0),bool> = true >
	inline void insert( Block block, std::size_t block_col, block_index_t source_index )
	{
		// we won't ever actually hit these if reserve was properly called
		//if( m_block_finders.size() <= block_col ) { const std::lock_guard<std::mutex> lock(m_finder_resize_mtx); m_block_finders.resize(block_col+1); }
		if( m_block_buffer.size() <= block_col ) { const std::lock_guard<std::mutex> lock(m_buffer_resize_mtx); m_block_buffer.resize(block_col+1); }

		m_block_buffer[block_col].emplace_back( block, source_index );
		if( m_block_buffer[block_col].size() == Cache )	{ this->flush(block_col); }
	}

	// write-thru version
	template< typename Block, std::size_t W = Cache, std::enable_if_t<(W==0),bool> = false >
	inline void insert( Block block, std::size_t block_col, block_index_t source_index )
	{
		if( m_block_finders.size() <= block_col ) { const std::lock_guard<std::mutex> lock(m_finder_resize_mtx); m_block_finders.resize(block_col+1); }
		m_block_finders[block_col].insert( block )->value.push_back(source_index);
	}

	inline block_type get( std::size_t block_col, block_index_t block_index ) const
	{
		return (m_block_storage_ptr->size() > block_col) ? ( (*m_block_storage_ptr)[block_col].size() > block_index ? (*m_block_storage_ptr)[block_col][block_index] : block_type() ) : block_type();
	}

	inline block_storage_ptr get_block_storage() { if(!m_extracted) { this->extract(); } return m_block_storage_ptr; }
	inline block_accounting_ptr get_block_accounting() { if(!m_extracted) { this->extract(); } return m_block_accounting_ptr; }

	// cached version -- flush all
	template< std::size_t W = Cache, std::enable_if_t<(W!=0),bool> = true >
	inline void flush() // flush all
	{
		if( m_block_finders.size() < m_block_buffer.size() ) { const std::lock_guard<std::mutex> lock(m_finder_resize_mtx); m_block_finders.resize(m_block_buffer.size()); }
		for( auto dest_and_src: apegrunt::zip_range( m_block_finders, m_block_buffer ) )
		{
			using boost::get;
			auto& tree = get<0>(dest_and_src);
			for( const auto& pair: get<1>(dest_and_src) )
			{
				tree.insert(pair.first)->value.push_back(pair.second);
			}
			get<1>(dest_and_src).clear();
		}
	}

	// cached version -- flush single
	template< std::size_t W = Cache, std::enable_if_t<(W!=0),bool> = true >
	inline void flush( std::size_t col ) // flush only col
	{
		if( m_block_finders.size() < m_block_buffer.size() ) { const std::lock_guard<std::mutex> lock(m_finder_resize_mtx); m_block_finders.resize(m_block_buffer.size()); }

		auto& tree = m_block_finders[col];
		for( const auto& pair: m_block_buffer[col] )
		{
			tree.insert(pair.first)->value.push_back(pair.second);
		}
		m_block_buffer[col].clear();
	}

	// write-thru version
	template< std::size_t W = Cache, std::enable_if_t<(W==0),bool> = false >
	inline void flush() { }

	void extract()
	{
		if( m_block_finders.empty() && m_block_buffer.empty() ) { m_extracted = true; return; }

		this->flush();

		auto& blocks = *m_block_storage_ptr;
		auto& accounting = *m_block_accounting_ptr;

		blocks.reserve( m_block_finders.size() );
		accounting.reserve( m_block_finders.size() );

		m_statistics.ntrees = m_block_finders.size();

		for( auto& tree: m_block_finders )
		{
			m_statistics.nkeys += tree.nkeys();
			m_statistics.nunique += tree.nunique_keys();
			m_statistics.nnodes += tree.nnodes();
			m_statistics.minodes = std::min( m_statistics.minodes, tree.nnodes() );
			m_statistics.manodes = std::max( m_statistics.manodes, tree.nnodes() );
			m_statistics.bytesize += tree.bytesize();

			accounting.emplace_back(); blocks.emplace_back();
			auto& a = accounting.back(); auto& b = blocks.back();
			auto fb( [&a, &b, this] (auto node) mutable { /*m_statistics.ipbn_distrib(node->value.size()); m_statistics.ipbk_distrib(node->value.storage().size());*/ a.emplace_back( std::move(node->value) ); b.emplace_back(std::move(node->key)); } );
			tree.for_each( fb );
		}
		m_extracted = true;
		m_block_finders.clear(); // clean-up
	}

	inline void statistics( std::ostream *out=nullptr ) const
	{
		//static const std::size_t N=8;
		if( !out ) { return; }

		if( !m_extracted ) { /*this->extract();*/ return; }

		{
			if( m_statistics.ntrees == 0 )
			{
				*out << "apegrunt: B|3: parse tree is empty\n"; // out->flush();
			}
			else
			{
				*out << "apegrunt: B|3:"
						<< " parse tree size = " << apegrunt::memory_string( m_statistics.bytesize )
						<< " (" << double(m_statistics.bytesize)/double(m_statistics.nkeys)*8 << " bits/stored object)"
						<< " | nodes = " << m_statistics.nnodes << " [min/mean/max per tree = "
						<< m_statistics.minodes << "/" << std::size_t(double(m_statistics.nnodes)/double(m_statistics.ntrees)) << "/" << m_statistics.manodes << "]"
						//<< " global<" << N << ">(" << this->nblocksglobal<N>() << ")"
/*
						<< " index(mean=" << boost::accumulators::distribution_mean(m_statistics.ipbn_distrib)
						<< ",std=" << boost::accumulators::distribution_std(m_statistics.ipbn_distrib)
						<< ",bins=" << boost::accumulators::distribution_bincount(m_statistics.ipbn_distrib)
						<< ")"
						<< " bm(mean=" << boost::accumulators::distribution_mean(m_statistics.ipbk_distrib)
						<< ",std=" << boost::accumulators::distribution_std(m_statistics.ipbk_distrib)
						<< ",bins=" << boost::accumulators::distribution_bincount(m_statistics.ipbk_distrib)
						<< ")"
*/
						<< "\n"
				;
/*
				auto csv_file = apegrunt::get_unique_ofstream( "bm_index_distribution.csv" );
				*csv_file->stream() << apegrunt::accumulators::csv(acc::distribution(m_statistics.ipbn_distrib));
				//*csv_file->stream() << apegrunt::accumulators::csv(apegrunt::accumulators::distribution_cumulative(m_statistics.ipb_distrib));
				auto csv_file2 = apegrunt::get_unique_ofstream( "bm_block_distribution.csv" );
				*csv_file2->stream() << apegrunt::accumulators::csv(acc::distribution(m_statistics.ipbk_distrib));
*/
/*
				auto& accounting = *m_block_accounting_ptr;
				for( auto& col: accounting ) { for( auto& index: col ) {
						for( auto& e: index.storage() )
						{
							std::cout << "payload=" << e.m_payload << " pos=" << e.m_pos << " bitfield=" << e.m_bitfield << std::endl;
						}
				} }
*/
			}
		}
	}

private:
	block_storage_ptr m_block_storage_ptr;
	block_accounting_ptr m_block_accounting_ptr;
	block_finder_storage_t m_block_finders;
	std::vector< std::vector< std::pair<internal_block_type,std::size_t> > > m_block_buffer;
	bool m_extracted;
	std::mutex m_buffer_resize_mtx;
	std::mutex m_finder_resize_mtx;

	struct statistics
	{
		statistics()
		: ntrees(0), nkeys(0), nunique(0), nnodes(0), minodes(std::numeric_limits<std::size_t>::max() ), manodes(0), bytesize(0) { }
//		  ipbn_distrib(acc::tag::distribution::binwidth=1), ipbk_distrib(acc::tag::distribution::binwidth=1) { }
		std::size_t ntrees, nkeys, nunique, nnodes, minodes, manodes, bytesize;
//		acc::accumulator_set<double, acc::stats<acc::tag::std(acc::from_distribution),acc::tag::distribution_bincount> > ipbn_distrib;
//		acc::accumulator_set<double, acc::stats<acc::tag::std(acc::from_distribution),acc::tag::distribution_bincount> > ipbk_distrib;
	} m_statistics;
};
template< std::size_t BlockSize >
struct Subset_block_structure_builder
{
	static constexpr std::size_t N=BlockSize;
	using index_storage_t = std::vector< std::vector<std::size_t> >;
	using const_iterator = typename index_storage_t::const_iterator;

	Subset_block_structure_builder() = delete;
	explicit Subset_block_structure_builder( std::size_t blockn ) : m_nfill(0), m_blockn(blockn) { m_src_posns.emplace_back(); }

	inline bool add( std::size_t pos )
	{
		if( m_nfill == N ) { return false; } // we're full

		if( !m_src_posns.back().empty() && apegrunt::get_block_index( m_src_posns.back().back() ) != apegrunt::get_block_index(pos) ) { m_src_posns.emplace_back(); }
		m_src_posns.back().emplace_back(pos); ++m_nfill;

		return true;
	}

	inline std::size_t size() const { return m_nfill; }
	inline std::size_t nblocks() const { return m_src_posns.size(); }
	inline std::size_t block_index() const { return m_blockn; }

	inline const_iterator begin() const { return m_src_posns.begin(); }
	inline const_iterator end() const { return m_src_posns.end(); }

	std::size_t m_nfill;
	const std::size_t m_blockn;

	index_storage_t m_src_posns;
};

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
	using block_accounting_container_t = typename base_type::block_accounting_container_t;
	using block_accounting_t = typename base_type::block_accounting_t;
	using block_accounting_ptr = typename base_type::block_accounting_ptr;

	using block_container_t = typename base_type::block_container_t;
	using block_storage_t = typename base_type::block_storage_t;
	using block_storage_ptr = typename base_type::block_storage_ptr;

	using block_adder_t = Block_adder<state_t,char,32>;
	using block_adder_ptr = std::shared_ptr< block_adder_t >;

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

	Alignment_impl_block_compressed_storage()
	: base_type(),
	  m_rows(),
	  m_block_storage( std::make_shared<block_storage_t>() ),
	  m_block_accounting( std::make_shared<block_accounting_t>() ),
	  m_block_adder( std::make_shared<block_adder_t>( m_block_storage, m_block_accounting ) )
	{
	}

	~Alignment_impl_block_compressed_storage() = default;

	Alignment_impl_block_compressed_storage( const my_type& other )
		: base_type( other.id_string() ),
		  m_rows( other.m_rows ),
		  m_block_storage( other.m_block_storage ),
		  m_block_accounting( other.m_block_accounting ),
		  m_block_adder( other.m_block_adder )
	{
		//std::cout << "Alignment_impl_block_compressed_storage( const my_type& other )" << std::endl;
	}

	Alignment_impl_block_compressed_storage( my_type&& other ) noexcept
		: base_type( other.id_string() ),
		  m_rows( std::move(other.m_rows) ),
		  m_block_storage( std::move(other.m_block_storage) ),
		  m_block_accounting( std::move(other.m_block_accounting) ),
		  m_block_adder( std::move(other.m_block_adder) )
	{
		//std::cout << "Alignment_impl_block_compressed_storage( my_type&& other )" << std::endl;
	}

	Alignment_ptr<state_t> clone() const
	{
		return make_Alignment_ptr( *this );
	}

    inline const_iterator cbegin() const { return const_iterator( std::make_unique<const_iterator_impl>( m_rows.cbegin() ) ); }
    inline const_iterator cend() const { return const_iterator( std::make_unique<const_iterator_impl>( m_rows.cend() ) ); }

    inline const_iterator begin() const { return this->cbegin(); }
    inline const_iterator end() const { return this->cend(); }

    inline iterator begin() { return iterator( std::make_unique<iterator_impl>( m_rows.begin() ) ); }
    inline iterator end() { return iterator( std::make_unique<iterator_impl>( m_rows.end() ) ); }

    inline value_type operator[]( std::size_t index ) const { return m_rows[index]; }

    iterator erase( iterator first, iterator last )
    {
    	using std::cbegin; using std::cend;
    	const auto beg = std::find( cbegin(m_rows), cend(m_rows), *first );
    	const auto end = std::find( beg, cend(m_rows), *last );
    	return iterator( std::make_unique<iterator_impl>( m_rows.erase(beg,end) ) );
    }

    inline std::size_t size() const { return m_rows.size(); }
    inline std::size_t n_loci() const { return m_nloci; }

    inline const std::type_info& type() const { return typeid(my_type); }

    inline frequencies_ptr frequencies() const
    {
    	m_cache_frequencies_mutex.lock();
    	if( !m_frequencies ) { this->cache_column_frequencies(); }
    	m_cache_frequencies_mutex.unlock();
    	return m_frequencies;
    }

    inline w_frequencies_ptr w_frequencies() const
    {
    	m_cache_w_frequencies_mutex.lock();
    	if( !m_w_frequencies ) { this->cache_column_w_frequencies(); }
    	m_cache_w_frequencies_mutex.unlock();
    	return m_w_frequencies;
    }

    inline Alignment_subscript_proxy< StateVector_ptr<state_t> > subscript_proxy() const { return Alignment_subscript_proxy< StateVector_ptr<state_t> >( &m_rows ); }

	block_accounting_ptr get_block_accounting() const
	{
		return m_block_accounting;
	}

	block_storage_ptr get_block_storage() const
	{
		return m_block_storage;
	}

	statecount_block_storage_ptr get_statecount_blocks() const
	{
		m_cache_statecount_and_statepresence_blocks_mutex.lock();
		if( !m_statecount_blocks ) { this->cache_statecount_and_statepresence_blocks(); }
		m_cache_statecount_and_statepresence_blocks_mutex.unlock();

		return m_statecount_blocks;
	}

	statepresence_block_storage_ptr get_statepresence_blocks() const
	{
		m_cache_statecount_and_statepresence_blocks_mutex.lock();
		if( !m_statepresence_blocks ) { this->cache_statecount_and_statepresence_blocks(); }
		m_cache_statecount_and_statepresence_blocks_mutex.unlock();

		return m_statepresence_blocks;
	}

	statepresence_block_storage_ptr get_statepresence_blocks_wo_gaps() const
	{
		m_cache_statecount_and_statepresence_blocks_mutex.lock(); // shares mutex with it's buddies
		if( !m_statepresence_blocks_wo_gaps ) { this->cache_statecount_and_statepresence_blocks(); }
		m_cache_statecount_and_statepresence_blocks_mutex.unlock();

		return m_statepresence_blocks_wo_gaps;
	}

	statepresence_block_storage_ptr get_gappresence_blocks() const
	{
		m_cache_statecount_and_statepresence_blocks_mutex.lock(); // shares mutex with it's buddies
		if( !m_gappresence_blocks ) { this->cache_statecount_and_statepresence_blocks(); }
		m_cache_statecount_and_statepresence_blocks_mutex.unlock();

		return m_gappresence_blocks;
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

		*out << "apegrunt: statistics for alignment \"" << this->id_string() << "\":\n";

		const std::size_t n_loci = this->n_loci();
		const std::size_t n_seqs = m_rows.size();
		const std::size_t n_block_columns = m_block_storage->size();

		std::size_t n_blocks = 0;
		std::size_t min_blocks_per_column = m_block_storage->empty() ? 0 : m_block_storage->front().size();
		std::size_t max_blocks_per_column = 0;
		std::vector<std::size_t> col_size; col_size.reserve(m_block_storage->size());

		for( const auto& block_col: *m_block_storage )
		{
			const std::size_t nb = block_col.size();
			col_size.push_back( nb );

			min_blocks_per_column = std::min(min_blocks_per_column,nb);
			max_blocks_per_column = std::max(max_blocks_per_column,nb);
			n_blocks += nb;
		}

		auto nth = col_size.begin()+col_size.size()/2 - (col_size.size() & 0x1); std::nth_element( col_size.begin(), nth, col_size.end() );
		const std::size_t n_col_median = *nth;
		const std::size_t n_col_mean = n_blocks / n_block_columns;

		std::size_t compressed_mem(0), compressed_index_mem(0), compressed_blocks_mem(0);

		for( auto&& as_bs: apegrunt::zip_range(*m_block_accounting, *m_block_storage) )
		{
			using boost::get;
			auto& as = get<0>(as_bs);
			auto& bs = get<1>(as_bs);

			compressed_blocks_mem += bs.size()*sizeof(block_type);

			for( const auto& a: as )
			{
				compressed_index_mem += bytesize(a);
			}
		}

		compressed_mem += compressed_index_mem+compressed_blocks_mem;

		// Get state statistics
		using std::begin; using std::end;
		std::vector<std::size_t> state_stats(number_of_states<state_t>::value, 0);

		//std::size_t b(0);
		for( const auto& block: *(this->get_statecount_blocks()) )
		{
			//for( auto count: block ) { count != 0 && ++state_stats[ count-1 ]; }
			for( std::size_t i=0; i < N; ++i ) { block[i] != 0 && ++state_stats[ block[i]-1 ]; } // need to guard for negative indices, as the last block may contain zero-valued entries
			//for( std::size_t i=0; i < N; ++i ) { if( block[i] == 1 ) { std::cout << "block=" << b << " pos=" << b*StateBlock_size+i << " i=" << i << std::endl; } }
			//++b;
		}

		std::size_t npoly = 0;
		for( std::size_t i=1; i < state_stats.size(); ++i ) { npoly += state_stats[i]; }
		// end Get state statistics

		*out << "apegrunt: alignment has " << n_seqs << " samples and " << n_loci << " positions\n"; // (" << onucs << " outcomes in total)\n";
		*out << "apegrunt: state distribution [#states:#positions] =";
		std::size_t statesum(0);
		for( std::size_t i=0; i < state_stats.size(); ++i )
		{
			if( state_stats[i] != 0 )
			{
				 statesum += state_stats[i];
				*out << " " << i+1 << ":" << state_stats[i];
			}
		}
		*out << " | total:" << statesum << " | polymorphic:" << npoly << "\n";
		*out << "apegrunt: #stateblocks = " << n_blocks << " (min/mean/median/max per column = " << min_blocks_per_column << "/" << n_col_mean << "/" << n_col_median << "/" << max_blocks_per_column << ")\n";
		out->flush();

		const std::size_t nvariables = n_loci*n_seqs; // # of uncompressed nucleotides
		const std::size_t uncompressed_mem = nvariables*sizeof(state_t);

		*out << "apegrunt: uncompressed size = " << apegrunt::memory_string(uncompressed_mem) << "\n";
		*out << "apegrunt: compressed size = " << apegrunt::memory_string(compressed_mem) << " (" << apegrunt::memory_string(compressed_blocks_mem) << " for " << n_blocks << " blocks and " << apegrunt::memory_string(compressed_index_mem) << " for index)\n";
		*out << "apegrunt: compression ratio = " << double(uncompressed_mem)/double(compressed_mem)
			<< " (" << double(compressed_mem)/double(nvariables)*8 << " bits/variable)\n";
		//out->flush();

		m_block_adder->statistics( out ); out->flush();
	}

private:
	using iterator_impl = apegrunt::iterator::Alignment_iterator_impl_block_compressed_storage< StateVector_ptr<state_t> >;
	using const_iterator_impl = apegrunt::iterator::Alignment_const_iterator_impl_block_compressed_storage< StateVector_ptr<state_t> >;

	// The data row proxies and the actual data matrix should in general not
	// be modified during the lifetime of a const Alignment..
	std::vector< StateVector_ptr<state_t> > m_rows;
	block_storage_ptr m_block_storage;
	block_accounting_ptr m_block_accounting;
	block_adder_ptr m_block_adder; // should always be declared *after* m_block_storage

	// ..but all caches may be modified and are therefore declared mutable
	mutable frequencies_ptr m_frequencies;
	mutable w_frequencies_ptr m_w_frequencies;
	mutable statecount_block_storage_ptr m_statecount_blocks;
	mutable statecount_block_storage_ptr m_statepresence_blocks;
	mutable statecount_block_storage_ptr m_statepresence_blocks_wo_gaps;
	mutable statecount_block_storage_ptr m_gappresence_blocks;

	mutable std::mutex m_cache_frequencies_mutex;
	mutable std::mutex m_cache_w_frequencies_mutex;
	mutable std::mutex m_cache_statecount_and_statepresence_blocks_mutex;

	std::size_t m_nloci = 0;

	void cache_column_frequencies() const
	{
		m_frequencies = std::make_shared<frequencies_t>( this->n_loci() );
		auto& frequencies = *m_frequencies;

		const auto block_accounting_ptr = this->get_block_accounting(); // keep shared_ptr alive
		const auto blocks_ptr = this->get_block_storage(); // keep shared_ptr alive
	    const auto& block_accounting = *block_accounting_ptr;
		const auto& blocks = *blocks_ptr;

		const std::size_t n_loci = this->n_loci(); // number of columns in the alignment
	    const std::size_t n_loci_per_block = apegrunt::StateBlock_size;
	    const std::size_t last_block_size = apegrunt::get_last_block_size(n_loci);
	    //const std::size_t last_block_index = get_last_block_index(n_loci);
	    //const std::size_t number_of_blocks = apegrunt::get_number_of_blocks(this->n_loci());
	    const std::size_t number_of_blocks = blocks.size();
	    const std::size_t last_block_index = number_of_blocks != 0 ? number_of_blocks-1 : 0;

		for( std::size_t n_block=0; n_block < number_of_blocks; ++n_block )
		{
			const auto n_end = ( n_block == last_block_index ? last_block_size : n_loci_per_block );
			using fq_acc_type = apegrunt::Column_frequency_accumulator<typename frequency_t::value_type, number_of_states<state_t>::value, apegrunt::StateBlock_size>;
			fq_acc_type column_frequency_accumulator( frequencies[n_block*n_loci_per_block].data(), n_end );
			column_frequency_accumulator.setZero();

			for( auto&& bs_as: apegrunt::zip_range(blocks[n_block], block_accounting[n_block]) )
			{
				using boost::get;
				column_frequency_accumulator.accumulate( get<0>(bs_as), get<1>(bs_as) );
			}
		}
	}

	void cache_column_w_frequencies() const
	{
		m_w_frequencies = std::make_shared<w_frequencies_t>( this->n_loci() );

	    const std::size_t n_loci = this->n_loci(); // number of columns in the alignment
	    const std::size_t n_loci_per_block = apegrunt::StateBlock_size;
	    const std::size_t last_block_size = get_last_block_size(n_loci);
	    const std::size_t last_block_index = get_last_block_index(n_loci);
	    const std::size_t number_of_blocks = apegrunt::get_number_of_blocks(this->n_loci());

		const auto block_accounting_ptr = this->get_block_accounting(); // keep shared_ptr alive
		const auto blocks_ptr = this->get_block_storage(); // keep shared_ptr alive
	    const auto& block_accounting = *block_accounting_ptr;
		const auto& blocks = *blocks_ptr;

		std::vector<typename StateVector<state_t>::weight_type > weights; weights.reserve(this->size());
		for( const auto& sequence: m_rows ) { weights.push_back(sequence->weight()); }

		auto& frequencies = *m_w_frequencies;

		for( std::size_t n_block=0; n_block < number_of_blocks; ++n_block )
		{
			const auto n_end = ( n_block == last_block_index ? last_block_size : n_loci_per_block );
			using fq_acc_type = apegrunt::Column_frequency_accumulator<typename w_frequency_t::value_type, number_of_states<state_t>::value, apegrunt::StateBlock_size>;
			fq_acc_type column_frequency_accumulator( frequencies[n_block*n_loci_per_block].data(), n_end );
			column_frequency_accumulator.setZero();
			const auto& sequence_blocks = blocks[n_block];
			const auto& indices = block_accounting[n_block];

			for( std::size_t block_index=0; block_index < indices.size(); ++block_index )
			{
				column_frequency_accumulator.accumulate( sequence_blocks[block_index], indices[block_index], weights.data() );
			}
		}
	}

	void cache_statecount_and_statepresence_blocks() const
	{
		const auto blocks_ptr = this->get_block_storage(); // keep shared_ptr alive
		const auto& blocks = *blocks_ptr;

		m_statecount_blocks = std::make_shared< statecount_block_storage_t >(); m_statecount_blocks->reserve( blocks.size() );
		auto& statecount_blocks = *m_statecount_blocks;

		m_statepresence_blocks = std::make_shared< statepresence_block_storage_t >(); m_statepresence_blocks->reserve( blocks.size() );
		auto& statepresence_blocks = *m_statepresence_blocks;

		m_statepresence_blocks_wo_gaps = std::make_shared< statepresence_block_storage_t >(); m_statepresence_blocks_wo_gaps->reserve( blocks.size() );
		auto& statepresence_blocks_wo_gaps = *m_statepresence_blocks_wo_gaps;

		m_gappresence_blocks = std::make_shared< statepresence_block_storage_t >(); m_gappresence_blocks->reserve( blocks.size() );
		auto& gappresence_blocks = *m_gappresence_blocks;

		// some temp definitions
		statecount_block_t ones( statecount_t(1) );
		statepresence_block_t gaps( 1 << std::size_t(apegrunt::gap_state<state_t>::value) ); // fill gap mask; this will work as long as gaps are in predefined positions
		statepresence_block_t all_but_gaps( ~(1 << std::size_t(apegrunt::gap_state<state_t>::value)) );

		// create statepresence masks
		for( const auto& block_column: blocks )
		{
			statecount_block_t statepresence_mask( statecount_t(0) );
			//std::cout << " at onset: \"" << statepresence_mask << "\"" << std::endl;
			for( const auto& block: block_column )
			{
				statepresence_mask = statepresence_mask | ( ones << block );
			}
			statepresence_blocks.push_back( statepresence_mask );
			statepresence_blocks_wo_gaps.push_back( statepresence_mask & all_but_gaps );
			gappresence_blocks.push_back( statepresence_mask & gaps );
			statecount_blocks.push_back( apegrunt::popcnt_per_element(statepresence_mask) );
		}

		// zero out padding elements of the last column block
		for( std::size_t i = apegrunt::get_last_block_size(this->n_loci()); i < N; ++i ) { statepresence_blocks.back()[i] = 0; statecount_blocks.back()[i] = 0; }
	}

	/// Parser interface

	statevector_t* get_new_sequence( const std::string& id_string )
	{
		m_rows.emplace_back( make_StateVector_ptr<statevector_t>( m_block_adder, m_rows.size(), id_string ) );
		return static_cast<statevector_t*>(m_rows.back().get());
	}

	void set_nloci( std::size_t nloci ) { m_nloci = nloci; }

	void finalize() { m_block_adder->extract(); }

	// Allow parser access to private members
	ALIGNMENT_PARSER_GRAMMAR_FRIENDS(my_type)

	friend base_type; // allow base class access to private members of this class

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

