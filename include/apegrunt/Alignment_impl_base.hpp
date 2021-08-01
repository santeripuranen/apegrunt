/** @file Alignment_impl_base.hpp

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

#ifndef APEGRUNT_ALIGNMENT_IMPL_BASE_HPP
#define APEGRUNT_ALIGNMENT_IMPL_BASE_HPP

#include <numeric> // for std::accumulate
#include <vector>
#include <memory>  // for std::shared_ptr

#include <boost/serialization/serialization.hpp>
#include <boost/serialization/nvp.hpp>

#include "Alignment.h"
#include "StateVector.h"
#include "StateVector_mutator.hpp"
#include "StateVector_utility.hpp"
#include "Apegrunt_utility.hpp"
#include "Loci_parsers.hpp"

namespace apegrunt {

template< typename AlignmentT, typename StateT >
class Alignment_impl_base : public Alignment< StateT > //Alignment< typename AlignmentT::state_t >
{
public:
	using state_t = StateT;
	//using state_t = typename AlignmentT::state_t;
	using const_iterator = apegrunt::iterator::Alignment_const_iterator<state_t>;
	using iterator = apegrunt::iterator::Alignment_iterator<state_t>;
	using value_type = typename const_iterator::value_type;

	using frequency_t = typename Alignment<state_t>::frequency_t;
	using frequencies_t = typename Alignment<state_t>::frequencies_t;
	using frequencies_ptr = typename Alignment<state_t>::frequencies_ptr;

	using w_frequency_t = typename Alignment<state_t>::w_frequency_t;
	using w_frequencies_t = typename Alignment<state_t>::w_frequencies_t;
	using w_frequencies_ptr = typename Alignment<state_t>::w_frequencies_ptr;

	using distance_matrix_t = typename Alignment<state_t>::distance_matrix_t;
	using distance_matrix_ptr = typename Alignment<state_t>::distance_matrix_ptr;

	using block_index_t = typename Alignment<state_t>::block_index_t;
	using block_accounting_container_t = typename Alignment<state_t>::block_accounting_container_t;
	using block_accounting_t = typename Alignment<state_t>::block_accounting_t;
	using block_accounting_ptr = typename Alignment<state_t>::block_accounting_ptr;

	using block_container_t = typename Alignment<state_t>::block_container_t;
	using block_storage_t = typename Alignment<state_t>::block_storage_t;
	using block_storage_ptr = typename Alignment<state_t>::block_storage_ptr;

	using block_type = typename Alignment<state_t>::block_type;

	using block_indices_t = typename Alignment<state_t>::block_indices_t;
	using block_indices_ptr = typename Alignment<state_t>::block_indices_ptr;

	using statecount_t = typename Alignment<state_t>::statecount_t;
	using statecount_block_t = typename Alignment<state_t>::statecount_block_t;
	using statecount_block_storage_t = typename Alignment<state_t>::statecount_block_storage_t;
	using statecount_block_storage_ptr = typename Alignment<state_t>::statecount_block_storage_ptr;

	using statepresence_t = typename Alignment<state_t>::statepresence_t;
	using statepresence_block_t = typename Alignment<state_t>::statepresence_block_t;
	using statepresence_block_storage_t = typename Alignment<state_t>::statepresence_block_storage_t;
	using statepresence_block_storage_ptr = typename Alignment<state_t>::statepresence_block_storage_ptr;

	Alignment_impl_base() = default;
	virtual ~Alignment_impl_base() override = default;

	Alignment_impl_base( const std::string& id_string ) : m_id_string(id_string) { }

	const std::string& id_string() const { return m_id_string; }
	void set_id_string( const std::string& id_string ) { m_id_string = id_string; }
	void set_id_string( std::string&& id_string ) { m_id_string = std::move(id_string); }

	Alignment_ptr<state_t> subset( Loci_ptr positions, std::ostream *out=nullptr ) const
	{
		using boost::get;

		if( positions->size() == 0 ) { return Alignment_ptr<StateT>(); }

		auto subset_alignment = std::make_shared<AlignmentT>();
		subset_alignment->set_id_string( this->id_string() + ( positions->id_string().empty() ? "" : "."+positions->id_string() ) );
		subset_alignment->set_loci_translation( apegrunt::combine(this->get_loci_translation(), positions) );
		subset_alignment->set_n_original_positions( this->n_original_positions() );

		for( const auto sequence: *this )
		{
			auto subset_sequence = StateVector_mutator<typename AlignmentT::statevector_t>( subset_alignment->get_new_sequence( sequence->id_string() ) );
			subset_sequence.set_weight( sequence->weight() ); // transfer weights

			for( const auto locus_index: positions )
			{
				subset_sequence( (*sequence)[locus_index] );
			}
		}

		return subset_alignment;
	}

    Alignment_ptr<state_t> subsample( const Loci_ptr samples, std::ostream *out=nullptr ) const
    {
		using boost::get;

		if( samples->size() == 0 ) { return Alignment_ptr<StateT>(); }

		auto subsample_alignment = std::make_shared<AlignmentT>();
		subsample_alignment->set_id_string( this->id_string() + ( samples->id_string().empty() ? "" : "."+samples->id_string() ) );
		subsample_alignment->set_loci_translation( this->get_loci_translation() );
		subsample_alignment->set_n_original_positions( this->n_original_positions() );

		for( auto seqindex: *samples )
		{
			auto sequence = (*this)[seqindex];
			auto new_sequence = StateVector_mutator<typename AlignmentT::statevector_t>( subsample_alignment->get_new_sequence( sequence->id_string() ) );
			//new_sequence.set_weight( sequence->weight() ); // not transferring weights here, since the set of samples changes
			for( const auto state: sequence ) {	new_sequence(state); }
		}

		return subsample_alignment;
    }

	typename w_frequency_t::value_type effective_size() const
	{
		using real_t = typename w_frequency_t::value_type;
		return std::accumulate( this->cbegin_impl(), this->cend_impl(),real_t(0), [=](auto sum, const auto seq) { return sum += real_t(seq->multiplicity())*real_t(seq->weight()); } );
	}

	void fuse_duplicates()
	{
		// find pair-wise duplicate sequences, move multiplicity contribution to the first and set the second to zero
		using boost::get;
		for( auto sv_i: *static_cast<cast_t>(this) )
		{
			if( sv_i->multiplicity() == 0 ) { continue; }
			for( auto sv_j = this->begin_impl(); *sv_j != sv_i; ++sv_j )
			{
				if( 0 == (*sv_j)->multiplicity() ) { continue; }

				//if( *sv_i == *(*sv_j) )
				if( identical(sv_i,*sv_j) )
				{
					sv_i->set_multiplicity( sv_i->multiplicity() + (*sv_j)->multiplicity() );
					(*sv_j)->set_multiplicity( 0 ); // set to null == tag for removal (below)
				}
			}
		}

		// partition by multiplicity -- all null go last
		auto first_null = std::stable_partition( this->begin_impl(), this->end_impl(), [](const auto& sv){ return sv->multiplicity() > 0; } );

		// remove all duplicate entries (those that have null multiplicity)
		this->erase_impl( first_null, this->end_impl() );
	}

    inline void set_loci_translation( Loci_ptr translation_table ) { m_loci_translation_table = translation_table; }
    inline Loci_ptr get_loci_translation() const
    {
    	if( !m_loci_translation_table )
    	{
    		std::vector< std::size_t > loci; loci.reserve(this->n_loci_impl());
    		for( std::size_t i=0; i < this->n_loci_impl(); ++i ) { loci.push_back(i); }
    		m_loci_translation_table = make_Loci_list(std::move(loci),0);
    	}

    	return m_loci_translation_table;
    }

    void statistics( std::ostream *out ) const
    {
    	// dummy implementation
    }

    inline block_accounting_ptr get_block_accounting() { return block_accounting_ptr(); } // default implementation returns empty shared_ptr
    inline block_storage_ptr get_block_storage() { return block_storage_ptr(); } // default implementation returns empty shared_ptr
    inline block_indices_ptr get_block_indices() { return block_indices_ptr(); } // default implementation returns empty shared_ptr
    inline statecount_block_storage_ptr get_statecount_blocks() { return statecount_block_storage_ptr(); } // default implementation returns empty shared_ptr
    inline statepresence_block_storage_ptr get_statepresence_blocks() { return statepresence_block_storage_ptr(); } // default implementation returns empty shared_ptr
    inline statepresence_block_storage_ptr get_statepresence_blocks_wo_gaps() { return statepresence_block_storage_ptr(); } // default implementation returns empty shared_ptr
    inline statepresence_block_storage_ptr get_gappresence_blocks() { return statepresence_block_storage_ptr(); } // default implementation returns empty shared_ptr

    inline std::size_t n_original_positions() const { return m_n_original_positions; }
    inline void set_n_original_positions( std::size_t npositions ) { m_n_original_positions = npositions; }

private:
	using derived_type = AlignmentT;
	using cast_t = derived_type* const;
	using const_cast_t = const derived_type* const;
	using const_ref_cast_t = const derived_type&;

	std::string m_id_string;
	mutable Loci_ptr m_loci_translation_table;
	std::size_t m_n_original_positions;

	virtual Alignment_ptr<state_t> clone_impl() const { return static_cast<const_cast_t>(this)->clone(); }
	//virtual Alignment_ptr move_impl() const { return static_cast<const_cast_t>(this)->move(); }

	iterator begin_impl() override { return static_cast<cast_t>(this)->begin(); }
	iterator end_impl() override { return static_cast<cast_t>(this)->end(); }

	const_iterator cbegin_impl() const override { return static_cast<const_cast_t>(this)->cbegin(); }
	const_iterator cend_impl() const override { return static_cast<const_cast_t>(this)->cend(); }

    value_type square_bracket_operator_impl( std::size_t index ) const override { return (*static_cast<const_cast_t>(this))[index]; }

    Alignment_ptr<state_t> subset_impl( const Loci_ptr positions, std::ostream *out=nullptr ) const override { return static_cast<const_cast_t>(this)->subset(positions,out); }
    Alignment_ptr<state_t> subsample_impl( const Loci_ptr samples, std::ostream *out=nullptr ) const override { return static_cast<const_cast_t>(this)->subsample(samples,out); }

    std::size_t size_impl() const override { return static_cast<const_cast_t>(this)->size(); }
    typename w_frequency_t::value_type effective_size_impl() const override { return static_cast<const_cast_t>(this)->effective_size(); }

    const std::string& id_string_impl() const override { return static_cast<const_cast_t>(this)->id_string(); };
    void set_id_string_impl( const std::string& id_string ) override { static_cast<cast_t>(this)->set_id_string(id_string); }

    std::size_t n_loci_impl() const override { return static_cast<const_cast_t>(this)->n_loci(); }
    //std::size_t get_index_of_impl( const StateVector_ptr& query ) const override { return static_cast<const_cast_t>(this)->get_index_of(query); }
    std::size_t n_original_positions_impl() const override { return static_cast<const_cast_t>(this)->n_original_positions(); }
    void set_n_original_positions_impl( std::size_t npositions ) override { static_cast<cast_t>(this)->set_n_original_positions( npositions ); }

    frequencies_ptr frequencies_impl() const override { return static_cast<const_cast_t>(this)->frequencies(); }
    w_frequencies_ptr w_frequencies_impl() const override { return static_cast<const_cast_t>(this)->w_frequencies(); }

    distance_matrix_ptr distance_matrix_impl() const override { return static_cast<const_cast_t>(this)->distance_matrix(); }

    const std::type_info& type_impl() const override { return static_cast<const_cast_t>(this)->type(); }

    void set_loci_translation_impl( Loci_ptr translation_table ) override { return static_cast<cast_t>(this)->set_loci_translation( translation_table ); }
    Loci_ptr get_loci_translation_impl() const override { return static_cast<const_cast_t>(this)->get_loci_translation(); }

    void fuse_duplicates_impl() override { static_cast<cast_t>(this)->fuse_duplicates(); }

    Alignment_subscript_proxy< StateVector_ptr<state_t> > subscript_proxy_impl() const override { return static_cast<const_ref_cast_t>(*this).subscript_proxy(); }

    void statistics_impl( std::ostream *out ) const override { static_cast<const_cast_t>(this)->statistics( out ); }

    block_accounting_ptr get_block_accounting_impl() const override { return static_cast<const_cast_t>(this)->get_block_accounting(); }
    block_storage_ptr get_block_storage_impl() const override { return static_cast<const_cast_t>(this)->get_block_storage(); }

    block_indices_ptr get_block_indices_impl() const override { return static_cast<const_cast_t>(this)->get_block_indices(); }

    statecount_block_storage_ptr get_statecount_blocks_impl() const override { return static_cast<const_cast_t>(this)->get_statecount_blocks(); }
    statepresence_block_storage_ptr get_statepresence_blocks_impl() const override { return static_cast<const_cast_t>(this)->get_statepresence_blocks(); }
    statepresence_block_storage_ptr get_statepresence_blocks_wo_gaps_impl() const override { return static_cast<const_cast_t>(this)->get_statepresence_blocks_wo_gaps(); }
    statepresence_block_storage_ptr get_gappresence_blocks_impl() const override { return static_cast<const_cast_t>(this)->get_gappresence_blocks(); }

    iterator erase_impl( iterator first, iterator last ) { return static_cast<cast_t>(this)->erase(first,last); }

    friend class boost::serialization::access;
	/// boost.serialization interface.
	template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
	//	ar & BOOST_SERIALIZATION_NVP(m_alignment);
    }
};

} // namespace apegrunt

#endif // APEGRUNT_ALIGNMENT_IMPL_BASE_HPP

