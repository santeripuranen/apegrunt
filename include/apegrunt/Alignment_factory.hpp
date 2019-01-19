/** @file Alignment_factory.hpp

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

#ifndef APEGRUNT_ALIGNMENT_FACTORY_HPP
#define APEGRUNT_ALIGNMENT_FACTORY_HPP

#include <vector>
#include <algorithm>
#include <memory>

#include "Alignment_forward.h"
#include "Alignment_interface.hpp"

#include "Loci.h"
#include "Loci_parsers.hpp"

#include "Apegrunt_utility.hpp"
#include "StateVector_mutator.hpp"

namespace apegrunt {

template< typename StateVectorT >
struct conditional_insert
{
	conditional_insert( StateVectorT *statevector, const std::vector<bool>& accept_list  )
	: m_statevector_ptr( statevector ),
	  m_accept(accept_list.cbegin()),
	  m_accept_end(accept_list.cend())
	{

	}
	~conditional_insert() { }

	template< typename StateT >
	void operator()( StateT state )
	{
		if( m_accept != m_accept_end )
		{
			*m_accept && m_statevector_ptr->push_back( state );
			++m_accept;
		}
	}

	StateVectorT *m_statevector_ptr;
	std::vector<bool>::const_iterator m_accept;
	const std::vector<bool>::const_iterator m_accept_end;

};

template< typename AlignmentT >
class Alignment_factory
{
public:
	Alignment_factory() = default;
	~Alignment_factory() = default;

	template< typename StateT >
	Alignment_ptr<StateT> operator()( Alignment_ptr<StateT> alignment, const Loci_ptr accept_list ) const
	{
		using boost::get;
		using std::cbegin; using std::cend;

		auto new_alignment = std::make_shared<AlignmentT>();
		new_alignment->set_id_string( alignment->id_string()+ ( accept_list->id_string().empty() ? "" : "."+accept_list->id_string() ) );
		new_alignment->set_loci_translation( combine(alignment->get_loci_translation(), accept_list) );
		new_alignment->set_n_original_positions( alignment->n_original_positions() );

		//std::cout << " {create new alignment"; std::cout.flush();

		//std::size_t i = 0;
		for( const auto sequence: alignment )
		{
			auto new_sequence = StateVector_mutator<typename AlignmentT::statevector_t>( new_alignment->get_new_sequence( sequence->id_string() ) );
			new_sequence.set_weight( sequence->weight() ); // transfer weights

			for( const auto locus_index: accept_list )
			{
				new_sequence( (*sequence)[locus_index] );
			}
			/*
			auto pos = cbegin(sequence);
			const auto end = cend(sequence);
			std::size_t i=0;
			for( const auto locus_index: accept_list )
			{
				while( locus_index != i ) { ++i; if( end == ++pos ) { break; } }
				//std::cout << " (" << locus_index << "," << *pos << ")"; std::cout.flush();
				new_sequence( *pos );
				//new_sequence( (*sequence)[locus_index] );
			}
			*/
		}
		//new_alignment->get_block_accounting();
		//std::cout << "}" << std::endl;
		return new_alignment;
	}

	// include selected loci
	template< typename StateT >
	Alignment_ptr<StateT> include( Alignment_ptr<StateT> alignment, const Loci_ptr include_list ) const
	{
		return this->operator()( alignment, include_list );
	}

	// exclude selected loci
	template< typename StateT >
	Alignment_ptr<StateT> exclude( Alignment_ptr<StateT> alignment, const Loci_ptr exclude_list ) const
	{
		const auto include_list = alignment->get_loci_translation() - exclude_list; // set difference
		return this->operator()( alignment, include_list );
	}

	// copy selected samples
	template< typename StateT, typename IterableT >
	Alignment_ptr<StateT> copy_selected( Alignment_ptr<StateT> alignment, const IterableT accept_list, const std::string& sample_name="sample" ) const
	{
		using boost::get;

		auto new_alignment = std::make_shared<AlignmentT>();
		new_alignment->set_id_string( alignment->id_string()+"."+sample_name);
		new_alignment->set_loci_translation( alignment->get_loci_translation() );
		new_alignment->set_n_original_positions( alignment->n_original_positions());

		for( auto seqindex: *accept_list )
		{
			auto sequence = (*alignment)[seqindex];
			auto new_sequence = StateVector_mutator<typename AlignmentT::statevector_t>( new_alignment->get_new_sequence( sequence->id_string() ) );
			//new_sequence.set_weight( sequence->weight() ); // not transferring weights here, since the set of samples changes
			for( const auto state: sequence ) {	new_sequence(state); }
		}
		//new_alignment->get_block_accounting();
		return new_alignment;
	}

	template< typename StateT >
	Alignment_ptr<StateT> transpose( const Alignment_ptr<StateT> alignment ) const
	{
		auto new_alignment = std::make_shared<AlignmentT>();
		new_alignment->set_id_string( alignment->id_string()+".transposed" );
		//new_alignment->set_n_original_positions( alignment->n_original_positions());

		for( std::size_t pos = 0; pos < alignment->n_loci(); ++pos )
		{

			auto new_sequence = StateVector_mutator< typename AlignmentT::statevector_t >( new_alignment->get_new_sequence( std::to_string(pos) ) );

			for( const auto sequence: alignment )
			{
				new_sequence( (*sequence)[pos] );
			}
		}
		return new_alignment;
	}

};

template<>
class Alignment_factory< Alignment_impl_block_compressed_storage< StateVector_impl_block_compressed_alignment_storage<triallelic_state_t> > >
{
public:
	Alignment_factory() = default;
	~Alignment_factory() = default;

	using block_storage_t = StateVector_impl_block_compressed_alignment_storage<triallelic_state_t>;
	using alignment_t = Alignment_impl_block_compressed_storage< block_storage_t >;

	template< typename SourceStateT >
	Alignment_ptr<triallelic_state_t> operator()( Alignment_ptr<SourceStateT> alignment, std::vector< std::vector< State_holder<SourceStateT> > >&& accept_list ) const
	{
		using boost::get; using std::cbegin; using std::cend;

		auto new_alignment = std::make_shared<alignment_t>();
		new_alignment->set_id_string( alignment->id_string()+"."+"4-states" );
		new_alignment->set_loci_translation( alignment->get_loci_translation() );
		new_alignment->set_n_original_positions( alignment->n_original_positions());

		for( const auto sequence: alignment )
		{
			auto new_sequence = StateVector_mutator< block_storage_t >( new_alignment->get_new_sequence( sequence->id_string() ) );

			for( const auto& accept_and_locus: zip_range( accept_list, sequence ) )
			{
				if( get<1>(accept_and_locus) == nucleic_acid_state_t::GAP ) { new_sequence(triallelic_state_t::GAP); }
				else
				{
					auto state_itr = std::find( cbegin( get<0>(accept_and_locus) ), cend( get<0>(accept_and_locus) ), get<1>(accept_and_locus) );
					if( state_itr != cend( get<0>(accept_and_locus) ) )
					{
						const std::size_t n_states = get<0>(accept_and_locus).size();
						const std::size_t state_n = std::distance( cbegin( get<0>(accept_and_locus) ), state_itr );

						new_sequence( trinary_lut[state_n] ); // always map the second-most frequent allele to mid_allele
					}
				}
			}
		}
		//new_alignment->get_block_accounting();
		return new_alignment;
	}

	Alignment_ptr<triallelic_state_t> transpose( const Alignment_ptr<triallelic_state_t> alignment ) const
	{
		auto new_alignment = std::make_shared<alignment_t>();
		new_alignment->set_id_string( alignment->id_string()+".transposed" );

		for( std::size_t pos = 0; pos < alignment->n_loci(); ++pos )
		{

			auto new_sequence = StateVector_mutator< block_storage_t >( new_alignment->get_new_sequence( std::to_string(pos) ) );

			for( const auto sequence: alignment )
			{
				new_sequence( (*sequence)[pos] );
			}
		}
		return new_alignment;
	}

	template< typename StateT >
	Alignment_ptr<StateT> reorder_alignment( const Alignment_ptr<StateT> alignment, const Loci_ptr reorder_list )
	{
		auto new_alignment = std::make_shared<alignment_t>();
		new_alignment->set_id_string( alignment->id_string()+ ( !reorder_list->id_string().empty() ? "."+reorder_list->id_string() : "" ) );
		new_alignment->set_loci_translation( combine(alignment->get_loci_translation(), reorder_list) );

		for( const auto sequence: alignment )
		{
			auto new_sequence = StateVector_mutator< block_storage_t >( new_alignment->get_new_sequence( sequence->id_string() ) );
			//auto source_sequence = sequence->subscript_proxy();
			for( const auto locus_index: reorder_list )
			{
				new_sequence( (*sequence)[locus_index] );
			}
		}
		//new_alignment->get_block_accounting();
		return new_alignment;
	}

private:

	triallelic_state_t binary_lut[2] = { triallelic_state_t::major_allele, triallelic_state_t::minor_allele };
	triallelic_state_t trinary_lut[3] = { triallelic_state_t::major_allele, triallelic_state_t::mid_allele, triallelic_state_t::minor_allele };
};

} // namespace apegrunt

#endif // APEGRUNT_ALIGNMENT_FACTORY_HPP

