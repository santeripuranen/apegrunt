/** @file Alignment_interface.hpp
 
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
#ifndef APEGRUNT_ALIGNMENT_INTERFACE_HPP
#define APEGRUNT_ALIGNMENT_INTERFACE_HPP

#include <string>
#include <functional> // for std::hash
#include <memory> // for std::unique_ptr and std::make_unique
#include <vector>

#include <boost/serialization/serialization.hpp>
#include <boost/serialization/nvp.hpp>
//#include <boost/range/iterator_range_core.hpp>

#include "Alignment_forward.h"
#include "Alignment_iterator.h"
#include "StateVector_forward.h"
#include "Loci_forward.h"
#include "State_block.hpp"

#include "apegrunt/IntegerSequence_Hybrid_bitset_range.h"
#include "misc/Vector.h"

#include "aligned_allocator.hpp"

namespace apegrunt {

template< typename StateVectorT >
class Alignment_subscript_proxy
{
public:
	Alignment_subscript_proxy() { }
	~Alignment_subscript_proxy() { }

	Alignment_subscript_proxy( const std::vector<StateVectorT>* container ) : m_container(container) { }

	inline StateVectorT operator[]( std::size_t index ) const { return (*m_container)[index]; }

private:
	const std::vector<StateVectorT> *m_container;
};

//> The Alignment interface class
template< typename StateT >
class Alignment
{
public:
	using state_t = StateT;
	using real_t = double;

	using const_iterator = apegrunt::iterator::Alignment_const_iterator<state_t>;
	using iterator = apegrunt::iterator::Alignment_iterator<state_t>;
	using value_type = typename const_iterator::value_type;

	using frequency_t = std::array< std::size_t, number_of_states<state_t>::value >;
	using frequencies_t = std::vector<frequency_t>;
	using frequencies_ptr = std::shared_ptr< frequencies_t >;

	using w_frequency_t = std::array< real_t, number_of_states<state_t>::value >;
	using w_frequencies_t = std::vector<w_frequency_t>;
	using w_frequencies_ptr = std::shared_ptr< w_frequencies_t >;

	using distance_matrix_t = std::vector<uint32_t>;
	using distance_matrix_ptr = std::shared_ptr< distance_matrix_t >;

	using block_weight_t = std::array< real_t, number_of_states<state_t>::value >;
	using block_weights_t = std::vector< std::vector<block_weight_t> >;
	using block_weights_ptr = std::shared_ptr< block_weights_t >;

	using block_index_t = uint32_t; // uint16_t is only good for alignments up to 2^16 (=65536) samples
	//using block_accounting_container_t = std::vector<block_index_t>;
	using block_accounting_container_t = apegrunt::Apegrunt_bitset<block_index_t>;
	using block_accounting_t = std::vector< std::vector< block_accounting_container_t > >;
	using block_accounting_ptr = std::shared_ptr< block_accounting_t >;

	using block_indices_t = std::vector< std::size_t >;
	using block_indices_ptr = std::shared_ptr< block_indices_t >;

	using block_type = State_block< State_holder<state_t>, apegrunt::StateBlock_size >;
//	using compressed_block_type = Compressed_state_block< State_holder<state_t>, apegrunt::StateBlock_size >;

	using allocator_t = memory::AlignedAllocator<block_type,alignof(block_type)>;
	using block_storage_t = std::vector< std::vector< block_type, allocator_t > >;
	using block_storage_ptr = std::shared_ptr< block_storage_t >;

	using statecount_t = uint8_t;
	using statecount_block_t = apegrunt::Vector< statecount_t, apegrunt::StateBlock_size >; //std::array< statecount_t, apegrunt::StateBlock_size >;
	using statecount_block_storage_t = std::vector< statecount_block_t >;
	using statecount_block_storage_ptr = std::shared_ptr< statecount_block_storage_t >;

	using statepresence_t = uint8_t;
	using statepresence_block_t = apegrunt::Vector< statepresence_t, apegrunt::StateBlock_size >; //std::array< statepresence_t, apegrunt::StateBlock_size >;
	using statepresence_block_storage_t = std::vector< statepresence_block_t >;
	using statepresence_block_storage_ptr = std::shared_ptr< statepresence_block_storage_t >;

	Alignment() = default;
    virtual ~Alignment() = default; // enable derived classes to be destructed through Alignment_ptr

    Alignment_ptr<state_t> clone() const { return this->clone_impl(); } // clone the object, leaving the callee valid and unmodified.

    inline iterator begin() { return this->begin_impl(); }
    inline iterator end() { return this->end_impl(); }

    inline const_iterator cbegin() const { return this->cbegin_impl(); }
    inline const_iterator cend() const { return this->cend_impl(); }

    inline value_type operator[]( std::size_t index ) const { return this->square_bracket_operator_impl( index ); }

    //> Return the number of sequences contained in the alignment
    inline std::size_t size() const { return this->size_impl(); }

    //> Return the effective number of sequences contained in the alignment. Weights in sequence multiplicity (fused duplicates increase multiplicity), such that effective_size() >= size().
    inline typename w_frequency_t::value_type effective_size() const { return this->effective_size_impl(); }

    inline const std::string& id_string() const { return this->id_string_impl(); }
	inline void set_id_string( const std::string& id_string ) { this->set_id_string_impl(id_string); }

    //> Return the number of (least common denominator) columns in the contained alignment
    inline std::size_t n_loci() const { return this->n_loci_impl(); }
    inline std::size_t n_original_positions() const { return this->n_original_positions_impl(); }
    inline void set_n_original_positions( std::size_t npositions ) { this->set_n_original_positions_impl( npositions ); }

    inline frequencies_ptr frequencies() const { return this->frequencies_impl(); }
    inline w_frequencies_ptr w_frequencies() const { return this->w_frequencies_impl(); }

    inline distance_matrix_ptr distance_matrix() const { return this->distance_matrix_impl(); }

    inline const std::type_info& type() const { return this->type_impl(); }

    inline void set_loci_translation( Loci_ptr translation_table ) { return this->set_loci_translation_impl( translation_table ); }
    inline Loci_ptr get_loci_translation() { return this->get_loci_translation_impl(); }

	inline void fuse_duplicates() { this->fuse_duplicates_impl(); }

	inline Alignment_subscript_proxy< StateVector_ptr<state_t> > subscript_proxy() const { return this->subscript_proxy_impl(); }

	inline void statistics( std::ostream* out=nullptr ) const { this->statistics_impl( out ); }

	inline block_accounting_ptr get_block_accounting() const { return this->get_block_accounting_impl(); }
	inline block_storage_ptr get_block_storage() const { return this->get_block_storage_impl(); }
	inline block_weights_ptr get_block_weights() const { return this->get_block_weights_impl(); }

	inline block_indices_ptr get_block_indices() const { return this->get_block_indices_impl(); }

	inline statecount_block_storage_ptr get_statecount_blocks() { return this->get_statecount_blocks_impl(); }
	inline statepresence_block_storage_ptr get_statepresence_blocks() { return this->get_statepresence_blocks_impl(); }
	inline statepresence_block_storage_ptr get_statepresence_blocks_wo_gaps() { return this->get_statepresence_blocks_wo_gaps_impl(); }
	inline statepresence_block_storage_ptr get_gappresence_blocks() { return this->get_gappresence_blocks_impl(); }
private:

    virtual Alignment_ptr<state_t> clone_impl() const = 0;

    virtual iterator begin_impl() = 0;
    virtual iterator end_impl() = 0;

    virtual const_iterator cbegin_impl() const = 0;
    virtual const_iterator cend_impl() const = 0;

    virtual value_type square_bracket_operator_impl( std::size_t index ) const = 0;

    virtual std::size_t size_impl() const = 0;
    virtual typename w_frequency_t::value_type effective_size_impl() const = 0;

    virtual const std::string& id_string_impl() const = 0;
    virtual void set_id_string_impl( const std::string& id_string ) = 0;

    virtual std::size_t n_loci_impl() const = 0;
    virtual std::size_t n_original_positions_impl() const = 0;
    virtual void set_n_original_positions_impl( std::size_t npositions ) = 0;

    virtual frequencies_ptr frequencies_impl() const = 0;
    virtual w_frequencies_ptr w_frequencies_impl() const = 0;

    virtual distance_matrix_ptr distance_matrix_impl() const = 0;

    virtual const std::type_info& type_impl() const = 0;

    virtual void set_loci_translation_impl( Loci_ptr translation_table ) = 0;
    virtual Loci_ptr get_loci_translation_impl() = 0;

    virtual void fuse_duplicates_impl() = 0;

    virtual Alignment_subscript_proxy< StateVector_ptr<state_t> > subscript_proxy_impl() const = 0;

    virtual void statistics_impl( std::ostream *out ) const = 0;

    virtual block_accounting_ptr get_block_accounting_impl() const = 0;
    virtual block_storage_ptr get_block_storage_impl() const = 0;
    virtual block_weights_ptr get_block_weights_impl() const = 0;

    virtual block_indices_ptr get_block_indices_impl() const = 0;

	virtual statecount_block_storage_ptr get_statecount_blocks_impl() const = 0;
	virtual statepresence_block_storage_ptr get_statepresence_blocks_impl() const = 0;
	virtual statepresence_block_storage_ptr get_statepresence_blocks_wo_gaps_impl() const = 0;
	virtual statepresence_block_storage_ptr get_gappresence_blocks_impl() const = 0;
};

template< typename StateT > typename Alignment<StateT>::const_iterator begin( const Alignment<StateT>& alignment ) { return alignment.cbegin(); }
template< typename StateT > typename Alignment<StateT>::const_iterator end( const Alignment<StateT>& alignment ) { return alignment.cend(); }

template< typename StateT > typename Alignment<StateT>::const_iterator begin( const Alignment_ptr<StateT>& alignment ) { return alignment->cbegin(); }
template< typename StateT > typename Alignment<StateT>::const_iterator end( const Alignment_ptr<StateT>& alignment ) { return alignment->cend(); }

template< typename StateT > typename Alignment<StateT>::const_iterator cbegin( const Alignment<StateT>& alignment ) { return alignment.cbegin(); }
template< typename StateT > typename Alignment<StateT>::const_iterator cend( const Alignment<StateT>& alignment ) { return alignment.cend(); }

template< typename StateT > typename Alignment<StateT>::const_iterator cbegin( const Alignment_ptr<StateT>& alignment ) { return alignment->cbegin(); }
template< typename StateT > typename Alignment<StateT>::const_iterator cend( const Alignment_ptr<StateT>& alignment ) { return alignment->cend(); }

} // namespace apegrunt

#endif // APEGRUNT_ALIGNMENT_INTERFACE_HPP
