/** @file StateVector_interface.hpp
 
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

#ifndef APEGRUNT_STATEVECTOR_INTERFACE_HPP
#define APEGRUNT_STATEVECTOR_INTERFACE_HPP

#include <string>
#include <functional> // for std::hash
#include <memory> // for std::unique_ptr and std::make_unique

#include <boost/serialization/serialization.hpp>
#include <boost/serialization/nvp.hpp>
#include <boost/range/iterator_range_core.hpp>

#include "StateVector_forward.h"
#include "StateVector_iterator.h"
#include "State_block.hpp"

#include "IndexVector.h"

namespace apegrunt {

/*
template< typename StateT >
class StateVector_subscript_proxy
{
public:
	StateVector_subscript_proxy() { }
	~StateVector_subscript_proxy() { }

	StateVector_subscript_proxy( const std::vector<StateT>* container ) : m_container(container) { }

	inline StateT operator[]( std::size_t index ) const { return (*m_container)[index]; }

private:
	const std::vector<StateT> *m_container;
};
*/

// A compositable "Template Method" interface class
template< typename StateT >
class StateVector
{
public:

	using state_t = StateT;
	using my_type = StateVector<state_t>;
	using const_iterator = apegrunt::iterator::StateVector_iterator<state_t>;
	using iterator = apegrunt::iterator::StateVector_iterator<state_t>;

	//using const_iterator = apegrunt::iterator::StateVector_const_iterator;
	//using iterator = apegrunt::iterator::StateVector_iterator;
	using value_type = state_t; //typename const_iterator::value_type;
	using block_type = State_block<value_type,apegrunt::StateBlock_size>;
	//using frequencies_type = std::array< std::size_t, number_of_states<state_t>::value >;
	using weight_type = double;


	using block_index_t = uint32_t; // uint16_t is sufficient only for at most 2^16 (=65536) blocks
	using block_index_container_t = std::vector< block_index_t >;
	using block_index_container_ptr = std::shared_ptr< block_index_container_t >;
	//using block_index_container_t = IndexVector< block_index_t >;

	StateVector() = default;
    virtual ~StateVector() = default; // enable derived classes to be destructed through StateVector_ptr

    inline const_iterator cbegin() const { return this->cbegin_impl(); }
    inline const_iterator cend() const { return this->cend_impl(); }

    inline value_type operator[]( std::size_t index ) const { return this->subscript_operator_impl( index ); }
    inline bool operator==( const my_type& rhs ) const { return this->equal_to_operator_impl( rhs ); }

    inline block_type get_block( std::size_t index ) const { return this->get_block_impl( index ); }
    inline const block_index_container_ptr get_block_indices() const { return this->get_block_indices_impl(); }

    //inline const frequencies_type& frequencies() const { return this->frequencies_impl(); }

    inline std::size_t size() const { return this->size_impl(); }
    inline std::size_t bytesize() const { return this->bytesize_impl(); }

    inline std::size_t id() const { return this->id_impl(); }
    // no set_id(); id is an internal index defined by the parent Alignment at StateVector construction
    inline const std::string& id_string() const { return this->id_string_impl(); }
    inline void set_id_string( const std::string& id_string ) { this->set_id_string_impl(id_string); }

    inline std::size_t multiplicity() const { return this->multiplicity_impl(); }
	inline void set_multiplicity( std::size_t multiplicity=1 ) { this->set_multiplicity_impl(multiplicity); }

    inline weight_type weight() const { return this->weight_impl(); }
	inline void set_weight( weight_type weight=1 ) { this->set_weight_impl(weight); }

	inline std::size_t operator&&( const my_type& rhs ) const { return this->logical_AND_operator( rhs ); }
	inline std::vector<bool> operator&( const my_type& rhs ) const { return this->bitwise_AND_operator( rhs ); }

	inline bool is_similar_to( const my_type& rhs, std::size_t hamming_distance_threshold ) const { return this->is_similar_to_impl( rhs, hamming_distance_threshold ); }

	inline const std::type_info& type() const { return this->type_impl(); }

	//inline StateVector_subscript_proxy<value_type> subscript_proxy() const { return this->subscript_proxy_impl(); }

private:
    //virtual StateVector_ptr clone_impl() const = 0;
    //virtual StateVector_ptr move_impl() const = 0;

    //virtual const std::string& m_get_id_string() const = 0;

    virtual const_iterator cbegin_impl() const = 0;
    virtual const_iterator cend_impl() const = 0;

    virtual value_type subscript_operator_impl( std::size_t index ) const = 0;
    virtual bool equal_to_operator_impl( const my_type& rhs ) const = 0;

    virtual block_type get_block_impl( std::size_t index ) const = 0;
    virtual const block_index_container_ptr get_block_indices_impl() const = 0;

    //virtual const frequencies_type& frequencies_impl() const = 0;

    virtual std::size_t size_impl() const = 0;
    virtual std::size_t bytesize_impl() const = 0;

    virtual std::size_t id_impl() const = 0;
    virtual const std::string& id_string_impl() const = 0;
    virtual void set_id_string_impl( const std::string& id_string ) = 0;

    virtual std::size_t multiplicity_impl() const = 0;
	virtual void set_multiplicity_impl( std::size_t multiplicity ) = 0;

    virtual double weight_impl() const = 0;
	virtual void set_weight_impl( double weight ) = 0;

	virtual std::size_t logical_AND_operator( const my_type& rhs ) const = 0;
	virtual std::vector<bool> bitwise_AND_operator( const my_type& rhs ) const = 0;

	virtual bool is_similar_to_impl( const my_type& rhs, std::size_t hamming_distance_threshold ) const = 0;

	virtual const std::type_info& type_impl() const = 0;

	//virtual StateVector_subscript_proxy<value_type> subscript_proxy_impl() const = 0;
};

template< typename StateT >
std::ostream& operator<< ( std::ostream& os, const StateVector<StateT>* statevector )
{
	os << statevector->id_string();
	return os;
}

template< typename StateT >
std::ostream& operator<< ( std::ostream& os, const StateVector_ptr<StateT>& statevector )
{
	os << statevector->id_string();
	return os;
}

template< typename StateT > typename StateVector<StateT>::const_iterator begin( const StateVector<StateT>& statevector ) { return statevector.cbegin(); }
template< typename StateT > typename StateVector<StateT>::const_iterator end( const StateVector<StateT>& statevector ) { return statevector.cend(); }

template< typename StateT > typename StateVector<StateT>::const_iterator begin( const StateVector_ptr<StateT>& statevector ) { return statevector->cbegin(); }
template< typename StateT > typename StateVector<StateT>::const_iterator end( const StateVector_ptr<StateT>& statevector ) { return statevector->cend(); }

template< typename StateT > typename StateVector<StateT>::const_iterator cbegin( const StateVector<StateT>& statevector ) { return statevector.cbegin(); }
template< typename StateT > typename StateVector<StateT>::const_iterator cend( const StateVector<StateT>& statevector ) { return statevector.cend(); }

template< typename StateT > typename StateVector<StateT>::const_iterator cbegin( const StateVector_ptr<StateT>& statevector ) { return statevector->cbegin(); }
template< typename StateT > typename StateVector<StateT>::const_iterator cend( const StateVector_ptr<StateT>& statevector ) { return statevector->cend(); }

template< typename StateT >
bool are_similar( const StateVector<StateT>* lhs, const StateVector<StateT>* rhs, std::size_t hamming_distance_threshold ) { return lhs->is_similar_to( *rhs, hamming_distance_threshold ); }

template< typename StateT >
bool are_similar( const StateVector_ptr<StateT> lhs, const StateVector_ptr<StateT> rhs, std::size_t hamming_distance_threshold ) { return lhs->is_similar_to( *rhs, hamming_distance_threshold ); }

/*
std::vector<bool> operator&( const StateVector& lhs, const StateVector& rhs )
{
	using boost::get;

	const std::size_t N = std::min(lhs.size(), rhs.size());
	std::vector<bool> boolvec; boolvec.reserve(N);
	for( auto pair: zip_range( lhs, rhs ) )
	{
		boolvec.push_back( get<0>(pair) == get<1>(pair) );
	}
	return boolvec;
}
*/
} // namespace apegrunt

#endif // APEGRUNT_STATEVECTOR_INTERFACE_HPP
