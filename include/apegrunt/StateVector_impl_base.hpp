/** @file StateVector_impl_base.hpp

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

#ifndef APEGRUNT_STATEVECTOR_IMPL_BASE_HPP
#define APEGRUNT_STATEVECTOR_IMPL_BASE_HPP

#include "StateVector.h"

namespace apegrunt {

template< typename StateVectorT, typename StateT >
class StateVector_impl_base : public StateVector< StateT >
{
public:
	using state_t = StateT;
	using const_iterator = typename StateVector<state_t>::const_iterator;
	using iterator = typename StateVector<state_t>::iterator;
	using value_type = state_t; //typename const_iterator::value_type;
	using block_type = typename StateVector<state_t>::block_type;
	//using frequencies_type = typename StateVector<state_t>::frequencies_type;
	using weight_type = typename StateVector<state_t>::weight_type;

	using block_index_t = typename StateVector<state_t>::block_index_t;
	using block_index_container_t = typename StateVector<state_t>::block_index_container_t;
	using block_index_container_ptr = typename StateVector<state_t>::block_index_container_ptr;

	enum { block_size = block_type::N };

	StateVector_impl_base() : m_multiplicity(1), m_weight(1) { }
	virtual ~StateVector_impl_base() = default;

	//StateVector_impl_base( const std::string& id_string, std::size_t multiplicity=1, weight_type weight=1 ) : m_multiplicity(multiplicity), m_weight(weight), m_id_string(id_string) { }
	StateVector_impl_base( std::size_t id ) : m_id(id), m_multiplicity(1), m_weight(1), m_id_string() { }
	StateVector_impl_base( std::size_t id, const std::string& id_string ) : m_id(id), m_multiplicity(1), m_weight(1), m_id_string(id_string) { }

	inline std::size_t bytesize() const { return this->size_impl()*sizeof(state_t); }

	inline std::size_t id() const { return m_id; }
	inline void set_id( std::size_t id ) { m_id = id; }
	inline const std::string& id_string() const { return m_id_string; }
	inline void set_id_string( const std::string& id_string ) { m_id_string = id_string; }

	std::size_t multiplicity() const { return m_multiplicity; }
	void set_multiplicity( std::size_t multiplicity=1 ) { m_multiplicity = multiplicity; } // zero means we're dead

	weight_type weight() const { return m_weight; }
	void set_weight( weight_type weight=1 ) { m_weight = weight; }

	//bool operator==( const StateVectorT& other ) const { std::cout << "StateVector_impl_base::operator==( const derived_type& other )" << std::endl; return static_cast<const_ref_cast_t>(*this).operator==(other); }

	// default implementation; used when derived_type does not implement operator==( const StateVector& other )
	bool operator==( const StateVector<state_t>& rhs ) const
	{
	   	using boost::get;
		//std::cout << "StateVector_impl_base: operator==( const StateVector& rhs )" << std::endl;
		if( static_cast<const_ref_cast_t>(*this).size() == rhs.size() )
		{
			for( auto state_pair: zip_range( *this, rhs ) )
			{
				if( get<0>(state_pair) != get<1>(state_pair) ) { return false; }
			}
		}
		return true;
	}

    inline block_type get_block( std::size_t index ) const
    {
    	auto block = block_type();

    	const auto beg = index*block_size;
    	if( beg+block_size-1 < this->size() )
    	{
    		for( std::size_t i=0; i<block_size; ++i ) { block[i] = (*static_cast<const_cast_t>(this))[beg+i]; }
    	}
    	else
    	{
    		const auto end = this->size() % block_size;
    		for( std::size_t i=0; i<end; ++i ) { block[i] = (*static_cast<const_cast_t>(this))[beg+i]; }
    	}
    	return block;
    }

	inline std::size_t operator&&( const StateVector<state_t>& rhs ) const
	{
   		std::cout << "StateVector_impl_base::operator&&( const StateVector& other ) -- FALLBACK" << std::endl;

   		using boost::get;

		std::size_t N = 0;
		for( auto pair: zip_range( static_cast<const_ref_cast_t>(*this), rhs ) )
		{
			get<0>(pair) == get<1>(pair) && ++N;
		}
		return N;
	}

	inline std::vector<bool> operator&( const StateVector<state_t>& rhs ) const
	{
    	std::cout << "StateVector_impl_base::operator&( const StateVector& rhs ) -- FALLBACK" << std::endl;

		using boost::get;

		const std::size_t N = std::min(static_cast<const_ref_cast_t>(*this).size(), rhs.size());
		std::vector<bool> boolvec; boolvec.reserve(N);
		for( auto pair: zip_range( static_cast<const_ref_cast_t>(*this), rhs ) )
		{
			boolvec.push_back( get<0>(pair) == get<1>(pair) );
		}
		return boolvec;
	}

	inline bool is_similar_to( const StateVector<state_t>& rhs, std::size_t hamming_distance_threshold ) const
	{
    	std::cout << "StateVector_impl_base::is_similar_to( const StateVector& rhs ) -- FALLBACK" << std::endl;
   		using boost::get;

   		std::size_t Npos = std::min( this->size(), rhs.size() );
		std::size_t N(0), n(0);

		for( auto pair: zip_range( static_cast<const_ref_cast_t>(*this), rhs ) )
		{
			++n;
			get<0>(pair) == get<1>(pair) && ++N;
			if( N+(Npos-n) < hamming_distance_threshold ) { return false; }
		}

		return (N > hamming_distance_threshold); // ? true: false;
	}

private:
	using derived_type = StateVectorT;
	using cast_t = derived_type* const;
	using const_cast_t = const derived_type* const;
	using const_ref_cast_t = const derived_type&;

	/*const*/ std::size_t m_id;
	mutable std::size_t m_multiplicity; // the number of virtual sequences (identical, or near identical) that this sequence represents
	weight_type m_weight;
	std::string m_id_string;

	//StateVector_ptr clone_impl() const override { return static_cast<const_cast_t>(this)->clone(); }

	const_iterator cbegin_impl() const override { return const_iterator( static_cast<const_cast_t>(this)->cbegin() ); }
	const_iterator cend_impl() const override { return const_iterator( static_cast<const_cast_t>(this)->cend() ); }

	iterator begin_impl() { return iterator( static_cast<cast_t>(this)->begin() ); }
	iterator end_impl() { return iterator( static_cast<cast_t>(this)->end() ); }

    value_type subscript_operator_impl( std::size_t index ) const override { return (*static_cast<const_cast_t>(this))[index]; }
    bool equal_to_operator_impl( const StateVector<state_t>& rhs ) const override
    {
    	//std::cout << "StateVector_impl_base::equal_to_operator_impl( const StateVector& rhs )" << std::endl;
    	if( static_cast<const_ref_cast_t>(*this).type() == rhs.type() )
    	{
    		return static_cast<const_ref_cast_t>(*this) == static_cast<const_ref_cast_t>(rhs);
    	}
    	else
    	{
    		// redirect to default implementation in StateVector_impl_base, unless derived_type implements operator==( const StateVector& rhs )
    		//return static_cast<const_ref_cast_t>(*this) == rhs;
    		return *this == rhs;
    	}
    }

    block_type get_block_impl( std::size_t index ) const override { return static_cast<const_cast_t>(this)->get_block(index); }
    const block_index_container_ptr get_block_indices_impl() const { return static_cast<const_cast_t>(this)->get_block_indices(); }

    //const frequencies_type& frequencies_impl() const override { return static_cast<const_cast_t>(this)->frequencies(); }

    inline std::size_t size_impl() const override { return static_cast<const_cast_t>(this)->size(); }
    inline std::size_t bytesize_impl() const override { return static_cast<const_cast_t>(this)->bytesize(); }

    inline std::size_t id_impl() const override { return static_cast<const_cast_t>(this)->id(); }
    inline const std::string& id_string_impl() const override { return static_cast<const_cast_t>(this)->id_string(); }
    inline void set_id_string_impl( const std::string& id_string ) override { static_cast<cast_t>(this)->set_id_string(id_string); }

    std::size_t multiplicity_impl() const override { return static_cast<const_cast_t>(this)->multiplicity(); }
	void set_multiplicity_impl( std::size_t multiplicity ) override { static_cast<cast_t>(this)->set_multiplicity(multiplicity); }

    weight_type weight_impl() const override { return static_cast<const_cast_t>(this)->weight(); }
	void set_weight_impl( weight_type weight ) override { static_cast<cast_t>(this)->set_weight(weight); }

	std::size_t logical_AND_operator( const StateVector<state_t>& rhs ) const override
	{
    	if( static_cast<const_ref_cast_t>(*this).type() == rhs.type() )
    	{
    		return static_cast<const_ref_cast_t>(*this) && static_cast<const_ref_cast_t>(rhs);
    	}
    	else
    	{
    		// redirect to default implementation in StateVector_impl_base, unless derived_type implements operator&&( const StateVector& rhs )
    		//return static_cast<const_ref_cast_t>(*this) && rhs;
    		return *this && rhs;
    	}
	}

	std::vector<bool> bitwise_AND_operator( const StateVector<state_t>& rhs ) const override
	{
    	if( static_cast<const_ref_cast_t>(*this).type() == rhs.type() )
    	{
    		return static_cast<const_ref_cast_t>(*this) & static_cast<const_ref_cast_t>(rhs);
    	}
    	else
    	{
    		// redirect to default implementation in StateVector_impl_base, unless derived_type implements operator&( const StateVector& rhs )
    		//return static_cast<const_ref_cast_t>(*this) & rhs;
    		return *this & rhs;
    	}
	}

	bool is_similar_to_impl( const StateVector<state_t>& rhs, std::size_t hamming_distance_threshold ) const override
	{
    	if( static_cast<const_ref_cast_t>(*this).type() == rhs.type() )
    	{
    		return static_cast<const_ref_cast_t>(*this).is_similar_to( static_cast<const_ref_cast_t>(rhs), hamming_distance_threshold );
    	}
    	else
    	{
    		// redirect to default implementation in StateVector_impl_base, unless derived_type implements operator&&( const StateVector& rhs )
    		//return static_cast<const_ref_cast_t>(*this) && rhs;
    		return this->is_similar_to( rhs, hamming_distance_threshold );
    	}
	}

	const std::type_info& type_impl() const override { return static_cast<const_cast_t>(this)->type(); }

    //StateVector_subscript_proxy<value_type> subscript_proxy_impl() const override { return static_cast<const_ref_cast_t>(*this).subscript_proxy(); }

    friend class boost::serialization::access;
	/// boost.serialization interface.
	template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
		ar & BOOST_SERIALIZATION_NVP(m_id_string);
    }
};

} // namespace apegrunt

#endif // APEGRUNT_STATEVECTOR_IMPL_BASE_HPP

