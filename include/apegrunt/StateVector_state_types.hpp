/** @file StateVector_state_types.hpp

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

#ifndef APEGRUNT_STATEVECTOR_STATE_TYPES_HPP
#define APEGRUNT_STATEVECTOR_STATE_TYPES_HPP

#include <cctype>

#include <boost/serialization/serialization.hpp>
#include <boost/serialization/nvp.hpp>

#include <boost/spirit/include/qi_symbols.hpp>

namespace apegrunt {
// /*
template< typename StateT >
struct state_name { static const std::string value; };

template< typename StateT >
const std::string state_name<StateT>::value = "generic";

template< typename StateT >
StateT char_to_state( char symbol )
{
	std::cerr << "Call to dummy implementation of template function \"char_to_state( char symbol )\". Please implement the proper specialization." << std::endl;
	return {};
}

template< typename StateT >
char state_to_char( StateT state )
{
	std::cerr << "Call to dummy implementation of template function \"state_to_char( StateT state )\". Please implement the proper specialization." << std::endl;
	return '-';
}
// */
template< typename InputState, typename OutputState >
OutputState to_state( InputState state_or_symbol )
{
	std::cerr << "Call to dummy implementation of template function \"to_state( nucleotide )\". Please implement the proper specialization." << std::endl;
	std::cerr << "with InputState=" << state_name<InputState>::value << ", OutputState=" << state_name<OutputState>::value << std::endl;
	return {};
}

template< typename StateT >
char to_char( StateT state )
{
	std::cerr << "Call to dummy implementation of template function \"to_char( StateT state )\". Please implement the proper specialization." << std::endl;
	return '-';
}

template< typename StateT >
struct has_gap_state : std::true_type { };

template< typename StateT >
struct gap_state { static constexpr StateT value = StateT::GAP; };

template< typename StateT >
struct number_of_states { enum { N=0, value=0 }; };

//template< typename StateT >
//struct state_array;

// Nucleic acid states

enum struct nucleic_acid_state_t : uint8_t
{
	a=0,
	t=1,
	c=2,
	g=3,
	GAP=4
};

void swap( nucleic_acid_state_t& a, nucleic_acid_state_t& b ); // { nucleic_acid_state_t temp(b); b=a; a=temp; }

//template<> const std::string state_name<nucleic_acid_state_t>::value = "nucleic_acid_state_t";
template<> struct number_of_states<nucleic_acid_state_t> { enum { N=5, value=5 }; };
template<> nucleic_acid_state_t char_to_state<nucleic_acid_state_t>( char nucleotide );
template<> nucleic_acid_state_t to_state<char,nucleic_acid_state_t>( char nucleotide );
template<> nucleic_acid_state_t to_state<nucleic_acid_state_t,nucleic_acid_state_t>( nucleic_acid_state_t nucleotide );
template<> char state_to_char<nucleic_acid_state_t>( nucleic_acid_state_t state );
template<> char to_char<nucleic_acid_state_t>( nucleic_acid_state_t state );

/*
template<>
const std::string state_name<nucleic_acid_state_t>::value = "nucleic_acid_state_t";

template<> struct number_of_states<nucleic_acid_state_t> { enum { N=5, value=5 }; };

template<>
nucleic_acid_state_t char_to_state<nucleic_acid_state_t>( char nucleotide )
{
	using state_t = nucleic_acid_state_t;
	switch( std::tolower(nucleotide) )
	{
		case 'a' : return state_t::a; break;
		case 't' : return state_t::t; break;
		case 'c' : return state_t::c; break;
		case 'g' : return state_t::g; break;
		default : return state_t::GAP;
	}
}

template<>
nucleic_acid_state_t to_state<char,nucleic_acid_state_t>( char nucleotide )
{
	return char_to_state<nucleic_acid_state_t>(nucleotide);
}

template<>
nucleic_acid_state_t to_state<nucleic_acid_state_t,nucleic_acid_state_t>( nucleic_acid_state_t nucleotide )
{
	return nucleotide;
}

//template< >
//struct state_array<nucleic_acid_state_t>
//{
//	std::array< nucleic_acid_state_t, number_of_states<nucleic_acid_state_t>::value > states{{ nucleic_acid_state_t::a, nucleic_acid_state_t::t, nucleic_acid_state_t::c, nucleic_acid_state_t::g, nucleic_acid_state_t::GAP }};
//};

template<>
constexpr char state_to_char<nucleic_acid_state_t>( nucleic_acid_state_t state )
{
	using state_t = nucleic_acid_state_t;
	switch( state )
	{
		case state_t::a : return 'a'; break;
		case state_t::t : return 't'; break;
		case state_t::c : return 'c'; break;
		case state_t::g : return 'g'; break;
		default : return '-';
	}
}

template<>
constexpr char to_char<nucleic_acid_state_t>( nucleic_acid_state_t state )
{
	return state_to_char<nucleic_acid_state_t>(state);
}
*/

// Using chars as states

template<> struct number_of_states<char> { enum { N=256, value=256 }; };
/*
template<>
constexpr char char_to_state<char>( char nucleotide ) {	return nucleotide; }

template<>
constexpr char state_to_char<char>( char state ) { return state; }
*/
template<>
constexpr char to_state<char,char>( char nucleotide ) {	return nucleotide; }

template<>
constexpr char to_char<char>( char state ) { return state; }

template<>
struct gap_state<char> { static constexpr char value = '-'; };

template<>
struct gap_state<unsigned char> { static constexpr unsigned char value = '-'; };

// a wrapper/holder for the raw state type
template< typename StateT >
struct State_holder
{
	using state_t = StateT;
	using my_type = State_holder<state_t>;

	State_holder() = delete;
	//State_holder() : m_state(gap_state<state_t>::value) { }
	~State_holder() = default;

	//template< typename InputState >
    template< typename InputState, typename = std::enable_if_t<!std::is_same<InputState, state_t>::value> >
	State_holder( InputState state ) : m_state(to_state<InputState,StateT>(state)) { }
	State_holder( StateT state ) : m_state(state) { }
	//State_holder( char state ) : m_state( char_to_state<StateT>(state) ) { }

	//template< typename InputState >
    template< typename InputState, typename = std::enable_if_t<!std::is_same<InputState, state_t>::value> >
	inline my_type& operator=( InputState state ) { m_state = to_state<InputState,StateT>(state); return *this; }
	inline my_type& operator=( StateT state ) { m_state = state; return *this; }
	//inline my_type& operator=( char state ) { m_state = char_to_state<StateT>(state); return *this; }

    inline bool operator==( const my_type& rhs ) const { return m_state == rhs.m_state; }
    inline bool operator<( const my_type& rhs ) const { return m_state < rhs.m_state; }

    inline operator char() const { return to_char<StateT>( m_state ); }
    //inline operator char() const { return state_to_char<StateT>( m_state ); }

    template< typename Q = state_t, typename = std::enable_if_t<!std::is_same<Q, char>::value> >
    inline operator state_t() const { return m_state; }

    template< typename IntegerT, typename = typename std::enable_if_t< std::is_integral<IntegerT>::value && !std::is_same<IntegerT, char>::value > >
	inline operator IntegerT() const { return IntegerT( m_state ); }

    StateT m_state;

	static constexpr state_t get_gap() { return state_t::GAP; }

	friend class boost::serialization::access;
	/// boost.serialization interface.
	template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
		ar & BOOST_SERIALIZATION_NVP(m_state);
    }
};

template< typename StateT >
std::ostream& operator<< ( std::ostream& os, const State_holder<StateT>& state )
{
	os << char(state);
	return os;
}

template< typename StateT >
inline bool operator!=( const State_holder<StateT>& lhs, const State_holder<StateT>& rhs ) { return !(lhs == rhs); }
template< typename StateT >
inline bool operator>( const State_holder<StateT>& lhs, const State_holder<StateT>& rhs ) { return (rhs < lhs); }
template< typename StateT >
inline bool operator<=( const State_holder<StateT>& lhs, const State_holder<StateT>& rhs ) { return !(lhs > rhs); }
template< typename StateT >
inline bool operator>=( const State_holder<StateT>& lhs, const State_holder<StateT>& rhs ) { return !(lhs < rhs); }

template< typename StateT >
inline bool operator==( const State_holder<StateT>& lhs, const StateT& rhs ) { return (lhs == State_holder<StateT>(rhs) ); }
template< typename StateT >
inline bool operator<( const State_holder<StateT>& lhs, const StateT& rhs ) { return (lhs < State_holder<StateT>(rhs) ); }

template< typename StateT >
inline bool operator!=( const State_holder<StateT>& lhs, const StateT& rhs ) { return !(lhs == State_holder<StateT>(rhs) ); }
template< typename StateT >
inline bool operator>( const State_holder<StateT>& lhs, const StateT& rhs ) { return (lhs > State_holder<StateT>(rhs) ); }
template< typename StateT >
inline bool operator<=( const State_holder<StateT>& lhs, const StateT& rhs ) { return !(lhs > State_holder<StateT>(rhs) ); }
template< typename StateT >
inline bool operator>=( const State_holder<StateT>& lhs, const StateT& rhs ) { return !(lhs < State_holder<StateT>(rhs) ); }

template< typename StateT >
struct gap_state< State_holder<StateT> > { static constexpr StateT value = StateT::GAP; };

template<>
struct gap_state< State_holder<char> > { static constexpr char value = gap_state<char>::value; };

template<>
struct gap_state< State_holder<unsigned char> > { static constexpr unsigned char value = gap_state<unsigned char>::value; };

} // namespace apegrunt

#endif // APEGRUNT_STATEVECTOR_STATE_TYPES_HPP

