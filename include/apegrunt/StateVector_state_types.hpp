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
	std::cerr << "Call to dummy implementation of template function \"char_to_state( char nucleotide )\". Please implement the proper specialization." << std::endl;
	return StateT();
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
	return OutputState();
}

template< typename StateT >
char to_char( StateT state )
{
	std::cerr << "Call to dummy implementation of template function \"to_char( StateT state )\". Please implement the proper specialization." << std::endl;
	return '-';
}

template< typename StateT >
struct gap_state { static constexpr StateT value = StateT::GAP; };

template< typename StateT >
struct number_of_states { enum { N=0, value=0 }; };

template< typename StateT >
struct state_array;

// Nucleic acid states

enum struct nucleic_acid_state_t : uint8_t
{
	a=0,
	t=1,
	c=2,
	g=3,
	GAP=4
};

//template<>
//struct gap_state<nucleic_acid_state_t> { static constexpr nucleic_acid_state_t value = nucleic_acid_state_t::GAP; };

template<>
const std::string state_name<nucleic_acid_state_t>::value = "nucleic_acid_state_t";

void swap( nucleic_acid_state_t& a, nucleic_acid_state_t& b ) { nucleic_acid_state_t temp(b); b=a; a=temp; }

/*
static boost::spirit::qi::symbols< nucleic_acid_state_t, nucleic_acid_state_t > nucleic_acid_base_complement;
nucleic_acid_base_complement.add
	( nucleic_acid_state_t::a, nucleic_acid_state_t::t )
	( nucleic_acid_state_t::t, nucleic_acid_state_t::a )
	( nucleic_acid_state_t::c, nucleic_acid_state_t::g )
	( nucleic_acid_state_t::g, nucleic_acid_state_t::c )
;
*/
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

template< >
struct state_array<nucleic_acid_state_t>
{
	std::array< nucleic_acid_state_t, number_of_states<nucleic_acid_state_t>::value > states{{ nucleic_acid_state_t::a, nucleic_acid_state_t::t, nucleic_acid_state_t::c, nucleic_acid_state_t::g, nucleic_acid_state_t::GAP }};
};

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


// Major/minor nucleic acid states
/*
enum struct biallelic_state_t : uint8_t
{
	major_allele=0, // major allele
	minor_allele=1, // minor allele
	GAP=2
};

template<> struct number_of_states<biallelic_state_t> { enum { N=3, value=3 }; };

template<>
biallelic_state_t char_to_state<biallelic_state_t>( char nucleotide )
{
	using state_t = biallelic_state_t;
	switch( std::tolower(nucleotide) )
	{
		case '0' : return state_t::major_allele; break;
		case '1' : return state_t::minor_allele; break;
		default : return state_t::GAP;
	}
}

template<>
char state_to_char<biallelic_state_t>( biallelic_state_t state )
{
	using state_t = biallelic_state_t;
	switch( state )
	{
		case state_t::major_allele : return '0'; break;
		case state_t::minor_allele : return '1'; break;
		default : return '-';
	}
}
*/
// Major/mid/minor nucleic acid states

enum struct triallelic_state_t : uint8_t
{
	major_allele=0, // major allele
	mid_allele=1, // mid allele
	minor_allele=2, // minor allele
	GAP=3
};

//template<>
//struct gap_state<triallelic_state_t> { static constexpr triallelic_state_t value = triallelic_state_t::GAP; };

template<>
const std::string state_name<triallelic_state_t>::value = "triallelic_state_t";

template<> struct number_of_states<triallelic_state_t> { enum { N=4, value=4 }; };

template<>
triallelic_state_t char_to_state<triallelic_state_t>( char nucleotide )
{
	using state_t = triallelic_state_t;
	switch( std::tolower(nucleotide) )
	{
		case '0' : return state_t::major_allele; break;
		case '1' : return state_t::mid_allele; break;
		case '2' : return state_t::minor_allele; break;
		default : return state_t::GAP;
	}
}

template<>
triallelic_state_t to_state<nucleic_acid_state_t,triallelic_state_t>( nucleic_acid_state_t symbol )
{
	return char_to_state<triallelic_state_t>(to_char(symbol));
}

template<>
triallelic_state_t to_state<char,triallelic_state_t>( char symbol )
{
	return char_to_state<triallelic_state_t>(symbol);
}

template<>
constexpr char state_to_char<triallelic_state_t>( triallelic_state_t state )
{
	using state_t = triallelic_state_t;
	switch( state )
	{
		case state_t::major_allele : return '0'; break;
		case state_t::mid_allele : return '1'; break;
		case state_t::minor_allele : return '2'; break;
		default : return '-';
	}
}

template<>
constexpr char to_char<triallelic_state_t>( triallelic_state_t state )
{
	return state_to_char<triallelic_state_t>(state);
}

// Tetra-allelic states type

enum struct tetraallelic_state_t : uint8_t
{
	f0=0, // most frequent allele
	f1=1, // 2nd most frequent allele
	f2=2, // 3rd most frequent allele
	f3=3 // 4th most frequent allele
};

template<> struct number_of_states<tetraallelic_state_t> { enum { N=4, value=4 }; };

template<>
tetraallelic_state_t char_to_state<tetraallelic_state_t>( char symbol )
{
	using state_t = tetraallelic_state_t;
	switch( std::tolower(symbol) )
	{
		case '0' : return state_t::f0; break;
		case '1' : return state_t::f1; break;
		case '2' : return state_t::f2; break;
		case '3' : return state_t::f3; break;
	}
}

template<>
tetraallelic_state_t to_state<char,tetraallelic_state_t>( char symbol )
{
	return char_to_state<tetraallelic_state_t>(symbol);
}

template< >
struct state_array<tetraallelic_state_t>
{
	using state_t = tetraallelic_state_t;
	std::array< state_t, number_of_states<state_t>::value > states{{ state_t::f0, state_t::f1, state_t::f2, state_t::f3 }};
};

template<>
constexpr char state_to_char<tetraallelic_state_t>( tetraallelic_state_t state )
{
	using state_t = tetraallelic_state_t;
	switch( state )
	{
		case state_t::f0 : return '0'; break;
		case state_t::f1 : return '1'; break;
		case state_t::f2 : return '2'; break;
		case state_t::f3 : return '3'; break;
		default : return '-';
	}
}

template<>
constexpr char to_char<tetraallelic_state_t>( tetraallelic_state_t state )
{
	return state_to_char<tetraallelic_state_t>(state);
}

// Penta-allelic states type

enum struct pentaallelic_state_t : uint8_t
{
	f0=0, // most frequent allele
	f1=1, // 2nd most frequent allele
	f2=2, // 3rd most frequent allele
	f3=3, // 4th most frequent allele
	f4=4 // 5th most frequent allele
};

template<> struct number_of_states<pentaallelic_state_t> { enum { N=5, value=5 }; };

template<>
pentaallelic_state_t char_to_state<pentaallelic_state_t>( char symbol )
{
	using state_t = pentaallelic_state_t;
	switch( std::tolower(symbol) )
	{
		case '0' : return state_t::f0; break;
		case '1' : return state_t::f1; break;
		case '2' : return state_t::f2; break;
		case '3' : return state_t::f3; break;
		case '4' : return state_t::f4; break;
	}
}

template<>
pentaallelic_state_t to_state<char,pentaallelic_state_t>( char symbol )
{
	return char_to_state<pentaallelic_state_t>(symbol);
}

template< >
struct state_array<pentaallelic_state_t>
{
	using state_t = pentaallelic_state_t;
	std::array< state_t, number_of_states<state_t>::value > states{{ state_t::f0, state_t::f1, state_t::f2, state_t::f3, state_t::f4 }};
};

template<>
constexpr char state_to_char<pentaallelic_state_t>( pentaallelic_state_t state )
{
	using state_t = pentaallelic_state_t;
	switch( state )
	{
		case state_t::f0 : return '0'; break;
		case state_t::f1 : return '1'; break;
		case state_t::f2 : return '2'; break;
		case state_t::f3 : return '3'; break;
		case state_t::f4 : return '4'; break;
		default : return '-';
	}
}

template<>
constexpr char to_char<pentaallelic_state_t>( pentaallelic_state_t state )
{
	return state_to_char<pentaallelic_state_t>(state);
}


/*
// Amino acid states

enum struct amino_acid_state_t : uint8_t
{
	A=0, C=1, D=2,
	E=3, F=4, G=5,
	H=6, I=7, K=8,
	L=9, M=10, N=11,
	P=12, Q=13, R=14,
	S=15, T=16, V=17,
	W=18, Y=19, GAP=20
};

template<> struct number_of_states<amino_acid_state_t> { enum { N=21, value=21 }; };

template<>
amino_acid_state_t char_to_state<amino_acid_state_t>( char aminoacid )
{
	using state_t = amino_acid_state_t;
	switch( std::toupper(aminoacid) )
	{
		case 'A' : return state_t::A; break;
		case 'C' : return state_t::C; break;
		case 'D' : return state_t::D; break;
		case 'E' : return state_t::E; break;
		case 'F' : return state_t::F; break;
		case 'G' : return state_t::G; break;
		case 'H' : return state_t::H; break;
		case 'I' : return state_t::I; break;
		case 'K' : return state_t::K; break;
		case 'L' : return state_t::L; break;
		case 'M' : return state_t::M; break;
		case 'N' : return state_t::N; break;
		case 'P' : return state_t::P; break;
		case 'Q' : return state_t::Q; break;
		case 'R' : return state_t::R; break;
		case 'S' : return state_t::S; break;
		case 'T' : return state_t::T; break;
		case 'V' : return state_t::V; break;
		case 'W' : return state_t::W; break;
		case 'Y' : return state_t::Y; break;
		default : return state_t::GAP;
	}
}

template<>
char state_to_char<amino_acid_state_t>( amino_acid_state_t state )
{
	using state_t = amino_acid_state_t;
	switch( state )
	{
		case state_t::A : return 'A'; break;
		case state_t::C : return 'C'; break;
		case state_t::D : return 'D'; break;
		case state_t::E : return 'E'; break;
		case state_t::F : return 'F'; break;
		case state_t::G : return 'G'; break;
		case state_t::H : return 'H'; break;
		case state_t::I : return 'I'; break;
		case state_t::K : return 'K'; break;
		case state_t::L : return 'L'; break;
		case state_t::M : return 'M'; break;
		case state_t::N : return 'N'; break;
		case state_t::P : return 'P'; break;
		case state_t::Q : return 'Q'; break;
		case state_t::R : return 'R'; break;
		case state_t::S : return 'S'; break;
		case state_t::T : return 'T'; break;
		case state_t::V : return 'V'; break;
		case state_t::W : return 'W'; break;
		case state_t::Y : return 'Y'; break;
		default : return '-';
	}
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

	State_holder() : m_state(gap_state<state_t>::value) { }
	~State_holder() { } //= default;

	template< typename InputState >
	State_holder( InputState state ) : m_state(to_state<InputState,StateT>(state)) {}
	State_holder( StateT state ) : m_state(state) { }
	//State_holder( char state ) : m_state( char_to_state<StateT>(state) ) { }

	template< typename InputState >
	inline my_type& operator=( InputState state ) { m_state = to_state<InputState,StateT>(state); return *this; }
	inline my_type& operator=( StateT state ) { m_state = state; return *this; }
	//inline my_type& operator=( char state ) { m_state = char_to_state<StateT>(state); return *this; }

    inline bool operator==( const my_type& rhs ) const { return m_state == rhs.m_state; }
    inline bool operator<( const my_type& rhs ) const { return m_state < rhs.m_state; }

    inline operator char() const { return to_char<StateT>( m_state ); }
    //inline operator char() const { return state_to_char<StateT>( m_state ); }

    inline operator std::size_t() const { return std::size_t(m_state); }
    inline operator uint8_t() const { return uint8_t(m_state); }
	//inline StateT operator()() const { return m_state; }
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

