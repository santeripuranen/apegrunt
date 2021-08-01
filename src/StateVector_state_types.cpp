/** @file StateVector_state_types.cpp

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

#include "apegrunt/StateVector_state_types.hpp"

namespace apegrunt {

void swap( nucleic_acid_state_t& a, nucleic_acid_state_t& b ) { nucleic_acid_state_t temp(b); b=a; a=temp; }

template<> const std::string state_name<nucleic_acid_state_t>::value = "nucleic_acid_state_t";

//template<> struct number_of_states<nucleic_acid_state_t> { enum { N=5, value=5 }; };

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

template<>
char state_to_char<nucleic_acid_state_t>( nucleic_acid_state_t state )
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
char to_char<nucleic_acid_state_t>( nucleic_acid_state_t state )
{
	return state_to_char<nucleic_acid_state_t>(state);
}

// Using chars as states

//template<> struct number_of_states<char> { enum { N=256, value=256 }; };
/*
template<>
constexpr char char_to_state<char>( char nucleotide ) {	return nucleotide; }

template<>
constexpr char state_to_char<char>( char state ) { return state; }

template<>
constexpr char to_state<char,char>( char symbol ) {	return symbol; }

//template<>
//constexpr char to_char<char>( char state ) { return state; }

//template<>
//struct gap_state<char> { static constexpr char value = '-'; };

template<>
struct gap_state<unsigned char> { static constexpr unsigned char value = '-'; };

template<>
struct has_gap_state< State_holder<binary_state_t> > : std::false_type { };

template<> struct state_name< State_holder<binary_state_t> > { static constexpr auto value = "State_holder<binary_state_t>"; };

template<>
struct gap_state< State_holder<char> > { static constexpr char value = gap_state<char>::value; };

template<>
struct gap_state< State_holder<unsigned char> > { static constexpr unsigned char value = gap_state<unsigned char>::value; };
*/
} // namespace apegrunt
