/** @file Valeu_list_parser_grammar.hpp

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

#ifndef APEGRUNT_VALUE_LIST_PARSER_GRAMMAR_HPP
#define APEGRUNT_VALUE_LIST_PARSER_GRAMMAR_HPP

#define BOOST_SPIRIT_USE_PHOENIX_V3 1
//#define BOOST_RESULT_OF_USE_DECLTYPE
//#define BOOST_RESULT_OF_USE_DECLTYPE_WITH_TR1_FALLBACK 1

#include <iosfwd>
#include <memory> // for std::unique_ptr and std::shared_ptr
#include <string>

#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/qi_symbols.hpp>
#include <boost/spirit/include/phoenix.hpp>

#include "Parser_commons.hpp"
#include "ValueVector_forward.h"

namespace apegrunt {

namespace parsers {

namespace qi	= boost::spirit::qi;
namespace ascii	= boost::spirit::ascii;
namespace phx	= boost::phoenix;

template< typename IteratorT, typename ValueT >
struct ValueVector_parser_grammar
	: qi::grammar< IteratorT, ValueVector_ptr<ValueT>(), qi::locals< std::shared_ptr<ValueVector<ValueT> > >, ascii::space_type >
{
	using space_t = ascii::space_type;

	using value_t = ValueT;
	using value_vector_t = ValueVector<value_t>;
	using value_vector_ptr_t = std::shared_ptr<value_vector_t>;

	// Explicit member function pointer signatures; used for the static cast below, in order for function
	// pointer signatures to be properly resolved when compiling with gcc against boost >1.62.0
	typedef void (value_vector_t::*push_back)(const value_t& val);

	explicit ValueVector_parser_grammar( std::size_t base_index=1 )
		: ValueVector_parser_grammar::base_type(start, "ValueVector_parser_grammar")
	{
		using apegrunt::parsers::printer;

		single_entry = (qi::double_)[ phx::bind( static_cast<push_back>(&value_vector_t::push_back), qi::_r1, qi::_1 ) ]
		;

		// Main node
		start =	/* a simple "qi::_a = std::make_shared<T>()" would not be lazy,
				instead it's evaluated immediately and only once. lazy_make_shared evaluates
				lazily each time the parent parser is successfully called. */
				qi::eps[ qi::_a = apegrunt::parsers::lazy_make_shared<value_vector_t>()() ]
			 >> *( single_entry( phx::bind( &value_vector_ptr_t::get, qi::_a ) ) ) // pass a raw ptr
			 >> qi::eps[ qi::_val = qi::_a ]
		;
	}

	qi::rule<IteratorT, void(value_vector_t* const), space_t>
		single_entry
	;

	qi::rule< IteratorT, ValueVector_ptr<value_t>(), qi::locals< std::shared_ptr<value_vector_t> >, ascii::space_type > start;
};

} // namespace parsers

} // namespace apegrunt

#endif // APEGRUNT_VALUE_LIST_PARSER_GRAMMAR_HPP

