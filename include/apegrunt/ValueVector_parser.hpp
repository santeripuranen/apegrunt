/** @file ValueVector_parser.hpp

	Copyright (c) 2018 Santeri Puranen.

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

#ifndef APEGRUNT_VALUE_VECTOR_PARSER_HPP
#define APEGRUNT_VALUE_VECTOR_PARSER_HPP

#include <iosfwd> // for passing std::istream*
#include <string>

#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/support_istream_iterator.hpp>
#include <boost/iostreams/device/mapped_file.hpp> // for boost::iostreams::mapped_file
#include <boost/filesystem/operations.hpp> // includes boost/filesystem/path.hpp

#include "ValueVector_forward.h"
#include "ValueVector_parser_grammar.hpp"

namespace apegrunt {

namespace fs = boost::filesystem;
namespace ascii	= boost::spirit::ascii;
namespace spirit = boost::spirit;

template< typename ValueT=double >
ValueVector_ptr<ValueT> parse_ValueVector( const std::string& filename )
{
	using value_t = ValueT;
	using iterator_t = char const*;

	ValueVector_ptr< value_t > weights;

	fs::path filepath( filename );

	if( fs::exists( filepath ) && fs::is_regular_file(filepath) )
	{
		boost::iostreams::mapped_file mmap(filename.c_str(), boost::iostreams::mapped_file::readonly);
		iterator_t begin = mmap.const_data();
		iterator_t end = begin + mmap.size();

		apegrunt::parsers::ValueVector_parser_grammar< iterator_t, value_t > parser;
		const bool success = boost::spirit::qi::phrase_parse(begin, end, parser, ascii::space, weights);

		mmap.close();
		if( !success )
		{
			return {};
		}
	}

	return weights;
}

} // namespace apegrunt

#endif // APEGRUNT_VALUE_VECTOR_PARSER_HPP

