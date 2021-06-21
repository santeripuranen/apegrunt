/** @file Loci_parsers.hpp

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

#ifndef APEGRUNT_LOCI_PARSERS_HPP
#define APEGRUNT_LOCI_PARSERS_HPP

#include <iosfwd> // for passing std::istream*
#include <string>
#include <random> // C++11

#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/support_istream_iterator.hpp>
#include <boost/iostreams/device/mapped_file.hpp> // for boost::iostreams::mapped_file
#include <boost/filesystem/operations.hpp> // includes boost/filesystem/path.hpp

//#include "Loci.h"
#include "Loci_impl_default_storage.hpp"
#include "Loci_parser_grammar.hpp"

namespace apegrunt {

namespace fs = boost::filesystem;
namespace ascii	= boost::spirit::ascii;
namespace spirit = boost::spirit;

template<typename LociT=Loci_impl_default_storage<std::size_t> >
Loci_ptr combine( Loci_ptr lhs, Loci_ptr rhs )
{
/*
	std::cout << "combine Loci_ptrs \"" << lhs->id_string() << "\" and \"" << rhs->id_string() << "\"" << std::endl;
	std::cout << "  lhs = {" << lhs << " }" << std::endl;
	std::cout << "  rhs = {" << rhs << " }" << std::endl;
*/
	std::vector<std::size_t> combined; combined.reserve( rhs->size() );
	for( auto index: rhs )
	{
		combined.push_back( (*lhs)[index] );
	}
	return make_Loci_ptr<LociT>( std::move( combined ) );
}

template< typename LociT=Loci_impl_default_storage<std::size_t> >
Loci_ptr parse_Loci_list( const std::string& infilename, std::size_t base_index=0 )
{
	using iterator_t = char const*;

	Loci_ptr loci;

	fs::path filepath( infilename );

	if( fs::exists( filepath ) && fs::is_regular_file(filepath) )
	{
		boost::iostreams::mapped_file mmap(infilename.c_str(), boost::iostreams::mapped_file::readonly);
		iterator_t begin = mmap.const_data();
		iterator_t end = begin + mmap.size();

		apegrunt::parsers::Loci_parser_grammar< iterator_t, LociT > parser( base_index ); // 0==zero-based numbering, 1==one-based numbering in the input
		const bool success = boost::spirit::qi::phrase_parse(begin, end, parser, ascii::space, loci);
		if( success )
		{
			loci->set_id_string( filepath.stem().c_str() );
		}

		mmap.close();
	}
	return loci;
}

template< typename LociT=Loci_impl_default_storage<std::size_t> >
Loci_ptr make_Loci_list( std::ostringstream& loci_stream, std::size_t base_index=0 )
{
	using iterator_t = char const*;
	Loci_ptr loci;

	auto loci_string = loci_stream.str();

	iterator_t begin = loci_string.data();
	iterator_t end = begin + loci_string.size();

	apegrunt::parsers::Loci_parser_grammar< iterator_t, LociT > parser( base_index ); // 0==zero-based numbering, 1==one-based numbering in the input
	const bool success = boost::spirit::qi::phrase_parse(begin, end, parser, ascii::space, loci);
	if( success )
	{
		//loci->set_id_string( std::string("default") );
	}
	return loci;
}

template< typename LociT=Loci_impl_default_storage<std::size_t> >
Loci_ptr make_Loci_list( std::vector<std::size_t>&& loci_list, std::size_t base_index=0 )
{
	return make_Loci_ptr<LociT>(std::move(loci_list));
}
/*
Loci_ptr create_random_Loci_list(
		std::size_t n_loci,
		std::size_t sample_count, // number of samples
		std::size_t random_seed=0 // if random_seed==0, then generate a new random random seed )
	)
{
	// initialize random number generator
	if( 0 == random_seed )
	{
		std::random_device rd; // used only once -- generates random seed for mt generator below
		random_seed = rd();
	}
	std::mt19937_64 generator(random_seed);
	std::uniform_int_distribution<std::size_t> distribution( 0, n_loci-1 );
	auto grand = std::bind( distribution, generator );

	auto accept_list = std::vector<std::size_t>();
	accept_list.reserve( sample_count );

	for( std::size_t i=0; i < sample_count; ++i ) { accept_list.push_back( grand() ); }

	return make_Loci_list( std::move(accept_list), 0 );
}
*/
} // namespace apegrunt

#endif // APEGRUNT_LOCI_PARSERS_HPP

