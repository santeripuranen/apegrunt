/** @file Alignment_parsers.hpp

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

#ifndef APEGRUNT_ALIGNMENT_PARSERS_HPP
#define APEGRUNT_ALIGNMENT_PARSERS_HPP

#include <iosfwd> // for passing std::istream*
// #include <iterator> // for std::istream_iterator

#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/support_istream_iterator.hpp>
//#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/device/mapped_file.hpp> // for boost::iostreams::mapped_file
#include <boost/iostreams/stream.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filter/lzma.hpp>
#include <boost/filesystem/operations.hpp> // includes boost/filesystem/path.hpp

#include "Alignment.h"
#include "Alignment_parser_FASTA_grammar.hpp"

#include "Apegrunt_options.h"

namespace apegrunt {

namespace io = boost::iostreams;
namespace fs = boost::filesystem;
namespace ascii	= boost::spirit::ascii;
namespace spirit = boost::spirit;

template< typename AlignmentT >
Alignment_ptr< typename AlignmentT::state_t > parse_Alignment( std::istream* instream )
{
	using iterator_t = spirit::istream_iterator;
	//using iterator_t = std::istream_iterator<char>;
	Alignment_ptr< typename AlignmentT::state_t > alignment;
	if( instream && instream->good() )
	{
		instream->unsetf(std::ios::skipws); // No white space skipping!

		// wrap istream into iterator
		iterator_t first( *instream );
		iterator_t last;

		apegrunt::parsers::Alignment_parser_FASTA_grammar< iterator_t, AlignmentT > parser;
		bool success = boost::spirit::qi::phrase_parse(first, last, parser, ascii::space, alignment);
	}
	alignment->set_n_original_positions( Apegrunt_options::genome_size() != 0 ? Apegrunt_options::genome_size() : alignment->n_loci() );
	return alignment;
}

template< typename AlignmentT >
Alignment_ptr< typename AlignmentT::state_t > parse_Alignment( const std::string& infilename )
{
	using timer_t = std::chrono::steady_clock;

	Alignment_ptr< typename AlignmentT::state_t > alignment;

	fs::path filepath( infilename );
	if( fs::exists( filepath ) && fs::is_regular_file(filepath) )
	{
		auto start = timer_t::now();
		bool success(false);
		std::size_t bytes_read(0);

		// gzip or xz input; it works, but is *excruciatingly* slow
		if( filepath.extension() == ".gz" || filepath.extension() == ".xz" )
		{
			// phrase_parse wants forward_traversal_tag, but std::istream_iterator
			// only complies with single_pass_traversal_tag, so this won't compile.
			//using iterator_t = std::istream_iterator<std::size_t>;
			// use the boost::spirit equivalent instead
			using iterator_t = spirit::istream_iterator;

			io::mapped_file_source mmap(infilename.c_str());
			io::stream< io::mapped_file_source > source(mmap);

			// this would also do the job; it's just as fast as mapping in this case
			//io::file_source source( infilename.c_str() );

			// filtering_istream & istream_iterator are apparently a major source of overhead
			io::filtering_istream fs;
			if( filepath.extension() == ".gz" )	{ fs.push( io::gzip_decompressor{} ); }
			else { fs.push( io::lzma_decompressor{} ); }
			fs.push( source );

			iterator_t first(fs >> std::noskipws);
			iterator_t last;

			apegrunt::parsers::Alignment_parser_FASTA_grammar< iterator_t, AlignmentT > parser;
			success = boost::spirit::qi::phrase_parse(first, last, parser, ascii::space, alignment);

			bytes_read = mmap.size();
			mmap.close();
		}
		else // regular uncompressed input
		{
			using iterator_t = char const*;

			io::mapped_file mmap(infilename.c_str(), io::mapped_file::readonly);
			iterator_t first = mmap.const_data();
			iterator_t last = first + mmap.size();

			apegrunt::parsers::Alignment_parser_FASTA_grammar< iterator_t, AlignmentT > parser;
			success = boost::spirit::qi::phrase_parse(first, last, parser, ascii::space, alignment);

			bytes_read = mmap.size();
			mmap.close();
		}

		if( success )
		{
			alignment->set_id_string( filepath.stem().c_str() );
			alignment->set_n_original_positions( Apegrunt_options::genome_size() != 0 ? Apegrunt_options::genome_size() : alignment->n_loci() );

			if( apegrunt::Apegrunt_options::verbose() )
			{
				auto oldprecision = apegrunt::Apegrunt_options::get_out_stream()->precision();
				*apegrunt::Apegrunt_options::get_out_stream() << "apegrunt: parser data rate = " << std::fixed << std::setprecision(2) << (double(bytes_read)/double(std::size_t(1)<<20))/double(std::chrono::duration_cast<std::chrono::seconds>(timer_t::now()-start).count()) << " MiB/s" << std::endl;
				apegrunt::Apegrunt_options::get_out_stream()->precision(oldprecision);
			}
		}
	}
	return alignment;
}

} // namespace apegrunt

#endif // APEGRUNT_ALIGNMENT_PARSERS_HPP

