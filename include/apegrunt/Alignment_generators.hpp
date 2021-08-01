/** @file Alignment_generators.hpp

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

#ifndef APEGRUNT_ALIGNMENT_GENERATORS_HPP
#define APEGRUNT_ALIGNMENT_GENERATORS_HPP

#include <iosfwd> // for passing std::istream*
#include <algorithm> // for std::copy

//#include <boost/spirit/include/karma.hpp>
//#include <boost/spirit/include/support_istream_iterator.hpp>
//#include <boost/iostreams/device/mapped_file.hpp> // for boost::iostreams::mapped_file
//#include <boost/filesystem/operations.hpp> // includes boost/filesystem/path.hpp

#include "Alignment_forward.h"
#include "Alignment_interface.hpp"
//#include "Alignment_generator_FASTA_grammar.hpp"

namespace apegrunt {

namespace fs = boost::filesystem;
namespace ascii	= boost::spirit::ascii;
namespace spirit = boost::spirit;

template< typename StateT >
bool generate_Alignment( Alignment_ptr<StateT> alignment, std::ostream* outstream )
{
// this is temporary file output code, used until we have a proper generator implementation.
	if( outstream && outstream->good() )
	{
		for( const auto& sequence: alignment )
		{
			*outstream << ">" << sequence->id_string() << "\n";
			for( const auto& state: sequence ) { *outstream << to_char(state); }
			*outstream << "\n";
		}
		return true;
	}
	return false;
}

template< typename StateT >
bool generate_PhandangoCSV( Alignment_ptr<StateT> alignment, std::ostream* outstream )
{
	if( outstream && outstream->good() )
	{
		auto index_translation = alignment->get_loci_translation();
		const std::size_t base_index = apegrunt::Apegrunt_options::get_output_indexing_base();

		// create header
		*outstream << "\"locus\"";
		for( const auto& index : alignment->get_loci_translation() ) { *outstream << ",\"L" << index+base_index << ":o1\""; }
		*outstream << ",\"L" << *begin(alignment->get_loci_translation())+base_index << ":o1\"";
		*outstream << "\n";

		// the payload
		for( const auto& sequence: alignment )
		{
			*outstream << "\"" << sequence->id_string() << "\"";
			for( const auto& state: sequence ) { *outstream << ",\"" << state << "\""; }
			*outstream << ",\"" << *(begin(sequence)) << "\"";
			*outstream << "\n";
		}
		return true;
	}
	return false;
}
/*
template< typename StateT >
bool generate_AG_Alignment( apegrunt::Alignment_ptr<StateT> alignment, std::ostream* outstream )
{
// this is temporary file output code, used until we have a proper generator implementation.
	if( outstream && outstream->good() )
	{
		*outstream << "<alignment id=\"" << alignment->id_string() << "\">" << "\n";
		*outstream << "<instances>";
		std::hash<std::string> seqid_hasher;
		std::size_t sid=0;
		for( const auto& sequence: alignment )
		{
			*outstream << "<i id=\"" << sid << "\" w=\"" << sequence->weight() << "\">"
					<< sequence->id_string()
					<< "</i>"
			;
			++sid;
		}
		*outstream << "</samples>" << "\n";

		//const auto& alignment = *(alignment.get());
		const auto& blocks = *(alignment->get_block_storage());
		const auto& block_accounting = *(alignment->get_block_accounting());
		const auto& index_translation = *(alignment->get_loci_translation());
		const std::size_t base_index = apegrunt::Apegrunt_options::get_output_indexing_base();

		const auto n_loci = alignment->n_loci();
		const auto number_of_blocks = apegrunt::get_number_of_blocks(n_loci);
		const auto last_block_size = apegrunt::get_last_block_size(n_loci);
		const auto last_block_index = apegrunt::get_last_block_index(n_loci);

		*outstream << "<data>";
		for( std::size_t i=0; i < number_of_blocks; ++i )
		{
			const auto sidxend = i == last_block_index ? last_block_size : apegrunt::StateBlock_size;
			{ // the block 'header'
				*outstream << "<blocks p=\""; // p=positions
				std::size_t oldidx = 0, begidx = 0;
				*outstream << std::hex << index_translation[i*apegrunt::StateBlock_size]+base_index;
				for( std::size_t sidx=1; sidx < sidxend; ++sidx )
				{
					if( index_translation[i*apegrunt::StateBlock_size+sidx]+base_index > index_translation[i*apegrunt::StateBlock_size+oldidx]+base_index+1 )
					{
						if( oldidx != begidx )
						{
							*outstream << "-" << index_translation[i*apegrunt::StateBlock_size+oldidx]+base_index;
						}
						*outstream << "," << index_translation[i*apegrunt::StateBlock_size+sidx]+base_index;
						begidx = sidx;
					}
					oldidx = sidx;
				}
				*outstream << std::dec << "\">";
			}
			{ // the block 'payload'
				for( std::size_t nblock=0; nblock < block_accounting[i].size(); ++nblock )
				{
					*outstream << "<block i=\""; // i=instances
					const auto& indices = block_accounting[i][nblock];
					std::size_t oldpos = 0, begpos = 0;
					*outstream << std::hex << indices[oldpos];
					for( std::size_t pos=1; pos < indices.size(); ++pos )
					{
						if( indices[pos] > indices[oldpos]+1 )
						{
							if( oldpos != begpos )
							{
								*outstream << "-" << indices[oldpos];
							}
							*outstream << "," << indices[pos];
							begpos = pos;
						}
						oldpos = pos;
					}
					*outstream << std::dec << "\">";
					*outstream << blocks[i][nblock] << "</block>";
				}
				*outstream << "</blocks>" << "\n";
			}
		}
		*outstream << "</data>" << "\n";
		*outstream << "</alignment>" << "\n";
		return true;
	}
	return false;
}
*/
/*
// Boost Spirit Karma generator
template< typename StateT >
bool generate_Alignment( Alignment_ptr<StateT> alignment, std::ostream* outstream )
{
	using iterator_t = spirit::ostream_iterator;
	bool success = false;

	if( outstream && outstream->good() )
	{
		// wrap ostream into iterator
		iterator_t iterator( *outstream );
		//iterator_t end;

		apegrunt::generators::Alignment_generator_FASTA_grammar< iterator_t, StateT > generator;
		success = boost::spirit::karma::generate(iterator, generator, alignment);
	}
	return success;
}
*/
/*
// Boost Spirit Karma generator
template< typename StateT >
bool generate_Alignment( Alignment_ptr<StateT> alignment, const std::string& outfilename )
{
	using iterator_t = char const*;
	bool success = false;

	fs::path filepath( outfilename );

	if( !fs::exists( filepath ) ) //&& fs::is_regular_file(filepath) )
	{
		boost::iostreams::mapped_file mmap(outfilename.c_str(), boost::iostreams::mapped_file::readwrite);
		iterator_t iterator = mmap.data();

		apegrunt::generators::Alignment_generator_FASTA_grammar< iterator_t, StateT > generator;
		success = boost::spirit::karma::generate(iterator, generator, alignment);
	}

	return success;
}
*/

} // namespace apegrunt

#endif // APEGRUNT_ALIGNMENT_GENERATORS_HPP

