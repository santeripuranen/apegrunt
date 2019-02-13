/** @file Apegrunt.cpp
	Utility for manipulation and space-efficient storage of column-aligned categorical data.

	Copyright (c) 2017-2019 Santeri Puranen.

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

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <fstream>

#include "boost/program_options.hpp"
#include "boost/filesystem/operations.hpp" // includes boost/filesystem/path.hpp

#include "apegrunt/Apegrunt.h"
#include "apegrunt/Loci_generators.hpp"
#include "apegrunt/ValueVector_parser.hpp"
#include "misc/Stopwatch.hpp"

/*
 * The main program
 */
int main(int argc, char **argv)
{
	namespace po = boost::program_options;
	namespace fs = boost::filesystem;

	#ifndef APEGRUNT_NO_TBB // Threading with Threading Building Blocks
	tbb::task_scheduler_init tbb_task_scheduler(tbb::task_scheduler_init::deferred); // Threading task scheduler
	#endif // #ifndef APEGRUNT_NO_TBB

	using namespace apegrunt;
	using std::exit;

	using real_t = double;
	using default_state_t = apegrunt::nucleic_acid_state_t;
	using alignment_default_storage_t = apegrunt::Alignment_impl_block_compressed_storage< apegrunt::StateVector_impl_block_compressed_alignment_storage<default_state_t> >;

	// Parse command line options

	// Check the command line options
	po::options_description	all_options;
	apegrunt::Apegrunt_options apegrunt_options( &std::cout, &std::cerr );
	apegrunt_options.AddOptions( &all_options ); // add apegrunt library options

	po::variables_map options_map;

	// Catch the alignment file, even if not specified by a flag in the input
	po::positional_options_description popt;
	popt.add("alignmentfile", -1);

	std::ostringstream options_string;
	options_string << all_options;
	apegrunt_options.set_options_string( options_string.str() );

	try
	{
		po::store( po::command_line_parser(argc, argv).options(all_options).positional(popt).run(), options_map );
		po::notify(options_map);

		// Apegrunt options
		if( !apegrunt_options.CheckOptions(&options_map) )
		{
			exit(EXIT_FAILURE);
		}
		#ifndef APEGRUNT_NO_TBB // Threading with Threading Building Blocks
		Apegrunt_options::set_threads( Apegrunt_options::threads() > 0 ? Apegrunt_options::threads() : tbb_task_scheduler.default_num_threads() );
		//Apegrunt_options::threads() > 0 ? tbb_task_scheduler.initialize( Apegrunt_options::threads() ) : tbb_task_scheduler.initialize() ); // Threading task scheduler
		#endif // #ifndef APEGRUNT_NO_TBB

		// Apegrunt options
		apegrunt_options.set_verbose( Apegrunt_options::verbose() );
		if( !apegrunt_options.CheckOptions(&options_map) )
		{
			exit(EXIT_FAILURE);
		}
		apegrunt_options.set_threads( Apegrunt_options::threads() );

		std::cout << Apegrunt_options::s_get_version_string() << "\n"
				  << Apegrunt_options::s_get_copyright_notice_string() << "\n"
				  << std::endl;

		#ifndef APEGRUNT_NO_TBB // Threading with Threading Building Blocks
		Apegrunt_options::threads() > 0 ? tbb_task_scheduler.initialize( Apegrunt_options::threads() ) : tbb_task_scheduler.initialize(); // Threading task scheduler
		if( Apegrunt_options::verbose() )
		{
			*Apegrunt_options::get_out_stream()
				<< "Apegrunt: TBB interface version " << TBB_INTERFACE_VERSION << "\n"
				<< "Apegrunt: TBB runtime interface version " << tbb::TBB_runtime_interface_version() << "\n"
				<< "Apegrunt: TBB task scheduler is " << ( tbb_task_scheduler.is_active() ? "ACTIVE" : "INACTIVE" );
			if( tbb_task_scheduler.is_active() )
			{
				*Apegrunt_options::get_out_stream()
					<< ": using "
					//<< ( Apegrunt_options::threads() > 0 ? Apegrunt_options::threads() : tbb_task_scheduler.default_num_threads() )
					<< Apegrunt_options::threads()
					<< " threads"
				;
			}
			*Apegrunt_options::get_out_stream() << "\n" << std::endl;
		}
		//aracne::ARACNE_options::set_threads( Apegrunt_options::threads() > 0 ? Apegrunt_options::threads() : tbb_task_scheduler.default_num_threads() );
		#endif // #ifndef APEGRUNT_NO_TBB
	}

	catch( std::exception& e )
	{
		// probably an unknown option
		*Apegrunt_options::get_err_stream() << "Apegrunt error: " << e.what() << "\n\n";
		*Apegrunt_options::get_out_stream() << Apegrunt_options::s_get_usage_string() << std::endl;
		exit(EXIT_FAILURE);
	}
	catch(...)
	{
		*Apegrunt_options::get_err_stream() << "Apegrunt error: Exception of unknown type!\n\n";
		exit(EXIT_FAILURE);
	}

	// setup global timer
	stopwatch::stopwatch globaltimer( Apegrunt_options::verbose() ? Apegrunt_options::get_out_stream() : nullptr ); // for timing statistics
	globaltimer.start();

	apegrunt::Loci_ptr sample_list;

	stopwatch::stopwatch cputimer( Apegrunt_options::verbose() ? Apegrunt_options::get_out_stream() : nullptr ); // for timing statistics

	// get alignments
	auto alignments = apegrunt::get_alignments<default_state_t>( 1 );

	for( auto alignment: alignments )
	{
		cputimer.start();
		{
			auto alignment_file = apegrunt::get_unique_ofstream( alignment->id_string()+"."+apegrunt::size_string(alignment)+".agf" );
			if( Apegrunt_options::verbose() )
			{
				*Apegrunt_options::get_out_stream() << "Apegrunt: write AGF alignment to file " << alignment_file->name() << "\n";
			}
			generate_AG_Alignment( alignment, alignment_file->stream() );
		}
		cputimer.stop();
		if( Apegrunt_options::verbose() ) { cputimer.print_timing_stats(); *Apegrunt_options::get_out_stream() << "\n"; }
	}


	// Check if we have the compulsory input alignment
	if( alignments.empty() ) { exit(EXIT_FAILURE); }

// /*
	stopwatch::stopwatch steptimer( Apegrunt_options::verbose() ? Apegrunt_options::get_out_stream() : nullptr ); // for timing statistics
	if( Apegrunt_options::verbose() )
	{
		*Apegrunt_options::get_out_stream() << "Apegrunt: pre-process " << alignments.size() << " alignment" << (alignments.size() > 1 ? "s" : "") << ":\n" << std::endl;
	}
	steptimer.start();
	for( auto& alignment: alignments )
	{
		if( Apegrunt_options::verbose() )
		{
			*Apegrunt_options::get_out_stream() << "Apegrunt: alignment \"" << alignment->id_string() << "\":\n";
		}
		if( alignments.size() < 2 )
		{
			if( apegrunt::Apegrunt_options::has_includelist_filename() )
			{
				if( Apegrunt_options::verbose() )
				{
					*Apegrunt_options::get_out_stream() << "Apegrunt: get include list from file \"" << apegrunt::Apegrunt_options::get_includelist_filename() << "\"\n";
				}
				cputimer.start();
				auto include_list = apegrunt::parse_Loci_list( apegrunt::Apegrunt_options::get_includelist_filename(), apegrunt::Apegrunt_options::get_input_indexing_base() );
				cputimer.print_timing_stats();

				if( Apegrunt_options::verbose() )
				{
					*Apegrunt_options::get_out_stream() << "Apegrunt: include list has " << include_list->size() << " loci.\n";
					*Apegrunt_options::get_out_stream() << "Apegrunt: trim alignment based on include list \"" << include_list->id_string() << "\"\n";
				}
				cputimer.start();
				alignment = apegrunt::Alignment_factory< alignment_default_storage_t >().include( alignment, include_list );
				cputimer.stop();
				if( Apegrunt_options::verbose() ) { cputimer.print_timing_stats(); *Apegrunt_options::get_out_stream() << "\n"; }
			}

			if( apegrunt::Apegrunt_options::has_excludelist_filename() )
			{
				if( Apegrunt_options::verbose() )
				{
					*Apegrunt_options::get_out_stream() << "Apegrunt: get exclude list from file \"" << apegrunt::Apegrunt_options::get_excludelist_filename() << "\"\n";
				}
				cputimer.start();
				auto exclude_list = apegrunt::parse_Loci_list( apegrunt::Apegrunt_options::get_excludelist_filename(), apegrunt::Apegrunt_options::get_input_indexing_base() );
				cputimer.print_timing_stats();

				if( Apegrunt_options::verbose() )
				{
					*Apegrunt_options::get_out_stream() << "Apegrunt: exclude list has " << exclude_list->size() << " loci.\n";
					*Apegrunt_options::get_out_stream() << "Apegrunt: trim alignment based on exclude list \"" << exclude_list->id_string() << "\"\n";
				}
				cputimer.start();
				alignment = apegrunt::Alignment_factory< alignment_default_storage_t >().exclude( alignment, exclude_list );
				cputimer.stop();
				if( Apegrunt_options::verbose() ) { cputimer.print_timing_stats(); *Apegrunt_options::get_out_stream() << "\n"; }
			}
		}

		if( apegrunt::Apegrunt_options::filter_alignment() )
		{
	    	// apply filters as defined on the command line
			cputimer.start();
			if( Apegrunt_options::verbose() )
			{
				*Apegrunt_options::get_out_stream() << "Apegrunt: apply filter rules";
				Apegrunt_options::get_out_stream()->flush();
			}
			auto alignment_filter = apegrunt::Alignment_filter( apegrunt::Alignment_filter::ParameterPolicy::AQUIRE_GLOBAL );
			//auto alignment_filter = apegrunt::Alignment_filter( apegrunt::Alignment_filter::ParameterPolicy::FILTER_SNPS );
			alignment = alignment_filter.operator()<alignment_default_storage_t>( alignment );
			cputimer.stop();
			if( Apegrunt_options::verbose() )
			{
				*Apegrunt_options::get_out_stream() << "\n";
				alignment->statistics( Apegrunt_options::get_out_stream() );
				cputimer.print_timing_stats(); *Apegrunt_options::get_out_stream() << "\n";
			}

		}

		if( sample_list )
		{
			if( Apegrunt_options::verbose() )
			{
				*Apegrunt_options::get_out_stream() << "Apegrunt: trim alignment samples based on sample list \"" << sample_list->id_string() << "\"\n";
			}
			cputimer.start();
			alignment = apegrunt::Alignment_factory< alignment_default_storage_t >().copy_selected( alignment, sample_list, sample_list->id_string() );
			cputimer.stop();
			if( Apegrunt_options::verbose() ) { cputimer.print_timing_stats(); }
 		}

		// get state frequency profile and output to file
		apegrunt::output_state_frequencies( alignment );

		// output the sample-sample Hamming distance matrix
		apegrunt::output_sample_distance_matrix( alignment );

		// assign sample weights (parse from file or determine automatically)
		apegrunt::cache_sample_weights( alignment );

		// output sample weights
		apegrunt::output_sample_weights( alignment );
	}

	steptimer.stop();
	if( Apegrunt_options::verbose() )
	{
		*Apegrunt_options::get_out_stream() << "Apegrunt: alignment pre-processing completed\n";
		steptimer.print_timing_stats();	*Apegrunt_options::get_out_stream() << "\n";
	}

	if( Apegrunt_options::verbose() )
	{
		*Apegrunt_options::get_out_stream() << "Apegrunt: analysis completed\n";
		globaltimer.stop(); globaltimer.print_timing_stats();
	}

	exit(EXIT_SUCCESS);
}
