/** @file Alignment_utility.hpp

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
#ifndef APEGRUNT_ALIGNMENT_UTILITY_HPP
#define APEGRUNT_ALIGNMENT_UTILITY_HPP

#include <vector>
#include <array>
#include <algorithm>
#include <sstream>
#include <iomanip> // std::setfill, std::setw
#include <random> // C++11
#include <type_traits> // std::true_type

#include "Apegrunt_utility.hpp"
#include "Apegrunt_IO_misc.hpp"
#include "Threshold_rule.hpp"

#include "Alignment_forward.h"
#include "StateVector_forward.h"
#include "StateVector_state_types.hpp"
#include "StateVector_utility.hpp"
#include "Alignment_StateVector_weights.hpp"

#include "Loci_parsers.hpp"
#include "Loci_generators.hpp"
#include "State_block.hpp"

#include "ValueVector_parser.hpp"

#include "aligned_allocator.hpp"

#include "accumulators/distribution.hpp"
#include "accumulators/distribution_std.hpp"
#include "accumulators/distribution_bincount.hpp"
#include "accumulators/distribution_ordered.hpp"
#include "accumulators/distribution_cumulative.hpp"
//#include "accumulators/distribution_generator_svg.hpp"
#include "accumulators/distribution_generator_csv.hpp"
#include "accumulators/accumulators_utility.hpp"

#include "graph/Graph.hpp"

#include "misc/Stopwatch.hpp"

namespace apegrunt {

template< typename StateT >
std::string size_string( Alignment_ptr<StateT> alignment ) { std::ostringstream label; label << "L" << alignment->n_loci() << "n" << alignment->size(); return label.str(); }

std::size_t norm2( const std::vector< std::size_t >& v )
{
	using std::pow;

	return std::accumulate( cbegin(v), cend(v), std::size_t(0), [=]( std::size_t sum, std::size_t x ) { return sum += x*x; } );
}

template< typename RealT >
RealT norm( const std::vector< std::size_t >& v )
{
	using real_t = RealT;
	using std::sqrt;
	auto n2 = norm2(v);
	return ( n2 != 0.0 ? sqrt( real_t(n2) ) : 0 );
}

template< typename StateT >
std::vector< Alignment_ptr<StateT> > get_alignments( std::size_t max_alignments=std::numeric_limits<std::size_t>::max() )
{
	using state_t = StateT;
	using alignment_default_storage_t = apegrunt::Alignment_impl_block_compressed_storage< apegrunt::StateVector_impl_block_compressed_alignment_storage<state_t> >;

	stopwatch::stopwatch cputimer( Apegrunt_options::verbose() ? Apegrunt_options::get_out_stream() : nullptr ); // for timing statistics

	std::vector< Alignment_ptr<state_t> > alignments;

	if( apegrunt::Apegrunt_options::has_alignment_filenames() )
	{
		if( apegrunt::Apegrunt_options::get_alignment_filenames().size() > max_alignments )
		{
			*Apegrunt_options::get_err_stream() << "apegrunt error: will not accept more than " << max_alignments << " alignment" << (max_alignments == 1 ? " " : "s") << "(got " << apegrunt::Apegrunt_options::get_alignment_filenames().size() << " filenames).\n\n";
			return alignments;
		}

		if( Apegrunt_options::verbose() )
		{
			*Apegrunt_options::get_out_stream() << "apegrunt: get " << apegrunt::Apegrunt_options::get_alignment_filenames().size() << " alignment" << (apegrunt::Apegrunt_options::get_alignment_filenames().size() > 1 ? "s" : "") << ":\n\n";
		}

		for( auto& filename: apegrunt::Apegrunt_options::get_alignment_filenames() )
		{
			if( Apegrunt_options::verbose() )
			{
				*Apegrunt_options::get_out_stream() << "apegrunt: parse alignment from file \"" << filename << "\"\n";
				Apegrunt_options::get_out_stream()->flush();
			}

			cputimer.start();
			auto alignment = apegrunt::parse_Alignment< alignment_default_storage_t >( filename );
			cputimer.stop();

			if( !alignment )
			{
				*Apegrunt_options::get_err_stream() << "apegrunt error: Could not get alignment from input file \"" << filename << "\"\n\n";
			}
			else
			{
				cputimer.print_timing_stats(); *Apegrunt_options::get_out_stream() << "\n";
				alignments.push_back( alignment ); // store the new alignment
				cputimer.start();
				if( Apegrunt_options::verbose() )
				{
					//*Apegrunt_options::get_out_stream() << "apegrunt: statistics for alignment \"" << alignment->id_string() << "\":\n";
					alignment->statistics( Apegrunt_options::get_out_stream() );
					//*Apegrunt_options::get_out_stream() << "done\n"; Apegrunt_options::get_out_stream()->flush();
				}
				cputimer.stop();
			}
			cputimer.print_timing_stats(); *Apegrunt_options::get_out_stream() << "\n";
		}
	}
	else
	{
		*Apegrunt_options::get_err_stream() << "apegrunt: no input files specified!\n\n";
		exit(EXIT_FAILURE);
	}
	if( apegrunt::Apegrunt_options::has_mappinglist_filename() )
	{
		if( Apegrunt_options::verbose() )
		{
			*Apegrunt_options::get_out_stream() << "apegrunt: get mapping list from file \"" << apegrunt::Apegrunt_options::get_mappinglist_filename() << "\"\n";
		}
		cputimer.start();
		auto mapping = apegrunt::parse_Loci_list( apegrunt::Apegrunt_options::get_mappinglist_filename(), apegrunt::Apegrunt_options::get_input_indexing_base() );
		cputimer.stop();

		auto alignment = alignments.front();

		if( alignment->n_loci() != mapping->size() )
		{
			*Apegrunt_options::get_err_stream() << "apegrunt: mapping list and alignment sizes are mismatched (" << mapping->size() << " and " << alignment->n_loci() << ")\n\n";
			cputimer.print_timing_stats();
			exit(EXIT_FAILURE);
		}
		else
		{
			if( Apegrunt_options::verbose() )
			{
				*Apegrunt_options::get_out_stream() << "apegrunt: associate mapping list with \"" << alignment->id_string() << "\"\n";
			}

			alignment->set_loci_translation( mapping );

			if( Apegrunt_options::verbose() )
			{
				cputimer.print_timing_stats(); *Apegrunt_options::get_out_stream() << "\n";
			}
		}
	}
/*
	// output alignments? I guess this is mostly useful for testing purposes
	if( apegrunt::Apegrunt_options::output_alignment() )
	{
		for( auto alignment: alignments )
		{
			apegrunt::output_alignment( alignment );
		}
	}

	// output column-wise frequencies?
	if( Apegrunt_options::output_state_frequencies() )
	{
		for( auto alignment: alignments )
		{
			// get state frequency profile and output to file
			apegrunt::output_state_frequencies( alignment );
		}
	}
*/
	return alignments;
}

// deal with sample weights
template< typename StateT, typename RealT=double >
void cache_sample_weights( Alignment_ptr<StateT> alignment )
{
	using real_t = RealT;
	using std::cbegin; using std::cend;

	stopwatch::stopwatch cputimer( Apegrunt_options::verbose() ? Apegrunt_options::get_out_stream() : nullptr ); // for timing statistics

	if( Apegrunt_options::reweight_samples() )
	{
		cputimer.start();

		if( Apegrunt_options::has_sample_weights_filename() )
		{
			if( Apegrunt_options::verbose() )
			{
				*Apegrunt_options::get_out_stream() << "apegrunt: get sample weights from file \"" << Apegrunt_options::get_sample_weights_filename() << "\"\n";
			}

			auto weights = apegrunt::parse_ValueVector<real_t>( Apegrunt_options::get_sample_weights_filename() );
			const auto n_size = alignment->size();
			using apegrunt::cbegin; using apegrunt::cend;
			const auto n_eff = std::accumulate( cbegin(weights), cend(weights), real_t(0) );

			for( auto seq_and_weight: apegrunt::zip_range( alignment, weights ) )
			{
				using boost::get;
				if( Apegrunt_options::rescale_sample_weights() )
				{
					// weights normalized and scaled to give neff == n_size
					get<0>(seq_and_weight)->set_weight( real_t(get<1>(seq_and_weight)) * real_t(n_size)/real_t(n_eff) );
				}
				else
				{
					// unnormalized and unscaled weights, neff != n_size
					get<0>(seq_and_weight)->set_weight( get<1>(seq_and_weight) );
				}
			}
		}
		else
		{
			if( Apegrunt_options::verbose() )
			{
				*Apegrunt_options::get_out_stream() << "apegrunt: calculate sample weights" << std::endl;
			}

			auto weights = calculate_weights( alignment );
			const auto n_size = alignment->size();
			const auto n_eff = std::accumulate( cbegin(weights), cend(weights), real_t(0) );

			for( auto seq_and_weight: zip_range( alignment, weights ) )
			{
				using boost::get;
				if( Apegrunt_options::rescale_sample_weights() )
				{
					// weights scaled to give neff == n_size
					get<0>(seq_and_weight)->set_weight( real_t(get<1>(seq_and_weight)) * real_t(n_size)/real_t(n_eff) );
				}
				else
				{
					// unscaled weights neff != n_size
					get<0>(seq_and_weight)->set_weight( get<1>(seq_and_weight) );
				}
			}
		}

		cputimer.stop();
		if( Apegrunt_options::verbose() )
		{
			*Apegrunt_options::get_out_stream() << "apegrunt: effective sample size = " << alignment->effective_size() << "\n";
			cputimer.print_timing_stats(); *Apegrunt_options::get_out_stream() << "\n";
		}
	}
}

template< typename StateT, typename RealT=double >
void output_sample_weights( Alignment_ptr<StateT> alignment )
{
	// output weights
	if( Apegrunt_options::output_sample_weights() )
	{
		// init timer
		stopwatch::stopwatch cputimer( Apegrunt_options::verbose() ? Apegrunt_options::get_out_stream() : nullptr ); // for timing statistics

		// output weights
		auto weights_file = get_unique_ofstream( alignment->id_string()+"."+apegrunt::size_string(alignment)+".weights" );
		auto& weights_stream = *weights_file->stream();
		weights_stream << std::scientific;
		weights_stream.precision(8);

		//auto weights = apegrunt::get_weights( alignment );

		if( Apegrunt_options::verbose() )
		{
			*Apegrunt_options::get_out_stream() << "apegrunt: write sample weights to file \"" << alignment->effective_size() << "\"\n";
			Apegrunt_options::get_out_stream()->flush();
		}
		cputimer.start();
		for( auto sequence: alignment ) { weights_stream << sequence->weight() << "\n"; }
		weights_stream.close();
		cputimer.stop();
		if( Apegrunt_options::verbose() )
		{
			cputimer.print_timing_stats();
			*Apegrunt_options::get_out_stream() << "\n";
		}
	}
}

template< typename StateT >
void output_alignment( Alignment_ptr<StateT> alignment, bool output_indices=true )
{
	// init timer
	stopwatch::stopwatch cputimer( Apegrunt_options::verbose() ? Apegrunt_options::get_out_stream() : nullptr ); // for timing statistics

	cputimer.start();

	// output alignment
	auto alignment_file = apegrunt::get_unique_ofstream( alignment->id_string()+"."+apegrunt::size_string(alignment)+"."+"fasta" );
	if( Apegrunt_options::verbose() )
	{
		*Apegrunt_options::get_out_stream() << "apegrunt: write alignment to file " << alignment_file->name() << "\n";
	}
	apegrunt::generate_Alignment( alignment, alignment_file->stream() );

	// output loci list
	if( output_indices )
	{
		auto locilist_file = apegrunt::get_unique_ofstream( alignment->id_string()+"."+apegrunt::size_string(alignment)+".loci"  );
		if( Apegrunt_options::verbose() )
		{
			*Apegrunt_options::get_out_stream() << "apegrunt: write original loci indices to file " << locilist_file->name() << "\n";
		}
		apegrunt::generate_Loci_list( alignment->get_loci_translation(), locilist_file->stream() );
	}

	cputimer.stop();
	if( Apegrunt_options::verbose() ) { cputimer.print_timing_stats(); *Apegrunt_options::get_out_stream() << "\n"; }
}


template< typename StateT >
std::vector< std::array< std::size_t, number_of_states<StateT>::N > > allele_occurrence( const Alignment_ptr<StateT> alignment )
{
	using state_t = StateT;
	enum { N=number_of_states<state_t>::N };
	using boost::get;

	std::vector< std::array<std::size_t,N> > occ_vec( alignment->n_loci(), {0} );

	for( const auto& sequence: alignment )
	{
		const std::size_t multiplicity = sequence->multiplicity();
		for( auto&& occ_and_state: zip_range(occ_vec, sequence) )
		{
			get<0>(occ_and_state)[ std::size_t( get<1>(occ_and_state) ) ] += multiplicity;
		}
		//std::cout << "seq=" << sequence->id_string() << std::endl;
	}
	return occ_vec;
}

using WeightedEdges = std::vector< apegrunt::Edge<true> >;

template< typename ContainerT, typename RealT >
WeightedEdges inplace_block_list_intersection_weights( const std::vector< ContainerT >& a, std::vector< ContainerT >& b, const std::vector<RealT>& sample_weights )
{
	using std::cbegin; using std::cend; using std::begin;

	using real_t = RealT;

	WeightedEdges intersections; intersections.reserve(a.size()*b.size());

	for( std::size_t i=0; i < a.size(); ++i )
	{
		auto ai( a[i] );
		for( std::size_t j=0; j < b.size() && !ai.empty(); ++j )
		{
			if( b[j].empty() ) { continue; }
			else
			{
				if( has_overlap(ai,b[j]) )
				{
					const auto result = apegrunt::inplace_set_differences_and_for_each_in_intersection( ai, b[j], apegrunt::gatherer<real_t>(sample_weights.data()) );
					if( result.sum != 0 )
					{
						intersections.emplace_back(i,j,result.sum);
					}
				}
			}
		}
	}

	return intersections;
}

template< typename ContainerT, typename RealT >
WeightedEdges block_list_intersection_weights( const std::vector< ContainerT >& a, const std::vector< ContainerT >& b, const std::vector<RealT>& sample_weights )
{
	auto _b(b);

	return inplace_block_list_intersection_weights(a,_b,sample_weights);
}

template< typename ContainerT, typename RealT >
WeightedEdges inplace_block_list_intersection_weights_dynamic( const std::vector< ContainerT >& a, std::vector< ContainerT >& b, const std::vector<RealT>& sample_weights )
{
	using std::cbegin; using std::cend; using std::begin;

	using real_t = RealT;

	WeightedEdges intersections; intersections.reserve(a.size()*b.size());

	for( std::size_t i=0; i < a.size(); ++i )
	{
		//auto& ai = a[i];
		auto ai( a[i] ); // make a copy
		for( std::size_t j=0; j < b.size() && !ai.empty(); ++j )
		{
			if( b[j].empty() ) { continue; }
			else
			{
				if( has_overlap(ai,b[j]) )
				{
					const auto result = apegrunt::inplace_set_differences_and_for_each_in_intersection( ai, b[j], apegrunt::gatherer<real_t>(sample_weights.data()) );
					if( result.sum != 0 )
					{
						intersections.emplace_back(i,j,result.sum);
					}
				}
			}
		}
	}

	return intersections;
}

template< typename ContainerT, typename RealT >
WeightedEdges block_list_intersection_weights_dynamic( const std::vector< ContainerT >& a, const std::vector< ContainerT >& b, const std::vector<RealT>& sample_weights )
{
	//auto _a(a);
	auto _b(b); // make a copy

	return inplace_block_list_intersection_weights_dynamic(a,_b,sample_weights);
}

template< typename StateT >
bool output_frequencies( const Alignment_ptr<StateT> alignment, std::ostream *os=nullptr, bool weighted=false )
{
	if( os && os->good() )
	{
		//os->precision(3); //os->width(8);
		const auto freq_vec = weighted ? weighted_columnwise_frequencies( alignment ) : columnwise_frequencies( alignment );
		for( const auto& freq: freq_vec ) // loop over all sequence positions
		{
			for( const auto& f: freq )
			{
				*os << f << " ";
			}
			*os << "\n";
		}
		return true;
	}
	return false;
}

template< typename StateT >
bool output_frequency_distribution( const Alignment_ptr<StateT> alignment, std::ostream *os=nullptr, bool weighted=false )
{
	namespace acc = boost::accumulators;
	using real_t = double;

	if( os && os->good() )
	{
		acc::accumulator_set<real_t, acc::stats<acc::tag::std(acc::from_distribution),acc::tag::distribution_bincount> >
		frequency_distribution( acc::tag::distribution::binwidth=0.01 );
		//os->precision(3); //os->width(8);
		const auto freq_vec = weighted ? weighted_columnwise_frequencies( alignment ) : columnwise_frequencies( alignment );
		frequency_distribution << freq_vec;
		*os << apegrunt::accumulators::csv(acc::distribution(frequency_distribution));
		return true;
	}
	return false;
}

template< typename StateT >
void output_state_frequencies( const Alignment_ptr<StateT> alignment )
{
	stopwatch::stopwatch cputimer( Apegrunt_options::verbose() ? Apegrunt_options::get_out_stream() : nullptr ); // for timing statistics

	if( Apegrunt_options::output_state_frequencies() )
	{
		{ // output un-weighted frequencies
			cputimer.start();
			auto frequencies_file = apegrunt::get_unique_ofstream( alignment->id_string()+"."+apegrunt::size_string(alignment)+".frequencies.txt" );
			if( Apegrunt_options::verbose() )
			{
				*Apegrunt_options::get_out_stream() << "apegrunt: write columnwise state frequencies to file \"" << frequencies_file->name() << "\"\n";
				Apegrunt_options::get_out_stream()->flush();
			}
			if( apegrunt::output_frequencies( alignment, frequencies_file->stream() ) )
			{
				cputimer.stop();
				if( Apegrunt_options::verbose() )
				{
					cputimer.print_timing_stats();
					*Apegrunt_options::get_out_stream() << "\n";
				}
			}
			else
			{
				cputimer.stop();
				if( Apegrunt_options::verbose() )
				{
					*Apegrunt_options::get_err_stream() << "apegrunt ERROR: unable to write file \"" << frequencies_file->name() << "\"\n" << std::endl;
				}
			}
		}
		{ // output un-weighted frequency distribution
			cputimer.start();
			auto distribution_file = apegrunt::get_unique_ofstream( alignment->id_string()+"."+apegrunt::size_string(alignment)+".frequency_distribution.txt" );
			if( Apegrunt_options::verbose() )
			{
				*Apegrunt_options::get_out_stream() << "apegrunt: write columnwise state frequency distribution to file \"" << distribution_file->name() << "\"\n";
				Apegrunt_options::get_out_stream()->flush();
			}
			if( apegrunt::output_frequency_distribution( alignment, distribution_file->stream() ) )
			{
				cputimer.stop();
				if( Apegrunt_options::verbose() )
				{
					cputimer.print_timing_stats();
					*Apegrunt_options::get_out_stream() << "\n";
				}
			}
			else
			{
				cputimer.stop();
				if( Apegrunt_options::verbose() )
				{
					*Apegrunt_options::get_err_stream() << "apegrunt ERROR: unable to write file \"" << distribution_file->name() << "\"\n" << std::endl;
				}
			}
		}

		if( apegrunt::Apegrunt_options::reweight_samples() )
		{
			{ // output weighted frequencies
				cputimer.start();
				auto frequencies_file = apegrunt::get_unique_ofstream( alignment->id_string()+"."+apegrunt::size_string(alignment)+".weighted_frequencies.txt" );
				if( Apegrunt_options::verbose() )
				{
					*Apegrunt_options::get_out_stream() << "apegrunt: write weighted columnwise state frequencies to file \"" << frequencies_file->name() << "\"\n";
				}
				if( apegrunt::output_frequencies( alignment, frequencies_file->stream(), true ) ) // true == output weighted frequencies
				{
					cputimer.stop();
					if( Apegrunt_options::verbose() )
					{
						cputimer.print_timing_stats();
						*Apegrunt_options::get_out_stream() << "\n";
					}
				}
				else
				{
					cputimer.stop();
					if( Apegrunt_options::verbose() )
					{
						*Apegrunt_options::get_err_stream() << "apegrunt ERROR: unable to write file \"" << frequencies_file->name() << "\"\n" << std::endl;
					}
				}
			}
			{ // output un-weighted frequency distribution
				cputimer.start();
				auto distribution_file = apegrunt::get_unique_ofstream( alignment->id_string()+"."+apegrunt::size_string(alignment)+".weighted_frequency_distribution.txt" );
				if( Apegrunt_options::verbose() )
				{
					*Apegrunt_options::get_out_stream() << "apegrunt: write columnwise state frequency distribution to file \"" << distribution_file->name() << "\"\n";
					Apegrunt_options::get_out_stream()->flush();
				}
				if( apegrunt::output_frequency_distribution( alignment, distribution_file->stream(), true ) ) // true == output weighted frequencies
				{
					cputimer.stop();
					if( Apegrunt_options::verbose() )
					{
						cputimer.print_timing_stats();
						*Apegrunt_options::get_out_stream() << "\n";
					}
				}
				else
				{
					cputimer.stop();
					if( Apegrunt_options::verbose() )
					{
						*Apegrunt_options::get_err_stream() << "apegrunt ERROR: unable to write file \"" << distribution_file->name() << "\"\n" << std::endl;
					}
				}
			}
		}
	}
}

template< typename StateT >
bool output_sample_distance_matrix_impl( const Alignment_ptr<StateT> alignment, std::ostream *os=nullptr )
{
	if( os && os->good() )
	{
		auto distance_matrix_ptr = alignment->distance_matrix();
		const auto& distance_matrix = *(distance_matrix_ptr.get());
		const auto n_samples = alignment->size();

		*os << "n=" << n_samples;
		for( std::size_t i = 0; i < n_samples; ++i )
		{
			//const auto sample_i = distance_matrix[i];
			for( std::size_t j = 0; j < i; ++j )
			{
				*os << distance_matrix[i*(i-1)/2+j] << " ";
			}
			*os << "\n";
		}

		return true;
	}
	return false;
}

template< typename StateT >
void output_sample_distance_matrix( const Alignment_ptr<StateT> alignment )
{
	stopwatch::stopwatch cputimer( Apegrunt_options::verbose() ? Apegrunt_options::get_out_stream() : nullptr ); // for timing statistics

	if( Apegrunt_options::output_sample_distance_matrix() )
	{
		cputimer.start();
		auto matrix_file = apegrunt::get_unique_ofstream( alignment->id_string()+"."+apegrunt::size_string(alignment)+".triangular_sample_distance_matrix" );

		if( Apegrunt_options::verbose() )
		{
			*Apegrunt_options::get_out_stream() << "apegrunt: write sample-sample Hamming distance matrix to file \"" << matrix_file->name() << "\"\n";
		}
		if( apegrunt::output_sample_distance_matrix_impl( alignment, matrix_file->stream() ) )
		{
			cputimer.stop();
			if( Apegrunt_options::verbose() )
			{
				cputimer.print_timing_stats();
				*Apegrunt_options::get_out_stream() << "\n";
			}
		}
		else
		{
			cputimer.stop();
			if( Apegrunt_options::verbose() )
			{
				*Apegrunt_options::get_err_stream() << "apegrunt ERROR: unable to write file \"" << matrix_file->name() << "\"\n" << std::endl;
			}
		}
	}
}

} // namespace apegrunt

#endif // APEGRUNT_ALIGNMENT_UTILITY_HPP
