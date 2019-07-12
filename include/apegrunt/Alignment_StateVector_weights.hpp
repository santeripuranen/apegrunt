/** @file Alignment_StateVector_weights.hpp

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

#ifndef APEGRUNT_ALIGNMENT_STATEVECTOR_WEIGHTS_HPP
#define APEGRUNT_ALIGNMENT_STATEVECTOR_WEIGHTS_HPP

#include "misc/Stopwatch.hpp"
#include "StateVector_forward.h"

namespace apegrunt {

template< typename StateT, typename RealT=double >
std::vector<RealT> calculate_weights( const apegrunt::Alignment_ptr<StateT> alignment )
{
	stopwatch::stopwatch cputimer(Apegrunt_options::get_out_stream()); // for timing statistics

	using real_t = RealT;

	const std::size_t n_seqs = alignment->size();
	auto ali = alignment->subscript_proxy();

	std::vector<std::size_t> int_weights; int_weights.reserve( n_seqs );
	for( const auto& seq_ptr: *alignment ) { int_weights.push_back( seq_ptr->multiplicity() ); }

	if( Apegrunt_options::verbose() ) { *Apegrunt_options::get_out_stream() << "apegrunt: calculate weights"; Apegrunt_options::get_out_stream()->flush(); }
	cputimer.start();
	std::size_t nonidentical_pairs = 0;

	//std::size_t sum = 0;
	for( std::size_t i = 0; i < n_seqs; ++i )
	{
		const auto seq_i = ali[i]; // (*alignment)[i];
		const auto seq_i_size = seq_i->size();
		const auto seq_i_multiplicity = seq_i->multiplicity();
		for( std::size_t j = 0; j < i; ++j )
		{
			if( seq_i->is_similar_to( *(ali[j]), Apegrunt_options::sample_reweighting_threshold()*double(std::min(seq_i_size,ali[j]->size())) ) )
			{
				int_weights[i] += ali[j]->multiplicity();
				int_weights[j] += seq_i_multiplicity;
			}
			else
			{
				++nonidentical_pairs;
			}
		}
	}
	cputimer.stop();
	if( Apegrunt_options::verbose() ) { cputimer.print_timing_stats(); }

	if( Apegrunt_options::verbose() )
	{
		const std::size_t threshold = Apegrunt_options::sample_reweighting_threshold() * double(alignment->n_loci());
		const auto total_pairs = (n_seqs*(n_seqs-1))/2;
		*Apegrunt_options::get_out_stream()
			<< "apegrunt: reweighting threshold=" << threshold << "\n"
			<< "apegrunt: # of non-identical pairs=" << nonidentical_pairs << " out of " << total_pairs << " (ratio=" << real_t(nonidentical_pairs)/real_t(total_pairs) << ")\n";
	}

	std::vector<real_t> weights; weights.reserve(n_seqs);

	using boost::get;
	for( auto&& int_weights_and_seq : zip_range(int_weights,*alignment) ) { weights.push_back( get<1>(int_weights_and_seq)->multiplicity() / real_t( get<0>(int_weights_and_seq) ) ); } // 1/weight * multiplicity
	//for( auto&& int_weights_and_seq : zip_range(int_weights,*alignment) ) { weights.push_back( ( ( get<1>(int_weights_and_seq)->multiplicity() / real_t( get<0>(int_weights_and_seq) ) ) / real_t(sum) ) * real_t(n_seqs) ); } // 1/weight * multiplicity

	return weights;
}

template< typename StateT, typename RealT=double >
std::vector<RealT> calculate_weights_w_cached_distance_matrix( apegrunt::Alignment_ptr<StateT> alignment )
{
	stopwatch::stopwatch cputimer(Apegrunt_options::get_out_stream()); // for timing statistics

	using real_t = RealT;

	auto ali = alignment->subscript_proxy();

	const real_t threshold = Apegrunt_options::sample_reweighting_threshold();
	const std::size_t n_seqs = alignment->size();

	std::vector<std::size_t> int_weights; int_weights.reserve( n_seqs );
	for( const auto& seq_ptr: *alignment ) { int_weights.push_back( seq_ptr->multiplicity() ); }

	real_t min_identity = 1.0;
	real_t max_identity = 0.0;

	if( Apegrunt_options::verbose() ) { *Apegrunt_options::get_out_stream() << "apegrunt: get distance matrix"; Apegrunt_options::get_out_stream()->flush(); }
	cputimer.start();
	const auto dmat = *(alignment->distance_matrix());
	cputimer.stop();
	if( Apegrunt_options::verbose() ) { cputimer.print_timing_stats(); }

	if( Apegrunt_options::verbose() ) { *Apegrunt_options::get_out_stream() << "apegrunt: calculate weights"; Apegrunt_options::get_out_stream()->flush(); }
	cputimer.start();
	std::size_t nonidentical_pairs = 0;
	//std::size_t sum = 0;
	for( std::size_t i = 0; i < n_seqs; ++i )
	{
		const auto seq_i = ali[i]; // (*alignment)[i];
		const auto seq_i_size = seq_i->size();
		const auto seq_i_multiplicity = seq_i->multiplicity();
		for( std::size_t j = 0; j < i; ++j )
		{

			const std::size_t n_elem = std::min( seq_i_size, ali[j]->size() );
			const auto identity = real_t( n_elem - dmat[ i*(i-1)/2+j ] ) / real_t(n_elem);
			min_identity = std::min( min_identity, identity );
			max_identity = std::max( max_identity, identity );

			if( identity < threshold )
			{
				++nonidentical_pairs;
			}
			else
			{
				//int_weights[i] += (*alignment)[j]->multiplicity();
				int_weights[i] += ali[j]->multiplicity();
				int_weights[j] += seq_i_multiplicity;
				//sum += ( ali[j]->multiplicity() + seq_i_multiplicity );
			}
		}
	}
	cputimer.stop();
	if( Apegrunt_options::verbose() ) { cputimer.print_timing_stats(); }

	if( Apegrunt_options::verbose() )
	{
		const auto total_pairs = (n_seqs*(n_seqs-1))/2;
		*Apegrunt_options::get_out_stream()
			<< "apegrunt: sequence identity min=" << min_identity << " max=" << max_identity << " reweighting threshold=" << threshold << "\n"
			<< "apegrunt: # of non-identical pairs=" << nonidentical_pairs << " out of " << total_pairs << " (ratio=" << real_t(nonidentical_pairs)/real_t(total_pairs) << ")\n";
	}

	std::vector<real_t> weights; weights.reserve(n_seqs);

	using boost::get;
	for( auto&& int_weights_and_seq : zip_range(int_weights,*alignment) ) { weights.push_back( get<1>(int_weights_and_seq)->multiplicity() / real_t( get<0>(int_weights_and_seq) ) ); } // 1/weight * multiplicity
	//for( auto&& int_weights_and_seq : zip_range(int_weights,*alignment) ) { weights.push_back( ( ( get<1>(int_weights_and_seq)->multiplicity() / real_t( get<0>(int_weights_and_seq) ) ) / real_t(sum) ) * real_t(n_seqs) ); } // 1/weight * multiplicity

	return weights;
}

template< typename StateT, typename RealT=double >
std::vector<RealT> get_weights( apegrunt::Alignment_ptr<StateT> alignment )
{
	using real_t = RealT;
	std::vector<real_t> weights; weights.reserve(alignment->size());

	for( const auto seq: alignment )
	{
		weights.push_back( seq->weight() );
	}

	return weights;
}

} // namespace apegrunt

#endif // APEGRUNT_ALIGNMENT_STATEVECTOR_WEIGHTS_HPP
