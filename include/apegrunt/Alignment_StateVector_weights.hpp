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

#include <memory> // for std::shared_ptr and std::make_shared
#include <atomic>
#include <execution>

#ifdef _OPENMP
#include <omp.h>
#endif

#ifndef NO_INTRINSICS
#include "misc/SIMD_intrinsics.h"
#endif // !NO_INTRINSICS

#include "Alignment_forward.h"
#include "Alignment_interface.hpp"
#include "Apegrunt_options.h"
#include "Apegrunt_utility.hpp"
#include "StateVector_state_types.hpp"

#include "misc/Stopwatch.hpp"
#include "StateVector_forward.h"
#include "StateVector_interface.hpp"

namespace apegrunt {

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
std::vector<RealT> calculate_weights_wo_cached_distance_matrix( apegrunt::Alignment_ptr<StateT> alignment )
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

	if( Apegrunt_options::verbose() ) { *Apegrunt_options::get_out_stream() << "apegrunt: calculate weights (w/o cached distance matrix)"; Apegrunt_options::get_out_stream()->flush(); }
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
			const auto identity = real_t( *seq_i && *(ali[j]) ) / real_t(n_elem); // bitwise-AND counts the number of identical positions
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
			<< "apegrunt: # of non-similar pairs=" << nonidentical_pairs << " out of " << total_pairs << " (ratio=" << real_t(nonidentical_pairs)/real_t(total_pairs) << ")\n";
	}

	std::vector<real_t> weights; weights.reserve(n_seqs);

	using boost::get;
	for( auto&& int_weights_and_seq : zip_range(int_weights,*alignment) ) { weights.push_back( get<1>(int_weights_and_seq)->multiplicity() / real_t( get<0>(int_weights_and_seq) ) ); } // 1/weight * multiplicity
	//for( auto&& int_weights_and_seq : zip_range(int_weights,*alignment) ) { weights.push_back( ( ( get<1>(int_weights_and_seq)->multiplicity() / real_t( get<0>(int_weights_and_seq) ) ) / real_t(sum) ) * real_t(n_seqs) ); } // 1/weight * multiplicity

	return weights;
}

template< typename StateT, typename RealT=double >
std::vector<RealT> calculate_weights( const apegrunt::Alignment_ptr<StateT> alignment )
{
	stopwatch::stopwatch cputimer(Apegrunt_options::get_out_stream()); // for timing statistics

	using real_t = RealT;

	auto ali = alignment->subscript_proxy();

	const std::size_t n_seqs = alignment->size();
	const std::size_t L = alignment->n_loci();
	//const std::size_t threshold = Apegrunt_options::sample_reweighting_threshold() * double(L);

	using boost::get;
	std::vector< std::atomic_size_t > int_weights( n_seqs );
	for( auto&& w_s: apegrunt::zip_range( int_weights, *alignment ) ) { get<0>(w_s) = get<1>(w_s)->multiplicity(); }

	const auto percent = ((n_seqs*(n_seqs-1))/2)/100; // integer division and type
	const auto output_frequency_mask = pow2_ceil_mask( percent*2 ); // make it a power-of-2 close to 2%
	std::atomic_size_t nonidentical_pairs(0);
	std::atomic_size_t pairs(0);

	auto print_statistics = [&]( std::size_t npairs ) {
		if( Apegrunt_options::verbose() )
		{
			*Apegrunt_options::get_out_stream() << "\rapegrunt: processing: " << std::setw(3) << std::setfill(' ') << std::size_t(double(npairs) / double(percent)) << "%";
			Apegrunt_options::get_out_stream()->flush();
		}
	};

	auto range_primer = [&]( const auto& pair ) {
		const auto seq_i = ali[0]; // (*alignment)[i];
		const auto seq_i_multiplicity = seq_i->multiplicity();

		for( auto first=pair.first; first != pair.second; ++first )
		{
			const auto j=first;
			if( seq_i->is_similar_to( *(ali[j]), Apegrunt_options::sample_reweighting_threshold()*double(L) ) )
			{
				int_weights[0] += ali[j]->multiplicity();
				int_weights[j] += seq_i_multiplicity;
			}
			else
			{
				++nonidentical_pairs;
			}
 			//++pairs;
			if( !(++pairs & output_frequency_mask) ) { print_statistics(pairs); }
 		}
	};

	auto range_worker = [&]( const auto& pair ) {
		for( auto first=pair.first; first != pair.second; ++first )
		{
			const auto i=first;
			const auto seq_i = ali[i];
			//const auto seq_i_size = seq_i->size();
			const auto seq_i_multiplicity = seq_i->multiplicity();
			for( std::size_t j = 0; j < i; ++j )
			{
				if( seq_i->is_similar_to( *(ali[j]), Apegrunt_options::sample_reweighting_threshold()*double(L) ) )
				{
					int_weights[i] += ali[j]->multiplicity();
					int_weights[j] += seq_i_multiplicity;
				}
				else
				{
					++nonidentical_pairs;
				}
				//++pairs
				if( !(++pairs & output_frequency_mask) ) { print_statistics(pairs); }
			}
		}
	};

	cputimer.start();

	print_statistics(pairs);

	// create *at most* as many tasks as there are sequences to process
#if __cplusplus < 201703L && !_OPENMP
	const std::size_t n(std::min(std::size_t(apegrunt::Apegrunt_options::threads()),n_seqs));
#else
	const std::size_t n(std::min(std::size_t(apegrunt::Apegrunt_options::threads()*2),n_seqs)); // if we have C++17 parallel algorithms, then twice as many tasks as there are threads
#endif

	// Caching done during the first call to is_similar_to() forces threads to wait for each other, which
	// skews the parallel work-load quite a bit. We mitigate this by priming the cache in a separate run.
	std::vector<std::pair<std::size_t,std::size_t> > primer_ranges; primer_ranges.reserve(n);
	auto first=0;
	for( auto i=0; i<n; ++i )
	{
		const auto last = std::size_t( double(n_seqs) * (double(i+1)/n) );
		primer_ranges.emplace_back( first, last );
		first = last;
	}
	assert(primer_ranges.back().second == n_seqs);

	// Partition the triangular matrix into (roughly) equally sized pieces.
	// This will help to balance (statically) the parallel work-load.
	std::vector<std::pair<std::size_t,std::size_t> > balanced_ranges; balanced_ranges.reserve(n);
	first=1;
	for( auto i=0; i<n; ++i )
	{
		const auto last = std::size_t( double(n_seqs) * std::sqrt(double(i+1)/n) );
		balanced_ranges.emplace_back( first, last );
		first = last;
	}
	assert(balanced_ranges.back().second == n_seqs);

	// The following evaluations are quite time-consuming, so run them in parallel, if possible
#if __cplusplus < 201703L
#ifdef _OPENMP
		omp_set_num_threads( apegrunt::Apegrunt_options::threads() );
		#pragma omp parallel for
		for( const auto& range: primer_ranges ) { range_primer(range); } // prime caches; compare sample 0 against all the others
		#pragma omp parallel for
		for( const auto& range: balanced_ranges ) { range_worker(range); } // compare samples [1,n_seq) against each other (the triangular part of a symmetric matrix)
#else // no _OPENMP
		//std::for_each(primer_ranges.begin(), primer_ranges.end(), range_primer);
		//std::for_each(balanced_ranges.begin(), balanced_ranges.end(), range_worker);
		{
			std::vector<std::thread> threads(n);
			using boost::get;
			for( auto t_r: apegrunt::zip_range(threads,primer_ranges) )
			{
				get<0>(t_r) = std::thread( range_primer, get<1>(t_r) );
			}
			for( auto& thread: threads ) { while( !thread.joinable() ); thread.join(); }

			for( auto t_r: apegrunt::zip_range(threads,balanced_ranges) )
			{
				get<0>(t_r) = std::thread( range_worker, get<1>(t_r) );
			}
			for( auto& thread: threads ) { while( !thread.joinable() ); thread.join(); }
		}
#endif // _OPENMP

#else // __cplusplus
	std::for_each(std::execution::par, primer_ranges.begin(), primer_ranges.end(), range_primer);
	std::for_each(std::execution::par, balanced_ranges.begin(), balanced_ranges.end(), range_worker);
#endif // __cplusplus

	print_statistics(pairs); // one last time for good measure and to get the counter to 100%

	cputimer.stop();
	if( Apegrunt_options::verbose() ) { cputimer.print_timing_stats(); }

	if( Apegrunt_options::verbose() )
	{
		const auto total_pairs = (n_seqs*(n_seqs-1))/2;
		*Apegrunt_options::get_out_stream()
			<< "apegrunt: # of non-identical pairs=" << nonidentical_pairs << " out of " << total_pairs << " (ratio=" << real_t(nonidentical_pairs)/real_t(total_pairs) << ")\n";
	}

	std::vector<real_t> weights; weights.reserve(n_seqs);

	using boost::get;
	for( auto&& seq_iw : zip_range(*alignment,int_weights) ) { weights.push_back( real_t( get<0>(seq_iw)->multiplicity() ) / real_t( get<1>(seq_iw) ) ); } // 1/weight * multiplicity
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
