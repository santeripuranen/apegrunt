/** @file Alignment_filter.hpp

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
#ifndef APEGRUNT_ALIGNMENT_FILTER_HPP
#define APEGRUNT_ALIGNMENT_FILTER_HPP

#include <vector>
#include <array>
#include <string>
#include <algorithm>

#include "Alignment_forward.h"
#include "StateVector_state_types.hpp"
#include "Apegrunt_options.h"
#include "Loci_parsers.hpp"

namespace apegrunt {

template< typename StateT, typename RealT=double >
std::vector< std::array< RealT, number_of_states<StateT>::N > > columnwise_frequencies( const Alignment<StateT>* alignment )
{
	//std::cout << " {get column frequencies list"; std::cout.flush();
	using real_t = RealT;
	enum { N=number_of_states<StateT>::N };
	using boost::get;

	std::vector< std::array<real_t,N> > freq_vec( alignment->n_loci() ); //, {0.0} );
	const auto neff = real_t( alignment->size() ); // number of samples (sequences)

	if( 0 < neff ) // protect against division by zero
	{
		const auto occ_vec = alignment->frequencies(); // unweighted column biases

		for( auto&& occ_and_freq: zip_range(occ_vec, freq_vec) ) // loop over all loci
		{
			std::transform( cbegin(get<0>(occ_and_freq)), cend(get<0>(occ_and_freq)), begin(get<1>(occ_and_freq)), [=](auto occ){ return real_t(occ)/neff; } );
		}
	}
	//std::cout << "}"; std::cout.flush();

	return freq_vec;
}

template< typename StateT, typename RealT=double >
std::vector< std::array< RealT, number_of_states<StateT>::N > > columnwise_frequencies( const Alignment_ptr<StateT> alignment )
{
	return apegrunt::columnwise_frequencies( alignment.get() );
}

template< typename StateT, typename RealT=double >
std::vector< std::array< RealT, number_of_states<StateT>::N > > weighted_columnwise_frequencies( const Alignment<StateT>* alignment, RealT pseudocount=0.0 )
{
	using real_t = RealT;
	enum { N=number_of_states<StateT>::N };
	using boost::get;

	std::vector< std::array<real_t,N> > freq_vec( alignment->n_loci() ); //, {0.0} );
	const auto neff = real_t( alignment->effective_size() ); // + real_t(N)*0.5; // effective number of samples (sequences)

	//std::cout << "apegrunt: neff=" << neff << std::endl; std::cout.flush();

	if( 0 < neff ) // protect against division by zero
	{
		//std::cout << "apegrunt: get w_frequencies"; std::cout.flush();
		const auto occ_vec = alignment->w_frequencies(); // weighted column biases

		//std::cout << "apegrunt: transform frequencies"; std::cout.flush();
		for( auto&& occ_and_freq: apegrunt::zip_range(occ_vec, freq_vec) ) // loop over all loci
		{
			std::transform( cbegin(get<0>(occ_and_freq)), cend(get<0>(occ_and_freq)), begin(get<1>(occ_and_freq)), [=](auto occ){ return (real_t(occ)+pseudocount)/(neff+pseudocount*N); } );
/*
			auto sum = std::accumulate( begin(get<1>(occ_and_freq)), end(get<1>(occ_and_freq)), real_t(0), [](auto sum, auto f){return sum+f;});
			if( std::abs(1-sum) > 1e-10 )
			{
				for( auto f: get<0>(occ_and_freq) ) { std::cout << " " << f; }
				std::cout << " sum=" << sum;
				for( auto f: get<1>(occ_and_freq) ) { std::cout << " " << f; }
				std::cout << std::endl;
				exit(0);
			}
*/
		}
	}

	//std::cout << "apegrunt: return " << freq_vec.size() << "-by-" << N << " frequencies" << std::endl; std::cout.flush();

	return freq_vec;
}

template< typename StateT, typename RealT=double >
std::vector< std::array< RealT, number_of_states<StateT>::N > > weighted_columnwise_frequencies( const Alignment_ptr<StateT> alignment, RealT pseudocount=0.0 )
{
	return apegrunt::weighted_columnwise_frequencies( alignment.get(), pseudocount );
}

template< typename StateT, typename RealT=double >
std::vector< RealT > weighted_column_entropies( const Alignment<StateT>* alignment, RealT pseudocount=0.0 )
{
	using real_t = RealT;
	using std::cbegin;
	enum { N=number_of_states<StateT>::N };

	//std::cout << "apegrunt: get weighted columnwise frequencies" << std::endl; std::cout.flush();
	auto frequencies = apegrunt::weighted_columnwise_frequencies( alignment, pseudocount );
	//std::cout << "apegrunt: get statecounts" << std::endl; std::cout.flush();
	auto statecounts = alignment->get_statecount_blocks();

	std::vector< real_t > entropies; entropies.reserve( alignment->n_loci() );

	//auto statecount_itr = cbegin(statecounts);
	//std::size_t i=0;
	//std::cout << "apegrunt: calculate entropies" << std::endl; std::cout.flush();
	for( auto& freq: frequencies )
	{
		// entropy
		entropies.push_back( -apegrunt::guard_sum_xlogx( apegrunt::make_Vector_view<real_t,N>( freq.data() ) ) );

		// normalized entropy
/*
		const std::size_t count = (*statecount_itr)[i];
		entropies.push_back( apegrunt::guard_sum_xlogx( apegrunt::make_Vector_view<real_t,N>( freq.data() ) ) / (apegrunt::xlogx(1.0/real_t(count))*real_t(count)) );

		if( ++i % apegrunt::StateBlock_size == 0 ) { i=0; ++statecount_itr; }


		for( const auto f: freq ) { std::cout << f << " "; }
		std::cout << entropies.back() << std::endl;
*/
	}
	//std::cout << "apegrunt: return entropies" << std::endl; std::cout.flush();
	return entropies;
}

template< typename StateT, typename RealT=double >
std::vector< RealT > weighted_column_entropies( const Alignment_ptr<StateT> alignment, RealT pseudocount=0.0 )
{
	return apegrunt::weighted_column_entropies( alignment.get(), pseudocount );
}

// primary template matches all state types
template< typename StateT, typename HasGapState = void >
struct Gap_rule
{
	using state_t = StateT;

	Gap_rule() = delete;
	Gap_rule( double threshold=1 ) { }
	~Gap_rule() = default;

	// always returns true
	inline bool operator()( const double* const val ) const { return true; }
};

// partial specialization for all state types that contain an actual gap state
template< typename StateT >
struct Gap_rule< StateT, typename std::enable_if<apegrunt::has_gap_state<StateT>::value>::type >
{
	using state_t = StateT;

	Gap_rule() = delete;
	Gap_rule( double threshold=1 ) : m_threshold(threshold) { }
	~Gap_rule() = default;

	inline bool operator()( const double* const val ) const { return *(val+int(gap_state<state_t>::value)) <= m_threshold; }
	//inline bool operator()( const double* const val ) const { return std::max( *(val+int(gap_state<state_t>::value)), *(val+int('N')) ) <= m_threshold; }

	const double m_threshold;
};

class Alignment_filter
{

public:

	enum ParameterPolicy
	{
		AQUIRE_GLOBAL=0,
		FILTER_SNPS=1,
		DEFERRED_AQUIRE=2
	};

	Alignment_filter( ParameterPolicy parameter_policy=ParameterPolicy::FILTER_SNPS )
	{
		switch( parameter_policy )
		{
		case ParameterPolicy::AQUIRE_GLOBAL:
			m_gap_threshold = Apegrunt_options::get_gap_frequency_threshold();
			m_maf_threshold = Apegrunt_options::get_minor_allele_frequency_threshold();
			m_state_rule = Apegrunt_options::get_allele_state_rule();
			break;
		case ParameterPolicy::FILTER_SNPS:
			m_gap_threshold = 1.0;
			m_maf_threshold = 0.0;
			m_state_rule = Threshold_rule<int>( std::string(">1") );
			break;
		default:
			m_gap_threshold = 1.0;
			m_maf_threshold = 0.0;
			m_state_rule = Threshold_rule<int>( std::string(">0") );
		}
	}
	~Alignment_filter() = default;
/*
	template< typename FilteredAlignmentT >
	Alignment_ptr< typename FilteredAlignmentT::state_t > operator()( const Alignment_ptr< typename FilteredAlignmentT::state_t > alignment ) const
	{
		using state_t = typename FilteredAlignmentT::state_t;
		if( Apegrunt_options::filter_alignment() )
		{
			//std::cout << "apegrunt: get filter list"; std::cout.flush();
			const auto accept_list = this->get_filter_list( alignment );
			//std::cout << " done" << std::endl;

			std::ostringstream id_stream;
			id_stream << "filtered";
			id_stream << "_ge" << std::setw(3) << std::setfill('0') << std::size_t(m_maf_threshold*100.) << "maf";
			id_stream << "_le" << std::setw(3) << std::setfill('0') << std::size_t(m_gap_threshold*100.) << "gf";
			id_stream << "_" << m_state_rule.safe_string() << "-lt0" << apegrunt::number_of_states<state_t>::value << "states";
			//id_stream << "_L" << accept_list->size() << "n" << alignment->size();

			accept_list->set_id_string( id_stream.str() );

			//std::cout << "apegrunt: create new alignment"; std::cout.flush();
			return Alignment_factory< FilteredAlignmentT >()( alignment, accept_list );
		}
		else
		{
			return alignment;
		}
	}
*/

	template< typename StateT >
	Loci_ptr get_filter_list( const Alignment<StateT>* alignment ) const
	{
		//std::cout << "Get filter list" << std::endl;
		using state_t = StateT;
		using std::begin; using std::end; using std::cbegin; using std::cend;

		//std::array<std::size_t,number_of_states<state_t>::value+1> statef; for( auto& s: statef ) { s=0; }

		std::vector<std::size_t> accept_list;

		//std::vector< std::pair<std::size_t,std::size_t> > accept;

		//std::cout << "apegrunt: get columnwise state frequencies"; std::cout.flush();
		const auto freq_vec = columnwise_frequencies( alignment );
		//std::cout << " done" << std::endl;

		//std::cout << "apegrunt: get gap rule"; std::cout.flush();
		const auto gap_rule = Gap_rule<state_t>( m_gap_threshold );
		//std::cout << " done" << std::endl;

		//std::size_t total_SNPs=0;
		std::size_t current_locus=0;

		// number of states (excluding gaps); states are ordered such that the 'GAP' state is always the last
		//const std::size_t n = apegrunt::number_of_states<state_t>::value - ( apegrunt::has_gap_state<state_t>::value ? 1 : 0 );
		const auto n = apegrunt::number_of_states<state_t>::value - std::size_t(apegrunt::has_gap_state<state_t>::value); // the latter is a bool
		//const std::size_t n = apegrunt::number_of_states<state_t>::value;

		//std::cout << "#states=" << number_of_states<state_t>::value << " | n=" << n << " gap=" << std::size_t(apegrunt::gap_state<state_t>::value) << std::endl;

		//std::cout << "\nGAP=" << std::size_t(apegrunt::gap_state<state_t>::value) << " (" << to_char(apegrunt::gap_state<state_t>::value) << ")\n";

		//std::cout << "apegrunt: scan frequencies"; std::cout.flush();
		for( const auto& freq: freq_vec ) // loop over all sequence positions
		{
			//const std::size_t n = freq.size() - ( apegrunt::has_gap_state<state_t>::value ? 1 : 0 ); // number of states (excluding gaps); states are ordered such that the 'GAP' state is always the last
			//std::cout << "#" << current_locus;
			//for( std::size_t i = 0; i < n; ++i ) { if(freq[i] != 0) { if(std::size_t(apegrunt::gap_state<state_t>::value) != i) { std::cout << "  " << freq[i]; } else {std::cout << " *" << freq[i]; } } }
			//std::cout << std::endl;

			std::size_t n_significant = 0;
			//std::size_t n_nz = 0;

			for( std::size_t i = 0; i < n; ++i )
			{	// exclude gaps
				freq[i] >= m_maf_threshold && ++n_significant; // mark state as significant if its frequency is at least maf_threshold
				//std::cout << "freq[" << i << "]="; std::cout.flush();
				//freq[i] > 0 && ++n_nz && m_maf_threshold <= freq[i] && ++n_significant; // mark state as significant if its frequency is at least maf_threshold
				//freq[i] > 0 && std::size_t(apegrunt::gap_state<state_t>::value) != i && /*std::size_t('N') != i &&*/ /* ++n_nz && */ !(m_maf_threshold > freq[i]) && ++n_significant; // mark state as significant if its frequency is at least maf_threshold
				//std::cout << freq[i] << std::endl;
			}
			//if( m_state_rule(n_significant) && n_nz < apegrunt::number_of_states<state_t>::value && (freq[std::size_t( gap_state<state_t>::value )] <= m_gap_threshold) ) // && n_nz == n_significant )
			//if( m_state_rule(n_significant) && (freq[std::size_t( gap_state<state_t>::value )] <= m_gap_threshold) ) // && n_nz == n_significant )
			//std::cout << " | " << n_significant; // << std::endl;
			if( m_state_rule(n_significant) && gap_rule( freq.data() ) ) // && n_nz == n_significant )
			{
				//std::cout << " -> ACCEPT (" << current_locus << ")";
				accept_list.push_back(current_locus);
				//accept.emplace_back(n_nz,current_locus);
			}
			//else { std::cout << " -> DROP"; }

			//std::cout << std::endl;

			//n_nz > 1 && ++total_SNPs;

			++current_locus;
		}
		//std::cout << " done" << std::endl;
		//std::sort( begin(accept), end(accept), std::less<>() );
		//accept_list.resize(accept.size());
		//std::transform( begin(accept), end(accept), begin(accept_list), [](auto pair){ return pair.second; } );


		auto filter_list = apegrunt::make_Loci_list(std::move(accept_list),0);

		std::ostringstream id_stream;
		id_stream << "filtered";
		id_stream << "_ge" << std::setw(3) << std::setfill('0') << std::size_t(m_maf_threshold*100.) << "maf";
		id_stream << "_le" << std::setw(3) << std::setfill('0') << std::size_t(m_gap_threshold*100.) << "gf";
		id_stream << "_" << m_state_rule.safe_string() << "-lt0" << apegrunt::number_of_states<state_t>::value << "states";
		//id_stream << "_L" << accept_list->size() << "n" << alignment->size();

		filter_list->set_id_string( id_stream.str() );

		//std::cout << " done" << std::endl;

		return filter_list;
	}

	template< typename StateT >
	Loci_ptr get_filter_list( const Alignment_ptr<StateT> alignment ) const
	{
		return this->get_filter_list( alignment.get() );
	}

	template< typename StateT >
	Loci_ptr get_entropy_filter_list( const Alignment<StateT>* alignment, double entropy_threshold, double pseudocount=0.0 ) const
	{
		//std::cout << "Get filter list" << std::endl;
		using std::begin; using std::end; using std::cbegin; using std::cend;

		std::vector<std::size_t> accept_list;

		using std::begin; using std::end;
		auto entropies = apegrunt::weighted_column_entropies( alignment, pseudocount );
		//const auto eend = std::stable_partition( begin(entropies), end(entropies), [&entropy_threshold](double e){ return e > entropy_threshold; } );
		//for( auto i=0; i< std::distance( begin(entropies), eend ); ++i ) { accept_list.push_back(i); }
		for( auto i=0; i<entropies.size(); ++i ) { if( entropies[i] > entropy_threshold ) { accept_list.push_back(i); } }
		auto filter_list = apegrunt::make_Loci_list(std::move(accept_list),0);

		std::ostringstream id_stream;
		id_stream << "entropy_gt" << std::setw(6) << std::setfill('0') << std::size_t(entropy_threshold*apegrunt::ipow(10,6));

		filter_list->set_id_string( id_stream.str() );

		//std::cout << " done" << std::endl;

		return filter_list;
	}

	template< typename StateT >
	Loci_ptr get_entropy_filter_list( const Alignment_ptr<StateT> alignment, double entropy_threshold, double pseudocount=0.0 ) const
	{
		return this->get_entropy_filter_list( alignment.get(), entropy_threshold, pseudocount );
	}


private:
	double m_gap_threshold;
	double m_maf_threshold;
	Threshold_rule<int> m_state_rule;

	struct index_and_state_count_association : public extend_comparison_operators
	{
		using my_type = index_and_state_count_association;

		index_and_state_count_association( const std::size_t& count, const std::size_t& index ) : state_count(count), locus_index(index) { }

		index_and_state_count_association( const my_type& other ) : state_count(other.state_count), locus_index(other.locus_index) { }

		inline bool operator==( const my_type& rhs ) const
		{
			return state_count == rhs.state_count;
		}

		inline bool operator<( const my_type& rhs ) const
		{
			return state_count < rhs.state_count;
		}

		std::size_t state_count;
		std::size_t locus_index;
	};

	template< typename StateT >
	std::array<std::size_t,number_of_states<StateT>::value+1> get_state_statistics( const Alignment_ptr<StateT> alignment ) const
	{
		using state_t = StateT;
		std::array<std::size_t,number_of_states<state_t>::value+1> statef; for( auto& s: statef ) { s=0; }
		const auto freq_vec = columnwise_frequencies( alignment );

		for( const auto& freq: freq_vec ) // loop over all sequence positions
		{
			std::size_t totaln_nz = 0;

			for( std::size_t i = 0; i < number_of_states<state_t>::value; ++i )
			{
				freq[i] > 0 && ++totaln_nz;
			}

			statef[totaln_nz] += 1;
		}

		return statef;
	}

	template< int N >
	void print_state_statistics( std::array<std::size_t,N> state_statistics, std::ostream *out=nullptr )
	{
		std::size_t positions = 0;
		std::size_t SNPs = 0;
		*out << "\n";
		for( std::size_t i=N; i != 0; --i )
		{
			*out << i << "-state positions = " << state_statistics[i] << std::endl;
			positions += state_statistics[i];
			i > 1 && SNPs += state_statistics[i];
		}
		*out << "SNPs / total positions: " << SNPs << " / " << positions << std::endl;
	}

};

} // namespace APEGRUNT_ALIGNMENT_FILTER_HPP

#endif // APEGRUNT_ALIGNMENT_FILTER_HPP
