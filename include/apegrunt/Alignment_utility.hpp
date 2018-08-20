/** @file Alignment_utility.hpp

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
#include "Threshold_rule.hpp"

#include "Alignment_forward.h"
#include "StateVector_forward.h"
#include "StateVector_state_types.hpp"
#include "StateVector_utility.hpp"
#include "Alignment_factory.hpp"

#include "Loci_parsers.hpp"
#include "State_block.hpp"

#include "aligned_allocator.hpp"

#include "accumulators/distribution_std.hpp"
#include "accumulators/distribution_bincount.hpp"
#include "accumulators/distribution_ordered.hpp"
#include "accumulators/distribution_cumulative.hpp"
//#include "accumulators/distribution_generator_svg.hpp"
#include "accumulators/distribution_generator_csv.hpp"

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

template< typename StateT, typename RealT=double >
std::vector< std::array< RealT, number_of_states<StateT>::N > > columnwise_frequencies( const Alignment_ptr<StateT> alignment )
{
	//std::cout << " {get column frequencies list"; std::cout.flush();
	using real_t = RealT;
	enum { N=number_of_states<StateT>::N };
	using boost::get;

	std::vector< std::array<real_t,N> > freq_vec( alignment->n_loci() ); //, {0.0} );
	const auto neff = real_t( alignment->size() ); // effective number of samples (sequences)

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
std::vector< std::array< RealT, number_of_states<StateT>::N > > weighted_columnwise_frequencies( const Alignment_ptr<StateT> alignment )
{
	using real_t = RealT;
	enum { N=number_of_states<StateT>::N };
	using boost::get;

	std::vector< std::array<real_t,N> > freq_vec( alignment->n_loci() ); //, {0.0} );
	const auto neff = real_t( alignment->effective_size() ); // effective number of samples (sequences)

	//std::cout << "apegrunt: neff=" << neff << std::endl;

	if( 0 < neff ) // protect against division by zero
	{
		const auto occ_vec = alignment->w_frequencies(); // weighted column biases

		for( auto&& occ_and_freq: zip_range(occ_vec, freq_vec) ) // loop over all loci
		{
			std::transform( cbegin(get<0>(occ_and_freq)), cend(get<0>(occ_and_freq)), begin(get<1>(occ_and_freq)), [=](auto occ){ return real_t(occ)/neff; } );
		}
	}

	return freq_vec;
}


// To use, inherit from extend_comparison_operators and
// define these public operators in your struct/class:
// inline bool T::operator==( const T& rhs ) const
// inline bool operator<( const my_type& rhs ) const
struct extend_comparison_operators { using enable_extended_comparison_operators = std::true_type; };

template< typename T, typename std::enable_if< std::is_same< typename T::enable_extended_comparison_operators,std::true_type >::type >::type >
inline bool operator!=( const T& lhs, const T& rhs ) { return !(lhs == rhs); }

template< typename T, typename std::enable_if< std::is_same< typename T::enable_extended_comparison_operators,std::true_type >::type >::type >
inline bool operator> ( const T& lhs, const T& rhs ) { return (rhs < lhs); }

template< typename T, typename std::enable_if< std::is_same< typename T::enable_extended_comparison_operators,std::true_type >::type >::type >
inline bool operator<=( const T& lhs, const T& rhs ) { return !(lhs > rhs); }

template< typename T, typename std::enable_if< std::is_same< typename T::enable_extended_comparison_operators,std::true_type >::type >::type >
inline bool operator>=( const T& lhs, const T& rhs ) { return !(lhs < rhs); }


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

	template< typename FilteredAlignmentT >
	Alignment_ptr< typename FilteredAlignmentT::state_t > operator()( const Alignment_ptr< typename FilteredAlignmentT::state_t > alignment ) const
	{
		if( Apegrunt_options::filter_alignment() )
		{
			//std::cout << "apegrunt: get filter list"; std::cout.flush();
			const auto accept_list = this->get_filter_list( alignment );
			//std::cout << " done" << std::endl;

			std::ostringstream id_stream;
			id_stream << "filtered";
			id_stream << "_ge" << std::setw(3) << std::setfill('0') << std::size_t(m_maf_threshold*100.) << "maf";
			id_stream << "_le" << std::setw(3) << std::setfill('0') << std::size_t(m_gap_threshold*100.) << "gf";
			id_stream << "_" << m_state_rule.safe_string() << "-lt04states";
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
/*
	template< typename ContainerT, typename T >
	struct back_inserter {
		back_inserter( ContainerT& container, const T& ) : m_container(container) { }
		template< typename T > T& operator()( const T& a ) { m_container.push_back(a.state)
	};
*/
	template< typename StateT >
	Loci_ptr get_filter_list( const Alignment_ptr<StateT> alignment ) const
	{
		//std::cout << " {get filter list"; std::cout.flush();
		using state_t = StateT;
		using std::begin; using std::end; using std::cbegin; using std::cend;

		std::array<std::size_t,number_of_states<state_t>::value> statef; for( auto& s: statef ) { s=0; }

		std::vector<index_and_state_count_association> proto_accept_list;
		std::vector<std::size_t> accept_list;
		//std::cout << "apegrunt: get columnwise frequencies"; std::cout.flush();
		const auto freq_vec = columnwise_frequencies( alignment );
		//std::cout << " done" << std::endl;

		std::size_t total_SNPs=0;
		std::size_t current_locus=0;

		for( const auto& freq: freq_vec ) // loop over all sequence positions
		{
			const std::size_t n = freq.size()-1; // number of states (excluding gaps); states are ordered such that the 'GAP' state is always the last

			std::size_t n_significant = 0;
			std::size_t n_nz = 0;
			std::size_t totaln_nz = 0;

			for( std::size_t i = 0; i < number_of_states<state_t>::value; ++i )
			{
				freq[i] > 0 && ++totaln_nz;
			}

			//std::cout << "[";
			for( std::size_t i = 0; i < n; ++i )
			{
				//std::cout << " " << freq[i];
				freq[i] > 0 && ++n_nz && m_maf_threshold <= freq[i] && ++n_significant; // mark state as significant if its frequency is at least maf_threshold
			}
			//std::cout << " ]\n";

			if( m_state_rule(n_significant) && n_nz < 4 && (freq[std::size_t( gap_state<state_t>::value )] <= m_gap_threshold) ) // && n_nz == n_significant )
			{
				accept_list.push_back(current_locus);
				//proto_accept_list.emplace_back(n_nz,current_locus);
			}

			n_nz > 1 && ++total_SNPs;

			statef[totaln_nz] += 1;

			++current_locus;
		}

		std::cout << "\n";
		for( std::size_t i=number_of_states<state_t>::value-1; i != 0; --i )
		{
			std::cout << "number of " << i+1 << "-state positions = " << statef[i] << std::endl;
		}
		//for( std::size_t i=0; i < statef.size(); ++i ) { std::cout << "apegrunt: " << i << " state = " << statef[i] << "\n"; }
/*
		accept_list.reserve(proto_accept_list.size());
		std::sort(begin(proto_accept_list),end(proto_accept_list));
		for( const auto& a: proto_accept_list ) { accept_list.push_back(a.locus_index); }
*/
		// move last element to first place
//		accept_list.insert(accept_list.begin(), accept_list[accept_list.size()-1]); accept_list.pop_back();

		//std::cout << "}"; std::cout.flush();
		return make_Loci_list(std::move(accept_list),0);
	}

};

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
bool output_sample_distance_matrix( const Alignment_ptr<StateT> alignment, std::ostream *os=nullptr )
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

template< typename AlignmentT, typename RealT >
Alignment_ptr< typename AlignmentT::state_t >
sequence_sample(
		const Alignment_ptr< typename AlignmentT::state_t >& alignment,
		RealT sample_fraction, // fraction [0,1] of sequences to keep
		std::size_t random_seed=0 // if random_seed==0, then generate a new random random seed
	)
{
	//using real_t = RealT;
	//using alignment_t = AlignmentT;
	//using state_t = typename alignment_t::state_t;

	// initialize random number generator
	if( 0 == random_seed )
	{
		std::random_device rd; // used only once -- generates random seed for mt generator below
		random_seed = rd();
	}
	std::mt19937_64 generator(random_seed);
	std::uniform_int_distribution<std::size_t> distribution( 0, alignment->size()-1 );
	auto grand = std::bind( distribution, generator );

	const std::size_t sample_count = alignment->size() * sample_fraction;

	auto accept_list = std::make_shared< std::vector<std::size_t> >();
	accept_list->reserve( sample_count );

	for( auto i=0; i < sample_count; ++i ) { accept_list->push_back( grand() ); }

	return Alignment_factory< AlignmentT >().copy_selected( alignment, std::move(accept_list) );
}

namespace detail
{
template< typename StateT, typename RealT >
struct state_frequency_pair
{
	using state_t = StateT;
	using real_t = RealT;
	using my_type = state_frequency_pair<state_t,real_t>;
	RealT freq;
	StateT state;

	state_frequency_pair( StateT s, RealT f ) : freq(f), state(s) { }

	state_frequency_pair( const my_type& rhs ) : freq(rhs.freq), state(rhs.state) { }
	state_frequency_pair( my_type&& rhs ) noexcept : freq(std::move(rhs.freq)), state(std::move(rhs.state)) { }

	my_type& operator=( const my_type& rhs ) { freq=rhs.freq; state=rhs.state; return *this; }
	my_type& operator=( my_type&& rhs ) noexcept { freq=std::move(rhs.freq); state=std::move(rhs.state); return *this; }

	bool operator==( const my_type& rhs ) const { return freq == rhs.freq; }
	bool operator<( const my_type& rhs ) const { return freq < rhs.freq; }
};

template< typename StateT, typename RealT >
inline bool operator!=( const state_frequency_pair<StateT,RealT>& lhs, const state_frequency_pair<StateT,RealT>& rhs ) { return !(lhs == rhs); }
template< typename StateT, typename RealT >
inline bool operator>( const state_frequency_pair<StateT,RealT>& lhs, const state_frequency_pair<StateT,RealT>& rhs ) { return (rhs < lhs); }
template< typename StateT, typename RealT >
inline bool operator<=( const state_frequency_pair<StateT,RealT>& lhs, const state_frequency_pair<StateT,RealT>& rhs ) { return !(lhs > rhs); }
template< typename StateT, typename RealT >
inline bool operator>=( const state_frequency_pair<StateT,RealT>& lhs, const state_frequency_pair<StateT,RealT>& rhs ) { return !(lhs < rhs); }

} // namespace detail

template< typename StateT >
Loci_ptr maf_based_alignment_reordering( const Alignment_ptr<StateT> alignment )
{
	//using state_t = StateT;
	using locus_frequency_pair = detail::state_frequency_pair< std::size_t, double >;
	using std::begin; using::std::end;

	std::vector< locus_frequency_pair > locus_reorder_list;
	const auto freq_vec = allele_frequency( alignment );
	std::vector<std::size_t> new_ordering; new_ordering.reserve(freq_vec.size());
	std::size_t current_locus=0;

	for( const auto& freq: freq_vec ) // loop over all sequence positions
	{
		locus_reorder_list.emplace_back( current_locus, freq[0] );
		++current_locus;
	}

	std::sort( begin(locus_reorder_list), end(locus_reorder_list), std::greater<locus_frequency_pair>() );
	for( const auto& lfpair: locus_reorder_list )
	{
		//std::cout << "idx=" << lfpair.state << " f=" << lfpair.freq << std::endl;
		new_ordering.push_back( lfpair.state );
	}

	return make_Loci_list(std::move(new_ordering),0);
}

template< typename StateT >
Loci_ptr state_flip_based_alignment_reordering( const Alignment_ptr<StateT> alignment )
{
	using locus_frequency_pair = detail::state_frequency_pair< std::size_t, std::size_t >;
	using std::begin; using::std::end;

	const std::size_t nloci = alignment->n_loci();
	std::vector<std::size_t> new_ordering; new_ordering.reserve(nloci);
	std::vector< locus_frequency_pair > locus_reorder_list; locus_reorder_list.reserve(nloci);

	for( std::size_t i=0; i < nloci; ++i )
	{
		bool first = true;
		State_holder<StateT> previous;
		std::size_t nflips = 0;
		for( const auto sequence: alignment )
		{
			const auto state = (*sequence)[i];
			if( first ) { previous = state; first = false; }
			else { if( state != previous ) { ++nflips; } }
		}
		locus_reorder_list.emplace_back( i, nflips );
	}
	std::sort( begin(locus_reorder_list), end(locus_reorder_list), std::greater<locus_frequency_pair>() );
	for( const auto& flippair: locus_reorder_list )
	{
		//std::cout << "idx=" << lfpair.state << " f=" << lfpair.freq << std::endl;
		new_ordering.push_back( flippair.state );
	}

	return make_Loci_list(std::move(new_ordering),0);
}

template< typename StateT >
Loci_ptr get_alignment_reordering( const Alignment_ptr<StateT> alignment )
{
	using locus_frequency_pair = detail::state_frequency_pair< std::size_t, std::size_t >;
	using std::begin; using::std::end;

	const std::size_t nloci = alignment->n_loci();
	std::vector<std::size_t> new_ordering; new_ordering.reserve(nloci);
	std::vector< locus_frequency_pair > locus_reorder_list; locus_reorder_list.reserve(nloci);

	for( std::size_t i=0; i < nloci; ++i )
	{
		bool first = true;
		State_holder<StateT> previous;
		std::size_t nflips = 0;
		for( const auto sequence: alignment )
		{
			const auto state = (*sequence)[i];
			if( first ) { previous = state; first = false; }
			else { if( state != previous ) { ++nflips; } }
		}
		locus_reorder_list.emplace_back( i, nflips );
	}
	std::sort( begin(locus_reorder_list), end(locus_reorder_list), std::greater<locus_frequency_pair>() );
	for( const auto& flippair: locus_reorder_list )
	{
		//std::cout << "idx=" << lfpair.state << " f=" << lfpair.freq << std::endl;
		new_ordering.push_back( flippair.state );
	}

	return make_Loci_list(std::move(new_ordering),0);
}

template< typename DestStateT, typename SrcStateT >
std::vector< std::vector< State_holder<SrcStateT> > > get_transform_list( const Alignment_ptr<SrcStateT> alignment, Loci_ptr& loci_list )
{
	using oldstate_t = SrcStateT;
	enum { Nold=number_of_states<oldstate_t>::N };
	using state_t = DestStateT;
	enum { N=number_of_states<state_t>::N };
	using sf_pair_t = detail::state_frequency_pair< State_holder<oldstate_t>, double >;

	const auto maf_threshold = Apegrunt_options::get_minor_allele_frequency_threshold();

	//const auto N = number_of_states<state_t>::N-1;
	//const auto Nold = number_of_states<oldstate_t>::N-1;

	std::vector<std::size_t> accept_list;
	std::vector< std::vector< State_holder<oldstate_t> > > transform_list;

	const auto freq_vec = columnwise_frequencies( alignment );

	std::size_t current_locus=0;
	std::size_t n_dropped = 0;
	std::size_t n_loci_with_dropped_alleles = 0;

//	std::array< state_t, number_of_states<state_t>::value > state_map;

	//std::vector<std::size_t> dropped_positions;

	//std::size_t locus=0;
	for( const auto& freq: freq_vec ) // loop over all sequence positions
	{
/*
		std::map<state_t,double> statecache;
		for( std::size_t i = 0; i < Nold-1; ++i )
		{
			if( freq[i] > 0 )
			{
				statecache[i].freq = freq[i];
			}
		}
*/
		//std::array< sf_pair_t, number_of_states<oldstate_t>::N > statecache = {{ {oldstate_t::a,0}, {oldstate_t::t,0}, {oldstate_t::c,0}, {oldstate_t::g,0}, {oldstate_t::GAP,0} }};
		std::array< sf_pair_t, Nold-1 > statecache = {{ {oldstate_t::a,freq[0]}, {oldstate_t::t,freq[1]}, {oldstate_t::c,freq[2]}, {oldstate_t::g,freq[3]} }};
		//for( std::size_t i = 0; i < Nold-1; ++i ) { statecache[i].freq = freq[i];/* std::cout << " " << State_holder<oldstate_t>(statecache[i].state) << "=" << freq[i]*100;*/ }/* std::cout << std::endl;*/

		using std::begin; using::std::end;
		std::sort( begin(statecache), end(statecache), std::greater<sf_pair_t>() );
/*
		std::cout << "post-sorting: ";
		for( std::size_t i = 0; i < Nold; ++i ) { std::cout << " \"" << State_holder<oldstate_t>(statecache[i].state) << "\"=" << statecache[i].freq*100; } std::cout << std::endl;
*/
		//if( statecache[Nold].freq > 0 ) { accept_list[current_locus].push_back( oldstate_t::GAP ); }


		const auto n_dropped_old = n_dropped;
		for( std::size_t i = N-1; i < Nold-1; ++i ) { statecache[i].freq > 0 && ++n_dropped; }

		if( n_dropped_old != n_dropped )
		{
			++n_loci_with_dropped_alleles;
			if( maf_threshold <= statecache[N-1].freq )
			{
				//dropped_positions.push_back( (*alignment->get_loci_translation())[transform_list.size()-1]+Apegrunt_options::get_input_indexing_base() );
				if( Apegrunt_options::verbose() )
				{
					auto out = Apegrunt_options::get_out_stream();
					*out << "apegrunt: Found more significant alleles (" << N-1+(n_dropped-n_dropped_old) << ") than we can store (" << N-1 << "); drop position " << transform_list.size() << " (original locus " << (*alignment->get_loci_translation())[transform_list.size()-1]+Apegrunt_options::get_input_indexing_base() << ") due to: ";
					out->precision(2);
					for( std::size_t i = N-1; i < Nold-1; ++i ) { *out << " \"" << statecache[i].state << "\"=" << statecache[i].freq*100; }
					*out << " in [";
					for( const auto & sc: statecache ) { *out << " \"" << sc.state << "\"=" << sc.freq*100; }
					//for(std::size_t i = 0; i < Nold; ++i ) { *out << " \"" << statecache[i].state << "\"(" << i << ")=" << statecache[i].freq*100; }
					*out << " ]\n";
				}
			}
		}
		else
		{
			transform_list.emplace_back(); // add empty element
			for( std::size_t i = 0; i < N-1; ++i )
			{
				if( statecache[i].freq > 0 )
				{
					transform_list.back().push_back( statecache[i].state );
					//accept_list[current_locus].push_back( statecache[i].state );
				}
			}
			accept_list.push_back( current_locus );
		}
		++current_locus;
	}

	if( n_dropped > 0 )
	{
		if( Apegrunt_options::verbose() )
		{
			*Apegrunt_options::get_out_stream() << "apegrunt: NOTE: " << n_loci_with_dropped_alleles << " loci will be dropped (" << n_dropped << " alleles in total)" << std::endl;
		}
	/*
		auto dropped_pos_file = apegrunt::get_unique_ofstream( alignment->id_string()+"."+apegrunt::size_string(alignment)+".dropped" );

		if( Apegrunt_options::verbose() )
		{
			*Apegrunt_options::get_out_stream() << "apegrunt: write dropped positions to file " << dropped_pos_file->name() << "\n";
		}
		auto out = dropped_pos_file->stream();
		for( const auto pos: dropped_positions )
		{
			*out << pos << "\n";
		}
	*/
	}

	//std::cout << "apegrunt: accept list has " << accept_list.size() << " entries." << std::endl;

	loci_list = make_Loci_list(std::move(accept_list),0);

	return transform_list;
}

template< typename DestStateT, typename SrcStateT >
Alignment_ptr<DestStateT> transform_alignment( const Alignment_ptr<SrcStateT> alignment )
{
	using alignment_t = Alignment_impl_block_compressed_storage< StateVector_impl_block_compressed_alignment_storage<SrcStateT> >;
	auto transformable_alignment = alignment;
	Loci_ptr accept_list;
	auto transform_list = get_transform_list<DestStateT>(alignment, accept_list);
	if( accept_list ) { transformable_alignment = Alignment_factory<alignment_t>()( alignment, accept_list ); }

	return Alignment_factory< Alignment_impl_block_compressed_storage< StateVector_impl_block_compressed_alignment_storage<DestStateT> > >()( transformable_alignment, std::move(transform_list) );
}

template< typename StateT >
Alignment_ptr<StateT> transpose_alignment( const Alignment_ptr<StateT> alignment )
{
	using alignment_t = Alignment_impl_block_compressed_storage< StateVector_impl_block_compressed_alignment_storage<StateT> >;

	return Alignment_factory< Alignment_impl_block_compressed_storage< StateVector_impl_block_compressed_alignment_storage<StateT> > >().transpose( alignment );
}



} // namespace apegrunt

#endif // APEGRUNT_ALIGNMENT_UTILITY_HPP
