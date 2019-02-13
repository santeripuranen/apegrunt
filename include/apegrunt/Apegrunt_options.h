	/** @file Apegrunt_options.h

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

#ifndef APEGRUNT_OPTIONS_H
#define APEGRUNT_OPTIONS_H

#include <iosfwd>
#include <vector>
#include <string>

// Boost includes
#include <boost/program_options.hpp>

#include "Threshold_rule.hpp"

namespace po = boost::program_options;

namespace apegrunt {

class Apegrunt_options_base
{
public:
	Apegrunt_options_base();
	Apegrunt_options_base( std::ostream *out, std::ostream *err=nullptr );
	virtual ~Apegrunt_options_base();

	void AddOptions( po::options_description *opdesc );
#ifdef APEGRUNT_STANDALONE_BUILD
	void AddOptions_standalone( po::options_description *opdesc );
#endif // APEGRUNT_STANDALONE_BUILD
	static bool CheckOptions( po::variables_map *varmap );
	static const std::string& s_get_copyright_notice_string();

	// Apegrunt standalone binary options
	static const std::string& s_get_usage_string();

	static const std::string& s_get_version_string();
	static const std::string& s_get_title_string();

	//> Test if textual output is desired. If true, then a call to get_out_stream() is guaranteed to return a valid (as in != null_ptr) ostream*.
	static bool verbose();
	static void set_verbose( bool value );

#ifdef APEGRUNT_STANDALONE_BUILD
	static const std::string& get_options_string();
	static void set_options_string( const std::string& options_string );
#endif // APEGRUNT_STANDALONE_BUILD

	//> Set an ostream. An invalid ostream* (as in "out->good() == false", will reset internal ostream ("ostream* == null_ptr").
	static void set_out_stream( std::ostream* out );
	static void set_err_stream( std::ostream* err );
	static std::ostream* get_out_stream();
	static std::ostream* get_err_stream();

	// Alignment IO options
	// Input
	static bool has_alignment_filenames();
	static const std::vector< std::string >& get_alignment_filenames();

	static bool has_includelist_filename();
	static const std::string& get_includelist_filename();

	static bool has_excludelist_filename();
	static const std::string& get_excludelist_filename();

	static bool has_mappinglist_filename();
	static const std::string& get_mappinglist_filename();

	static bool has_samplelist_filename();
	static const std::string& get_samplelist_filename();

	static bool has_sample_weights_filename();
	static const std::string& get_sample_weights_filename();

	static std::size_t get_input_indexing_base();
	static void set_input_indexing_base( std::size_t base_index );

	// Output control
	static bool output_state_frequencies();
	static bool output_sample_weights();
	static bool output_sample_distance_matrix();

	// Alignment processing
	static bool filter_alignment();
	static bool reweight_samples();
	static bool rescale_sample_weights();
	static int sample_alignment();

	static double get_minor_allele_frequency_threshold();
	static double get_gap_frequency_threshold();
	static Threshold_rule<int> get_allele_state_rule();
	static double sample_reweighting_threshold();
	static const std::vector<double>& sample_fractions();

	static std::size_t genome_size();
	static bool linear_genome();

	static std::size_t get_output_indexing_base();
	static void set_output_indexing_base( std::size_t base_index );

	static bool cuda();
	static void set_cuda( bool use_cuda=true );

	static int threads();
	static void set_threads( int nthreads );

	// Unused options
	static bool fuse_duplicates(); // functionality unused at the moment
	static bool optimize_column_order(); // not used at the moment

	static bool variable_penalty(); // not used at the moment
	static std::size_t distance_penalty_threshold(); // not used at the moment
	static double distance_penalty_shape(); // not used at the moment
	static double distance_penalty_scaling(); // not used at the moment

	static int get_begin_locus();
	static int get_end_locus();

private:
	static std::vector< std::string > s_alignment_file_names;
	static std::string s_includelist_file_name;
	static std::string s_excludelist_file_name;
	static std::string s_mappinglist_file_name;
	static std::string s_samplelist_file_name;
	static std::string s_sample_weights_file_name;

	//	static bool s_fuse_duplicates;
	static double s_minor_allele_frequency_threshold;
	static double s_gap_frequency_threshold;
	static Threshold_rule<int> s_allele_state_rule;
	static bool s_no_filter_alignment;
	static bool s_output_state_frequencies;
	static bool s_no_optimize_column_order;
	static int s_sample_alignment;
	static double s_sample_reweighting_threshold;
	static bool s_output_sample_weights;
	static bool s_no_sample_reweighting;
	static bool s_rescale_sample_weights;
	static bool s_output_sample_distance_matrix;

	static std::size_t s_genome_size;
	static bool s_linear_genome;
	static bool s_variable_penalty;
	static std::size_t s_distance_penalty_threshold;
	static double s_distance_penalty_shape;
	static double s_distance_penalty_scaling;

	static std::ostream *s_out;
	static std::ostream *s_err;
//	static std::string s_errfile_name;
//	static std::string s_outfile_name;
	static std::size_t s_output_indexing_base;
	static std::size_t s_input_indexing_base;
	static int s_begin_locus;
	static int s_end_locus;

	static int s_threads;
	static int s_nodes;
	static bool s_use_cuda;
	static bool s_verbose;

	static const std::string s_title_string;
	static const std::string s_version_string;
	static const std::string s_copyright_notice;
	static const std::string s_long_copyright_notice;
	static std::string s_options_string;
	static const std::string s_usage_string;
	void m_init();

	static void s_init_verbose( bool verbose );
#ifndef APEGRUNT_NO_TBB // Threading with Threading Building Blocks
	static void s_init_threads( int nthreads );
#endif // APEGRUNT_NO_TBB

#ifndef APEGRUNT_NO_CUDA
	static void s_init_use_cuda( bool use_cuda );
#endif // APEGRUNT_NO_CUDA
	static void s_init_fuse_duplicates( bool flag );
	static void s_init_minor_allele_frequency_threshold( double threshold );
	static void s_init_gap_frequency_threshold( double threshold );
	static void s_init_allele_state_rule( const std::string& rule_string );
	static void s_init_no_filter_alignment( bool flag );
	static void s_init_output_state_frequencies( bool flag );
	static void s_init_sample_alignment( int val );
	static void s_init_no_sample_reweighting( bool flag );
	static void s_init_rescale_sample_weights( bool flag );
	static void s_init_sample_reweighting_threshold( double threshold );
	static void s_init_output_sample_weights( bool flag );
	static void s_init_output_sample_distance_matrix( bool flag );
	static void s_init_genome_size( std::size_t npositions );
	static void s_init_linear_genome( bool flag );
	static void s_init_variable_penalty( bool flag );
	static void s_init_distance_penalty_threshold( std::size_t value );
	static void s_init_distance_penalty_shape( double value );
	static void s_init_distance_penalty_scaling( double value );
	static void s_init_output_indexing_base( std::size_t base_index );
	static void s_init_input_indexing_base( std::size_t base_index );
	static void s_init_begin_locus( int locus );
	static void s_init_end_locus( int locus );

	po::options_description
#ifdef APEGRUNT_STANDALONE_BUILD
		// Apegrunt standalone binary options
		m_standalone_options /*("apegrunt executable options")*/,
#endif // APEGRUNT_STANDALONE_BUILD

		// Apegrunt library options
		m_alignment_options /*("apegrunt alignment preprocessing options")*/,
		m_algorithm_options /*("apegrunt algorithm options")*/
	;
};

class Apegrunt_options : public Apegrunt_options_base
{
public:
	using Apegrunt_options_base::Apegrunt_options_base; // pull in base class constructor(s)
};

} // namespace apegrunt

#endif // APEGRUNT_OPTIONS_H
