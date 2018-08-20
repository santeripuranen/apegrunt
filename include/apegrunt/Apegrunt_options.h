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
	static bool CheckOptions( po::variables_map *varmap );
	static const std::string& s_get_copyright_notice_string();
#ifdef APEGRUNT_STANDALONE_BUILD
	static const std::string& s_get_usage_string();
#endif // APEGRUNT_STANDALONE_BUILD
	static const std::string& s_get_version_string();
	static const std::string& s_get_title_string();

	uint state() const;

	//> Set an ostream. An invalid ostream* (as in "out->good() == false", will reset internal ostream ("ostream* == null_ptr").
	static void set_out_stream( std::ostream* out );
	static void set_err_stream( std::ostream* err );
	static std::ostream* get_out_stream();
	static std::ostream* get_err_stream();

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

	static double get_minor_allele_frequency_threshold();
	static double get_gap_frequency_threshold();
	static Threshold_rule<int> get_allele_state_rule();

	// alignment processing
	static bool fuse_duplicates();
	static bool filter_alignment();
	static bool output_state_frequencies();
	static bool optimize_column_order();
	static int sample_alignment();
	static const std::vector<double>& sample_fractions();
	static double sample_reweighting_threshold();
	static bool output_sample_weights();
	static bool reweight_samples();

	static bool has_sample_weights_filename();
	static const std::string& get_sample_weights_filename();

	static bool rescale_sample_weights();
	static bool output_sample_distance_matrix();

	static bool linear_genome();
	static bool variable_penalty();
	static std::size_t distance_penalty_threshold();
	static double distance_penalty_shape();
	static double distance_penalty_scaling();

	//> Test if textual output is desired. If true, then a call to get_out_stream() is guaranteed to return a valid (as in != null_ptr) ostream*.
	static bool verbose();
	static void set_verbose( bool verbose=true );

	static std::size_t get_output_indexing_base();
	static void set_output_indexing_base( std::size_t base_index );
	static std::size_t get_input_indexing_base();
	static void set_input_indexing_base( std::size_t base_index );
	static int get_begin_locus();
	static int get_end_locus();

	static bool cuda();
	static void set_cuda( bool use_cuda=true );

	static int threads();
	static void set_threads( int nthreads );

private:
// all the member variables

	static uint s_state; // 1 for normal operation, 0 signals a wish to (gracefully) terminate process
	static bool s_verbose;

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

	static const std::string s_title_string;
#ifdef APEGRUNT_STANDALONE_BUILD
	static const std::string s_usage_string;
#endif // APEGRUNT_STANDALONE_BUILD
	static const std::string s_version_string;
	static const std::string s_copyright_notice;
	static const std::string s_long_copyright_notice;
	static std::string s_options_string;
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
		m_general_options /*("apegrunt general options")*/,
		m_parallel_options /*("apegrunt parallel options")*/,
#endif // APEGRUNT_STANDALONE_BUILD
		m_alignment_options /*("apegrunt alignment preprocessing options")*/,
		m_algorithm_options /*("apegrunt algorithm options")*/
	;
};

class Apegrunt_options : public Apegrunt_options_base
{
public:
	using Apegrunt_options_base::Apegrunt_options_base; // pull in base class constructor(s)
/*
	void AddOptions( po::options_description *opdesc );
	static bool CheckOptions( po::variables_map *varmap );

	static const std::string& s_get_copyright_notice_string();
#ifdef APEGRUNT_STANDALONE_BUILD
	static const std::string& s_get_usage_string();
#endif // APEGRUNT_STANDALONE_BUILD
	static const std::string& s_get_version_string();
	static const std::string& s_get_title_string();

	uint state() const;

	//> Set an ostream. An invalid ostream* (as in "out->good() == false", will reset internal ostream ("ostream* == null_ptr").
	static void set_out_stream( std::ostream* out );
	static void set_err_stream( std::ostream* err );
	static std::ostream* get_out_stream();
	static std::ostream* get_err_stream();

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

	static double get_minor_allele_frequency_threshold();
	static double get_gap_frequency_threshold();
	static Threshold_rule<int> get_allele_state_rule();

	// alignment processing
	static bool fuse_duplicates();
	static bool filter_alignment();
	static bool output_state_frequencies();
	static bool optimize_column_order();
	static int sample_alignment();
	static const std::vector<double>& sample_fractions();
	static double sample_reweighting_threshold();
	static bool output_sample_weights();
	static bool reweight_samples();

	static bool has_sample_weights_filename();
	static const std::string& get_sample_weights_filename();

	static bool rescale_sample_weights();
	static bool output_sample_distance_matrix();

	static bool linear_genome();
	static bool variable_penalty();
	static std::size_t distance_penalty_threshold();
	static double distance_penalty_shape();
	static double distance_penalty_scaling();

	//> Test if textual output is desired. If true, then a call to get_out_stream() is guaranteed to return a valid (as in != null_ptr) ostream*.
	static bool verbose();
	static void set_verbose( bool verbose=true );

	static std::size_t get_output_indexing_base();
	static void set_output_indexing_base( std::size_t base_index );
	static std::size_t get_input_indexing_base();
	static void set_input_indexing_base( std::size_t base_index );
	static int get_begin_locus();
	static int get_end_locus();

	static bool cuda();
	static void set_cuda( bool use_cuda=true );

	static int threads();
	static void set_threads( int nthreads );

private:
	static uint s_state; // 1 for normal operation, 0 signals a wish to (gracefully) terminate process
	static bool s_verbose;

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

	static const std::string s_title_string;
#ifdef APEGRUNT_STANDALONE_BUILD
	static const std::string s_usage_string;
#endif // APEGRUNT_STANDALONE_BUILD
	static const std::string s_version_string;
	static const std::string s_copyright_notice;
	static const std::string s_long_copyright_notice;
	static std::string s_options_string;
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
		m_general_options, // ("apegrunt general options"),
		m_parallel_options // ("apegrunt parallel options"),
#endif // APEGRUNT_STANDALONE_BUILD
		m_alignment_options, // ("apegrunt alignment preprocessing options"),
		m_algorithm_options // ("apegrunt algorithm options")
	;
*/
};

} // namespace apegrunt

#endif // APEGRUNT_OPTIONS_H
