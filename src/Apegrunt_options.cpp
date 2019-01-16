/** @file Apegrunt_options.cpp

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

#include "apegrunt/Apegrunt_options.h"
#include "apegrunt/Apegrunt_version.h"

namespace apegrunt {

std::ostream* Apegrunt_options_base::s_out = nullptr;
std::ostream* Apegrunt_options_base::s_err = nullptr;

uint Apegrunt_options_base::s_state = 1;
bool Apegrunt_options_base::s_verbose = false;

std::size_t Apegrunt_options_base::s_output_indexing_base = 1;
std::size_t Apegrunt_options_base::s_input_indexing_base = 1;

int Apegrunt_options_base::s_begin_locus=Apegrunt_options_base::s_input_indexing_base;
int Apegrunt_options_base::s_end_locus=-1;

//std::string Apegrunt_options_base::s_outfile_name = "apegrunt.log";
//std::string Apegrunt_options_base::s_errfile_name = "apegrunt.err";

#ifndef APEGRUNT_NO_TBB // Threading with Threading Building Blocks
int Apegrunt_options_base::s_threads = -1;
#else
int Apegrunt_options_base::s_threads = 1;
#endif // APEGRUNT_NO_TBB

#ifndef APEGRUNT_NO_MPI
int Apegrunt_options_base::s_nodes = -1;
#else
int Apegrunt_options_base::s_nodes = 1;
#endif // APEGRUNT_NO_MPI

#ifndef APEGRUNT_NO_CUDA
bool Apegrunt_options_base::s_use_cuda = true;
#else
bool Apegrunt_options_base::s_use_cuda = false;
#endif

#ifdef APEGRUNT_STANDALONE_BUILD
std::string Apegrunt_options_base::s_options_string;

const std::string Apegrunt_options_base::s_usage_string(
	  std::string("Usage: apegrunt") /*+ std::string(argv[0])*/ + " [options] <alignmentfile> [-o <outputfile>]\nOption '--help' will print a list of available options.\n"
);
#endif // APEGRUNT_STANDALONE_BUILD

const std::string Apegrunt_options_base::s_title_string(
	  std::string("Apegrunt: A library for parsing, processing and storing alignments of categorical variables.\n")
);

#define STRINGIFY(x) #x
#define TOSTRING(x) STRINGIFY(x)

const std::string Apegrunt_options_base::s_version_string(
	std::string("Apegrunt version ") + std::to_string(Apegrunt_version::s_MajorVersion) + "." + std::to_string(Apegrunt_version::s_MinorVersion) + "." + std::to_string(Apegrunt_version::s_SubminorVersion)
	+ " | revision " + TOSTRING(APEGRUNT_GIT_BRANCH) + "-" + TOSTRING(APEGRUNT_GIT_COMMIT_HASH) + " | " +
#ifdef __AVX2__
	"AVX2"
#elif __AVX__
	"AVX"
#elif __SSE2__
	"SSE2"
#else
	"generic"
#endif
	+ " build " + std::string(__DATE__) + " " + std::string(__TIME__)
);

const std::string Apegrunt_options_base::s_copyright_notice(
	std::string("Copyright (c) 2016-2018 Santeri Puranen\nLicensed under the GNU Affero General Public License version 3.\n\n")
	+ "THIS SOFTWARE IS PROVIDED \"AS IS\", WITHOUT WARRANTY OF ANY KIND."
);

const std::string Apegrunt_options_base::s_long_copyright_notice(
	std::string("Copyright (c) 2016-2017 Santeri Puranen\nLicensed under the GNU Affero General Public License version 3.\n\n")
	+ "This program is free software: you can redistribute it and/or modify\n"
	+ "it under the terms of the GNU Affero General Public License version 3 as\n"
	+ "published by the Free Software Foundation.\n\n"
	+ "This program is distributed in the hope that it will be useful,\n"
	+ "but WITHOUT ANY WARRANTY; without even the implied warranty of\n"
	+ "MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the\n"
	+ "GNU Affero General Public License for more details.\n\n"
	+ "You should have received a copy of the GNU Affero General Public License\n"
	+ "along with this program. If not, see <http://www.gnu.org/licenses/>."
);

// Alignment filtering
// /*
std::vector< std::string > Apegrunt_options_base::s_alignment_file_names;
std::string Apegrunt_options_base::s_includelist_file_name;
std::string Apegrunt_options_base::s_excludelist_file_name;
std::string Apegrunt_options_base::s_mappinglist_file_name;
std::string Apegrunt_options_base::s_samplelist_file_name;
std::string Apegrunt_options_base::s_sample_weights_file_name;
// */
double Apegrunt_options_base::s_minor_allele_frequency_threshold = 0.01; // default=0.01, i.e. the alleles below 1% are ignored
double Apegrunt_options_base::s_gap_frequency_threshold = 0.15; // default=0.15, i.e. positions with more that 15% gaps are ignored
Threshold_rule<int> Apegrunt_options_base::s_allele_state_rule(">1");

bool Apegrunt_options_base::s_no_filter_alignment = false;

bool Apegrunt_options_base::s_output_state_frequencies = false;

double Apegrunt_options_base::s_sample_reweighting_threshold = 0.9; // default=1.0, i.e. all elements should be identical for two sequences to be considered as one
bool Apegrunt_options_base::s_output_sample_weights = false;
bool Apegrunt_options_base::s_no_sample_reweighting = false;
bool Apegrunt_options_base::s_rescale_sample_weights = false;
bool Apegrunt_options_base::s_output_sample_distance_matrix = false;

std::size_t Apegrunt_options_base::s_genome_size = 0;
bool Apegrunt_options_base::s_linear_genome = false;
bool Apegrunt_options_base::s_variable_penalty = false;
std::size_t Apegrunt_options_base::s_distance_penalty_threshold = 500;
double Apegrunt_options_base::s_distance_penalty_shape = 5.0;
double Apegrunt_options_base::s_distance_penalty_scaling = 1.0;

bool Apegrunt_options_base::s_no_optimize_column_order = true; // permanently disabled

int Apegrunt_options_base::s_sample_alignment = 0; // if <1, then don't sample, else choose the indicated number of samples

Apegrunt_options_base::Apegrunt_options_base() { this->m_init(); }
Apegrunt_options_base::~Apegrunt_options_base() { }

Apegrunt_options_base::Apegrunt_options_base( std::ostream *out, std::ostream *err )
{
	Apegrunt_options_base::s_out = out;
	Apegrunt_options_base::s_err = err;
	this->m_init();
}

const std::string& Apegrunt_options_base::s_get_copyright_notice_string() { return s_copyright_notice; }
#ifdef APEGRUNT_STANDALONE_BUILD
const std::string& Apegrunt_options_base::s_get_usage_string() { return s_usage_string; }
#endif // APEGRUNT_STANDALONE_BUILD
const std::string& Apegrunt_options_base::s_get_version_string() { return s_version_string; }
const std::string& Apegrunt_options_base::s_get_title_string() { return s_title_string; }

uint Apegrunt_options_base::state() const { return s_state; }

void Apegrunt_options_base::set_out_stream( std::ostream *out ) { s_out = out->good() ? out : nullptr; }
void Apegrunt_options_base::set_err_stream( std::ostream *err ) { s_err = err->good() ? err : nullptr; }

std::ostream* Apegrunt_options_base::get_out_stream() { return s_out; }
std::ostream* Apegrunt_options_base::get_err_stream() { return s_err; }

bool Apegrunt_options_base::verbose() { return ( s_verbose && s_out ); } // be verbose only if we have a valid ostream*.
void Apegrunt_options_base::set_verbose( bool verbose ) { s_verbose = verbose; }

// alignment processing
bool Apegrunt_options_base::has_alignment_filenames() { return !s_alignment_file_names.empty(); }
const std::vector< std::string >& Apegrunt_options_base::get_alignment_filenames() { return s_alignment_file_names; }

bool Apegrunt_options_base::has_includelist_filename() { return !s_includelist_file_name.empty(); }
const std::string& Apegrunt_options_base::get_includelist_filename() { return s_includelist_file_name; }

bool Apegrunt_options_base::has_excludelist_filename() { return !s_excludelist_file_name.empty(); }
const std::string& Apegrunt_options_base::get_excludelist_filename() { return s_excludelist_file_name; }

bool Apegrunt_options_base::has_mappinglist_filename() { return !s_mappinglist_file_name.empty(); }
const std::string& Apegrunt_options_base::get_mappinglist_filename() { return s_mappinglist_file_name; }

bool Apegrunt_options_base::has_samplelist_filename() { return !s_samplelist_file_name.empty(); }
const std::string& Apegrunt_options_base::get_samplelist_filename() { return s_samplelist_file_name; }

bool Apegrunt_options_base::has_sample_weights_filename() { return !s_sample_weights_file_name.empty(); }
const std::string& Apegrunt_options_base::get_sample_weights_filename() { return s_sample_weights_file_name; }

//bool Apegrunt_options_base::fuse_duplicates() { return s_fuse_duplicates; }
double Apegrunt_options_base::get_minor_allele_frequency_threshold() { return s_minor_allele_frequency_threshold; }
double Apegrunt_options_base::get_gap_frequency_threshold() { return s_gap_frequency_threshold; }
Threshold_rule<int> Apegrunt_options_base::get_allele_state_rule() { return s_allele_state_rule; }
bool Apegrunt_options_base::filter_alignment() { return !s_no_filter_alignment; }
bool Apegrunt_options_base::output_state_frequencies() { return s_output_state_frequencies; }
bool Apegrunt_options_base::optimize_column_order() { return !s_no_optimize_column_order; }
int Apegrunt_options_base::sample_alignment() { return s_sample_alignment; }

double Apegrunt_options_base::sample_reweighting_threshold() { return s_sample_reweighting_threshold; }
bool Apegrunt_options_base::reweight_samples() { return !s_no_sample_reweighting; }
bool Apegrunt_options_base::rescale_sample_weights() { return s_rescale_sample_weights; }
bool Apegrunt_options_base::output_sample_weights() { return s_output_sample_weights; }
bool Apegrunt_options_base::output_sample_distance_matrix() { return s_output_sample_distance_matrix; }

std::size_t Apegrunt_options_base::genome_size() { return s_genome_size; }
bool Apegrunt_options_base::linear_genome() { return s_linear_genome; }
bool Apegrunt_options_base::variable_penalty() { return s_variable_penalty; }
std::size_t Apegrunt_options_base::distance_penalty_threshold() { return s_distance_penalty_threshold; }
double Apegrunt_options_base::distance_penalty_shape() { return s_distance_penalty_shape; }
double Apegrunt_options_base::distance_penalty_scaling() { return s_distance_penalty_scaling; }

std::size_t Apegrunt_options_base::get_output_indexing_base() { return s_output_indexing_base; }
void Apegrunt_options_base::set_output_indexing_base( std::size_t base_index ) { s_output_indexing_base = base_index; }
std::size_t Apegrunt_options_base::get_input_indexing_base() { return s_input_indexing_base; }
void Apegrunt_options_base::set_input_indexing_base( std::size_t base_index ) { s_input_indexing_base = base_index; }
int Apegrunt_options_base::get_begin_locus() { return s_begin_locus-s_input_indexing_base; }
int Apegrunt_options_base::get_end_locus() { return s_end_locus-s_input_indexing_base; }

bool Apegrunt_options_base::cuda() { return s_use_cuda; }
#ifndef APEGRUNT_NO_CUDA
void Apegrunt_options_base::set_cuda( bool use_cuda ) { s_use_cuda = use_cuda; }
#else
void Apegrunt_options_base::set_cuda( bool use_cuda ) { } // do nothing
#endif // APEGRUNT_NO_CUDA

int Apegrunt_options_base::threads() { return s_threads; }
#ifndef APEGRUNT_NO_TBB // Threading with Threading Building Blocks
void Apegrunt_options_base::set_threads( int nthreads ) { s_threads = nthreads; }
#else
void Apegrunt_options_base::set_threads( int nthreads ) { } // do nothing
#endif // APEGRUNT_NO_TBB

void Apegrunt_options_base::m_init()
{
	namespace po = boost::program_options;

#ifdef APEGRUNT_STANDALONE_BUILD
m_general_options.add_options()
		("outfile", po::value< std::string >( &Apegrunt_options_base::s_outfile_name )->default_value(Apegrunt_options_base::s_outfile_name), "Log filename.")
		("errfile", po::value< std::string >( &Apegrunt_options_base::s_errfile_name )->default_value(Apegrunt_options_base::s_errfile_name), "Error log filename.")
	;
	m_parallel_options.add_options()
#ifndef APEGRUNT_NO_TBB // Threading with Threading Building Blocks
		("threads,t", po::value< int >( &Apegrunt_options_base::s_threads )->default_value(Apegrunt_options_base::s_threads)->notifier(Apegrunt_options_base::s_init_threads), "Number of threads per shared memory node (-1=use all hardware threads that the OS/environment exposes).")
#endif // APEGRUNT_NO_TBB
#ifndef APEGRUNT_NO_CUDA
		("cuda", po::bool_switch( &Apegrunt_options_base::s_use_cuda )->default_value(Apegrunt_options_base::s_use_cuda)->notifier(Apegrunt_options_base::s_init_use_cuda), "Use CUDA devices, if available.")
#endif // APEGRUNT_NO_CUDA
	;
#endif // APEGRUNT_STANDALONE_BUILD

	m_alignment_options.add_options()
		("alignmentfile", po::value< std::vector< std::string > >( &Apegrunt_options_base::s_alignment_file_names )->composing(), "The input alignment filename(s). When two filenames are specified, only inter-alignment links will be probed for.")
		("include-list", po::value< std::string >( &Apegrunt_options_base::s_includelist_file_name ), "Name of file containing list of loci to include in analysis.")
		("exclude-list", po::value< std::string >( &Apegrunt_options_base::s_excludelist_file_name ), "Name of file containing list of loci to exclude from analysis.")
		("mapping-list", po::value< std::string >( &Apegrunt_options_base::s_mappinglist_file_name ), "Name of file containing loci mappings.")
		("sample-list", po::value< std::string >( &Apegrunt_options_base::s_samplelist_file_name ), "The sample filter list input filename.")

		("no-filter-alignment", po::bool_switch( &Apegrunt_options_base::s_no_filter_alignment )->default_value(Apegrunt_options_base::s_no_filter_alignment)->notifier(Apegrunt_options_base::s_init_no_filter_alignment), "Do not reduce the number of apegrunt input loci by applying MAF and GAP thresholds.")
		("maf-threshold", po::value< double >( &Apegrunt_options_base::s_minor_allele_frequency_threshold )->default_value(Apegrunt_options_base::s_minor_allele_frequency_threshold)->notifier(Apegrunt_options_base::s_init_minor_allele_frequency_threshold), "Minor state frequency threshold. Loci with less than 2 states above threshold are removed from alignment.")
		("gap-threshold", po::value< double >( &Apegrunt_options_base::s_gap_frequency_threshold )->default_value(Apegrunt_options_base::s_gap_frequency_threshold)->notifier(Apegrunt_options_base::s_init_gap_frequency_threshold), "Gap frequency threshold. Positions with a gap frequency above the threshold are excluded from the pair-analysis.")

//		("allele-state-rule", po::value< std::string >()/*->default_value( Apegrunt_options_base::s_allele_state_rule.str() )*/->notifier(Apegrunt_options_base::s_init_allele_state_rule), "Allele state filtering rule.")
		("output-state-frequencies", po::bool_switch( &Apegrunt_options_base::s_output_state_frequencies )->default_value(Apegrunt_options_base::s_output_state_frequencies)->notifier(Apegrunt_options_base::s_init_output_state_frequencies), "Write column-wise state frequencies to file.")

		("sample-alignment", po::value< int >( &Apegrunt_options_base::s_sample_alignment )->default_value(Apegrunt_options_base::s_sample_alignment)->notifier(Apegrunt_options_base::s_init_sample_alignment), "If >0 choose a random subset of samples. The number indicates sample count.")

		("no-sample-reweighting", po::bool_switch( &Apegrunt_options_base::s_no_sample_reweighting )->default_value(Apegrunt_options_base::s_no_sample_reweighting)->notifier(Apegrunt_options_base::s_init_no_sample_reweighting), "Turn sample reweighting off, i.e. do not try to correct for population structure.")
		("sample-reweighting-threshold", po::value< double >( &Apegrunt_options_base::s_sample_reweighting_threshold )->default_value(Apegrunt_options_base::s_sample_reweighting_threshold)->notifier(Apegrunt_options_base::s_init_sample_reweighting_threshold), "Fraction of identical positions required for two samples to be considered identical.")
		("rescale-sample-weights", po::bool_switch( &Apegrunt_options_base::s_rescale_sample_weights )->default_value(Apegrunt_options_base::s_rescale_sample_weights)->notifier(Apegrunt_options_base::s_init_rescale_sample_weights), "Rescale sample weights so that n(effective) == n.")
		("sample-weights",  po::value< std::string >( &Apegrunt_options_base::s_sample_weights_file_name ), "Name of file containing sample weights.")
		("output-sample-weights", po::bool_switch( &Apegrunt_options_base::s_output_sample_weights )->default_value(Apegrunt_options_base::s_output_sample_weights)->notifier(Apegrunt_options_base::s_init_output_sample_weights), "Output sample weights.")
		("output-sample-distance-matrix", po::bool_switch( &Apegrunt_options_base::s_output_sample_distance_matrix )->default_value(Apegrunt_options_base::s_output_sample_distance_matrix)->notifier(Apegrunt_options_base::s_init_output_sample_distance_matrix), "Output triangular sample-sample Hamming distance matrix.")

		("genome-size", po::value< std::size_t >( &Apegrunt_options_base::s_genome_size )->notifier(Apegrunt_options_base::s_init_genome_size), "Genome size, if different from input. Default = 0: detect size from input.")
		("linear-genome", po::bool_switch( &Apegrunt_options_base::s_linear_genome )->default_value(Apegrunt_options_base::s_linear_genome)->notifier(Apegrunt_options_base::s_init_linear_genome), "Assume linear genome in distance calculations.")
//		("variable-penalty", po::bool_switch( &Apegrunt_options_base::s_variable_penalty )->default_value(Apegrunt_options_base::s_variable_penalty)->notifier(Apegrunt_options_base::s_init_variable_penalty), "Use a column-wise varaible penalty function.")
//		("distance-penalty-threshold", po::value< std::size_t >( &Apegrunt_options_base::s_distance_penalty_threshold )->default_value(Apegrunt_options_base::s_distance_penalty_threshold)->notifier(Apegrunt_options_base::s_init_distance_penalty_threshold), "Distance threshold for the variable regularization/penalty function.")
//		("distance-penalty-shape", po::value< double >( &Apegrunt_options_base::s_distance_penalty_shape )->default_value(Apegrunt_options_base::s_distance_penalty_shape)->notifier(Apegrunt_options_base::s_init_distance_penalty_shape), "Distance penalty function shape parameter.")
//		("distance-penalty-scaling", po::value< double >( &Apegrunt_options_base::s_distance_penalty_scaling )->default_value(Apegrunt_options_base::s_distance_penalty_scaling)->notifier(Apegrunt_options_base::s_init_distance_penalty_scaling), "The penalty scaling factor for the adaptive penalty function. The maximum penalty value will be lambda_J times the scaling_factor (0.0 == flat response).")

		("input-indexing-base", po::value< std::size_t >( &Apegrunt_options_base::s_input_indexing_base )->default_value(Apegrunt_options_base::s_input_indexing_base)->notifier(Apegrunt_options_base::s_init_input_indexing_base), "Base index for input." )
		("output-indexing-base", po::value< std::size_t >( &Apegrunt_options_base::s_output_indexing_base )->default_value(Apegrunt_options_base::s_output_indexing_base)->notifier(Apegrunt_options_base::s_init_output_indexing_base), "Base index for output." )
//		("no-optimize-column-order", po::bool_switch( &Apegrunt_options_base::s_no_optimize_column_order )->default_value(Apegrunt_options_base::s_no_optimize_column_order), "Do not optimize data column order.")
	;
	m_algorithm_options.add_options()
		("begin", po::value< int >( &Apegrunt_options_base::s_begin_locus )->default_value(Apegrunt_options_base::s_begin_locus)->notifier(Apegrunt_options_base::s_init_begin_locus), "The first locus index to work on. Used to define a range.")
		("end", po::value< int >( &Apegrunt_options_base::s_end_locus )->default_value(Apegrunt_options_base::s_end_locus)->notifier(Apegrunt_options_base::s_init_end_locus), "The last locus index to work on (-1=end of input). Used to define a range.")
	;
}

void Apegrunt_options_base::AddOptions( boost::program_options::options_description *opdesc )
{
	namespace po = boost::program_options;
#ifdef APEGRUNT_STANDALONE_BUILD
	opdesc->add(m_general_options);
	opdesc->add(m_parallel_options);
#endif // APEGRUNT_STANDALONE_BUILD
	opdesc->add(m_alignment_options);
	opdesc->add(m_algorithm_options);
}

/// Check options stored in varmap. Return false if a fatal inconsistency is detected.
bool Apegrunt_options_base::CheckOptions( boost::program_options::variables_map *varmap )
{
	try
	{
		/*
		if( s_verbose && s_out )
		{
			*s_out << "apegrunt: being verbose." << std::endl;
		}
		*/

		if( !varmap->count("alignmentfile") && s_err )
		{
			*s_err << "apegrunt error: No alignment file specified!" << std::endl;
#ifdef APEGRUNT_STANDALONE_BUILD
			if( s_out )
			{
				*s_out << s_usage_string << std::endl;
			}
#endif // APEGRUNT_STANDALONE_BUILD
			return false;
		}

		if( varmap->count("include-list") && varmap->count("exclude-list") && s_err )
		{
			*s_err << "apegrunt error: options \"--include-list\" and \"--exclude-list\" are mutually exclusive! Please use only one of them." << std::endl;
			return false;
		}

	}
	catch( std::exception& e)
	{
		if( s_err )
		{
			// probably an unknown option
			*s_err << "apegrunt error: " << e.what() << "\n\n";
		}
		return false;
	}
	catch(...)
	{
		if( s_err )
		{
			*s_err << "apegrunt error: Exception of unknown type!\n\n";
		}
		return false;
	}

	return true;
}

void Apegrunt_options_base::s_init_verbose( bool flag )
{
	if( flag && s_out )
	{
		//*s_out << "apegrunt: begin verbose.\n";
	}
}

void Apegrunt_options_base::s_init_no_filter_alignment( bool flag )
{
	if( flag && s_verbose && s_out )
	{
		*s_out << "apegrunt: do not reduce the number of apegrunt analyzable input loci by applying MAF and GAP thresholds.\n";
	}
}

void Apegrunt_options_base::s_init_minor_allele_frequency_threshold( double threshold )
{
	if( s_verbose && s_out )
	{
		*s_out << "apegrunt: minor allele frequency threshold set to " << threshold << ".\n";
	}
}

void Apegrunt_options_base::s_init_gap_frequency_threshold( double threshold )
{
	if( s_verbose && s_out )
	{
		*s_out << "apegrunt: gap frequency threshold set to " << threshold << ".\n";
	}
}

void Apegrunt_options_base::s_init_output_state_frequencies( bool flag )
{
	if( flag && s_verbose && s_out )
	{
		*s_out << "apegrunt: output column-wise state frequencies to file.\n";
	}
}

void Apegrunt_options_base::s_init_no_sample_reweighting( bool flag )
{
	if( flag && s_verbose && s_out )
	{
		*s_out << "apegrunt: do not reweight samples (i.e. do not try to correct for population structure).\n";
	}
}

void Apegrunt_options_base::s_init_rescale_sample_weights( bool flag )
{
	if( flag && s_verbose && s_out )
	{
		*s_out << "apegrunt: scale sample weights so that n(effective) == n.\n";
	}
}

void Apegrunt_options_base::s_init_sample_reweighting_threshold( double threshold )
{
	if( s_verbose && s_out )
	{
		*s_out << "apegrunt: reweighting threshold set to " << threshold << ".\n";
	}
}

void Apegrunt_options_base::s_init_output_sample_weights( bool flag )
{
	if( flag && s_verbose && s_out )
	{
		*s_out << "apegrunt: output sample weights to file.\n";
	}
}

void Apegrunt_options_base::s_init_output_sample_distance_matrix( bool flag )
{
	if( flag && s_verbose && s_out )
	{
		*s_out << "apegrunt: output triangular sample-sample Hamming distance matrix to file.\n";
	}
}

void Apegrunt_options_base::s_init_genome_size( std::size_t npositions )
{
	if( npositions != 0 && s_verbose && s_out )
	{
		*s_out << "apegrunt: genome size set to " << npositions << ".\n";
	}
}

void Apegrunt_options_base::s_init_linear_genome( bool flag )
{
	if( flag && s_verbose && s_out )
	{
		*s_out << "apegrunt: assume linear genome for distance calculations.\n";
	}
}

void Apegrunt_options_base::s_init_variable_penalty( bool flag )
{
	if( flag && s_verbose && s_out )
	{
		*s_out << "apegrunt: use column-wise varible penalty function.\n";
	}
}

void Apegrunt_options_base::s_init_distance_penalty_threshold( std::size_t val )
{
	if( s_verbose && s_out )
	{
		*s_out << "apegrunt: distance penalty threshold set to " << val << ".\n";
	}
}

void Apegrunt_options_base::s_init_distance_penalty_shape( double val )
{
	if( s_verbose && s_out )
	{
		*s_out << "apegrunt: distance penalty function shape parameter set to " << val << ".\n";
	}
}

void Apegrunt_options_base::s_init_distance_penalty_scaling( double val )
{
	if( s_verbose && s_out )
	{
		*s_out << "apegrunt: distance penalty function scaling parameter set to " << val << ".\n";
	}
}

#ifndef APEGRUNT_NO_TBB // Threading with Threading Building Blocks
void Apegrunt_options_base::s_init_threads( int nthreads )
{
	if( s_verbose && s_out )
	{
		*s_out << "apegrunt: user requests " << nthreads << " compute threads.\n";
	}
}
#endif // APEGRUNT_NO_TBB

#ifndef APEGRUNT_NO_CUDA
void Apegrunt_options_base::s_init_use_cuda( bool use_cuda )
{
	if( s_verbose && s_out && s_use_cuda )
	{
		*s_out << "apegrunt: use CUDA if available.\n";
	}
}
#endif // APEGRUNT_NO_CUDA

void Apegrunt_options_base::s_init_fuse_duplicates( bool flag )
{
	if( s_verbose && s_out && flag )
	{
		*s_out << "apegrunt: fuse duplicate sequences.\n";
	}
}

void Apegrunt_options_base::s_init_sample_alignment( int val )
{
	if( s_verbose && s_out && val )
	{
		if( val > 0 )
		{
			*s_out << "apegrunt: choose a subset of " << val << " random samples from the alignment.\n";
		}
	}
}

void Apegrunt_options_base::s_init_output_indexing_base( std::size_t base_index )
{
	if( s_verbose && s_out )
	{
		*s_out << "apegrunt: output indexing will begin at " << base_index << ".\n";
	}
}

void Apegrunt_options_base::s_init_input_indexing_base( std::size_t base_index )
{
	if( s_verbose && s_out )
	{
		*s_out << "apegrunt: will assume that input indexing begins at " << base_index << ".\n";
	}
}

void Apegrunt_options_base::s_init_allele_state_rule( const std::string& rule_string )
{
	s_allele_state_rule.set_rule(rule_string);
	if( s_verbose && s_out )
	{
		*s_out << "apegrunt: filter for loci that have " << s_allele_state_rule << " states.\n";
	}
}

void Apegrunt_options_base::s_init_begin_locus( int locus ) { }
void Apegrunt_options_base::s_init_end_locus( int locus ) { }

} // namespace apegrunt

