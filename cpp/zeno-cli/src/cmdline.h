/** @file cmdline.h
 *  @brief The header file for the command line option parser
 *  generated by GNU Gengetopt version 2.23
 *  http://www.gnu.org/software/gengetopt.
 *  DO NOT modify this file, since it can be overwritten
 *  @author GNU Gengetopt */

#ifndef CMDLINE_H
#define CMDLINE_H

/* If we use autoconf.  */
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h> /* for FILE */

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#ifndef CMDLINE_PARSER_PACKAGE
/** @brief the program name (used for printing errors) */
#define CMDLINE_PARSER_PACKAGE "zeno"
#endif

#ifndef CMDLINE_PARSER_PACKAGE_NAME
/** @brief the complete program name (used for help and version) */
#define CMDLINE_PARSER_PACKAGE_NAME "zeno"
#endif

#ifndef CMDLINE_PARSER_VERSION
/** @brief the program version */
#define CMDLINE_PARSER_VERSION "5.3.1"
#endif

/** @brief Where the command line options are stored */
struct gengetopt_args_info
{
  const char *help_help; /**< @brief Print help and exit help description.  */
  const char *version_help; /**< @brief Print version and exit help description.  */
  char * config_file_arg;	/**< @brief Config file name.  */
  char * config_file_orig;	/**< @brief Config file name original value given at command line.  */
  const char *config_file_help; /**< @brief Config file name help description.  */
  char * input_file_arg;	/**< @brief Input file name.  */
  char * input_file_orig;	/**< @brief Input file name original value given at command line.  */
  const char *input_file_help; /**< @brief Input file name help description.  */
  char * csv_output_file_arg;	/**< @brief Write output in CSV format to the specified file in addition to displaying the regular output.  */
  char * csv_output_file_orig;	/**< @brief Write output in CSV format to the specified file in addition to displaying the regular output original value given at command line.  */
  const char *csv_output_file_help; /**< @brief Write output in CSV format to the specified file in addition to displaying the regular output help description.  */
  #if defined(HAVE_LONG_LONG) || defined(HAVE_LONG_LONG_INT)
  long long int num_walks_arg;	/**< @brief Number of walk-on-spheres walks to perform.  Fewer walks may be performed if another stopping condition is reached first.  */
  #else
  long num_walks_arg;	/**< @brief Number of walk-on-spheres walks to perform.  Fewer walks may be performed if another stopping condition is reached first.  */
  #endif
  char * num_walks_orig;	/**< @brief Number of walk-on-spheres walks to perform.  Fewer walks may be performed if another stopping condition is reached first original value given at command line.  */
  const char *num_walks_help; /**< @brief Number of walk-on-spheres walks to perform.  Fewer walks may be performed if another stopping condition is reached first help description.  */
  #if defined(HAVE_LONG_LONG) || defined(HAVE_LONG_LONG_INT)
  long long int num_interior_samples_arg;	/**< @brief Number of interior samples to take.  Fewer samples may be taken if another stopping condition is reached first.  */
  #else
  long num_interior_samples_arg;	/**< @brief Number of interior samples to take.  Fewer samples may be taken if another stopping condition is reached first.  */
  #endif
  char * num_interior_samples_orig;	/**< @brief Number of interior samples to take.  Fewer samples may be taken if another stopping condition is reached first original value given at command line.  */
  const char *num_interior_samples_help; /**< @brief Number of interior samples to take.  Fewer samples may be taken if another stopping condition is reached first help description.  */
  double max_rsd_capacitance_arg;	/**< @brief Perform walk-on-spheres walks until the relative standard deviation of the capacitance drops below this value.  This value may not be reached if another stopping condition is reached first.  Relative standard deviation is defined as (Standard_Deviation/Mean)*100%.  */
  char * max_rsd_capacitance_orig;	/**< @brief Perform walk-on-spheres walks until the relative standard deviation of the capacitance drops below this value.  This value may not be reached if another stopping condition is reached first.  Relative standard deviation is defined as (Standard_Deviation/Mean)*100% original value given at command line.  */
  const char *max_rsd_capacitance_help; /**< @brief Perform walk-on-spheres walks until the relative standard deviation of the capacitance drops below this value.  This value may not be reached if another stopping condition is reached first.  Relative standard deviation is defined as (Standard_Deviation/Mean)*100% help description.  */
  double max_rsd_polarizability_arg;	/**< @brief Perform walk-on-spheres walks until the relative standard deviation of the mean electric polarizability drops below this value. This value may not be reached if another stopping condition is reached first. Relative standard deviation is defined as (Standard_Deviation/Mean)*100%.  */
  char * max_rsd_polarizability_orig;	/**< @brief Perform walk-on-spheres walks until the relative standard deviation of the mean electric polarizability drops below this value. This value may not be reached if another stopping condition is reached first. Relative standard deviation is defined as (Standard_Deviation/Mean)*100% original value given at command line.  */
  const char *max_rsd_polarizability_help; /**< @brief Perform walk-on-spheres walks until the relative standard deviation of the mean electric polarizability drops below this value. This value may not be reached if another stopping condition is reached first. Relative standard deviation is defined as (Standard_Deviation/Mean)*100% help description.  */
  double max_rsd_volume_arg;	/**< @brief Take interior samples until the relative standard deviation of volume drops below this value.  This value may not be reached if another stopping condition is reached first.  Relative standard deviation is defined as (Standard_Deviation/Mean)*100%.  */
  char * max_rsd_volume_orig;	/**< @brief Take interior samples until the relative standard deviation of volume drops below this value.  This value may not be reached if another stopping condition is reached first.  Relative standard deviation is defined as (Standard_Deviation/Mean)*100% original value given at command line.  */
  const char *max_rsd_volume_help; /**< @brief Take interior samples until the relative standard deviation of volume drops below this value.  This value may not be reached if another stopping condition is reached first.  Relative standard deviation is defined as (Standard_Deviation/Mean)*100% help description.  */
  #if defined(HAVE_LONG_LONG) || defined(HAVE_LONG_LONG_INT)
  long long int min_num_walks_arg;	/**< @brief Minimum number of walk-on-spheres walks to perform when using max-rsd stopping conditions (default=1000).  */
  #else
  long min_num_walks_arg;	/**< @brief Minimum number of walk-on-spheres walks to perform when using max-rsd stopping conditions (default=1000).  */
  #endif
  char * min_num_walks_orig;	/**< @brief Minimum number of walk-on-spheres walks to perform when using max-rsd stopping conditions original value given at command line.  */
  const char *min_num_walks_help; /**< @brief Minimum number of walk-on-spheres walks to perform when using max-rsd stopping conditions help description.  */
  #if defined(HAVE_LONG_LONG) || defined(HAVE_LONG_LONG_INT)
  long long int min_num_interior_samples_arg;	/**< @brief Minimum number of interior samples to take when using max-rsd stopping conditions (default=10000).  */
  #else
  long min_num_interior_samples_arg;	/**< @brief Minimum number of interior samples to take when using max-rsd stopping conditions (default=10000).  */
  #endif
  char * min_num_interior_samples_orig;	/**< @brief Minimum number of interior samples to take when using max-rsd stopping conditions original value given at command line.  */
  const char *min_num_interior_samples_help; /**< @brief Minimum number of interior samples to take when using max-rsd stopping conditions help description.  */
  double max_run_time_arg;	/**< @brief Max time (in seconds) that the program is allowed to run.  If this time is reached, the computation will be stopped and the results computed so far will be displayed.  */
  char * max_run_time_orig;	/**< @brief Max time (in seconds) that the program is allowed to run.  If this time is reached, the computation will be stopped and the results computed so far will be displayed original value given at command line.  */
  const char *max_run_time_help; /**< @brief Max time (in seconds) that the program is allowed to run.  If this time is reached, the computation will be stopped and the results computed so far will be displayed help description.  */
  int num_threads_arg;	/**< @brief Number of threads to use  (default=Number of logical cores).  */
  char * num_threads_orig;	/**< @brief Number of threads to use  (default=Number of logical cores) original value given at command line.  */
  const char *num_threads_help; /**< @brief Number of threads to use  (default=Number of logical cores) help description.  */
  int seed_arg;	/**< @brief Seed for the random number generator  (default=Randomly set).  */
  char * seed_orig;	/**< @brief Seed for the random number generator  (default=Randomly set) original value given at command line.  */
  const char *seed_help; /**< @brief Seed for the random number generator  (default=Randomly set) help description.  */
  char * surface_points_file_arg;	/**< @brief Name of file for writing the surface points from Walk-on-Spheres.  */
  char * surface_points_file_orig;	/**< @brief Name of file for writing the surface points from Walk-on-Spheres original value given at command line.  */
  const char *surface_points_file_help; /**< @brief Name of file for writing the surface points from Walk-on-Spheres help description.  */
  char * interior_points_file_arg;	/**< @brief Name of file for writing the interior sample points.  */
  char * interior_points_file_orig;	/**< @brief Name of file for writing the interior sample points original value given at command line.  */
  const char *interior_points_file_help; /**< @brief Name of file for writing the interior sample points help description.  */
  const char *print_counts_help; /**< @brief Print statistics related to counts of hit points help description.  */
  const char *print_benchmarks_help; /**< @brief Print detailed RAM and timing information help description.  */
  #if defined(HAVE_LONG_LONG) || defined(HAVE_LONG_LONG_INT)
  long long int virial_steps_arg;	/**< @brief Number of Monte Carlo steps to perform to compute a virial coefficient.  */
  #else
  long virial_steps_arg;	/**< @brief Number of Monte Carlo steps to perform to compute a virial coefficient.  */
  #endif
  char * virial_steps_orig;	/**< @brief Number of Monte Carlo steps to perform to compute a virial coefficient original value given at command line.  */
  const char *virial_steps_help; /**< @brief Number of Monte Carlo steps to perform to compute a virial coefficient help description.  */
  int virial_coefficient_order_arg;	/**< @brief Order of virial coefficient to compute.  */
  char * virial_coefficient_order_orig;	/**< @brief Order of virial coefficient to compute original value given at command line.  */
  const char *virial_coefficient_order_help; /**< @brief Order of virial coefficient to compute help description.  */
  double virial_reference_diameter_arg;	/**< @brief Virial coefficient reference diameter.  */
  char * virial_reference_diameter_orig;	/**< @brief Virial coefficient reference diameter original value given at command line.  */
  const char *virial_reference_diameter_help; /**< @brief Virial coefficient reference diameter help description.  */
  double temperature_arg;	/**< @brief Temperature.  */
  char * temperature_orig;	/**< @brief Temperature original value given at command line.  */
  const char *temperature_help; /**< @brief Temperature help description.  */
  
  unsigned int help_given ;	/**< @brief Whether help was given.  */
  unsigned int version_given ;	/**< @brief Whether version was given.  */
  unsigned int config_file_given ;	/**< @brief Whether config-file was given.  */
  unsigned int input_file_given ;	/**< @brief Whether input-file was given.  */
  unsigned int csv_output_file_given ;	/**< @brief Whether csv-output-file was given.  */
  unsigned int num_walks_given ;	/**< @brief Whether num-walks was given.  */
  unsigned int num_interior_samples_given ;	/**< @brief Whether num-interior-samples was given.  */
  unsigned int max_rsd_capacitance_given ;	/**< @brief Whether max-rsd-capacitance was given.  */
  unsigned int max_rsd_polarizability_given ;	/**< @brief Whether max-rsd-polarizability was given.  */
  unsigned int max_rsd_volume_given ;	/**< @brief Whether max-rsd-volume was given.  */
  unsigned int min_num_walks_given ;	/**< @brief Whether min-num-walks was given.  */
  unsigned int min_num_interior_samples_given ;	/**< @brief Whether min-num-interior-samples was given.  */
  unsigned int max_run_time_given ;	/**< @brief Whether max-run-time was given.  */
  unsigned int num_threads_given ;	/**< @brief Whether num-threads was given.  */
  unsigned int seed_given ;	/**< @brief Whether seed was given.  */
  unsigned int surface_points_file_given ;	/**< @brief Whether surface-points-file was given.  */
  unsigned int interior_points_file_given ;	/**< @brief Whether interior-points-file was given.  */
  unsigned int print_counts_given ;	/**< @brief Whether print-counts was given.  */
  unsigned int print_benchmarks_given ;	/**< @brief Whether print-benchmarks was given.  */
  unsigned int virial_steps_given ;	/**< @brief Whether virial-steps was given.  */
  unsigned int virial_coefficient_order_given ;	/**< @brief Whether virial-coefficient-order was given.  */
  unsigned int virial_reference_diameter_given ;	/**< @brief Whether virial-reference-diameter was given.  */
  unsigned int temperature_given ;	/**< @brief Whether temperature was given.  */

} ;

/** @brief The additional parameters to pass to parser functions */
struct cmdline_parser_params
{
  int override; /**< @brief whether to override possibly already present options (default 0) */
  int initialize; /**< @brief whether to initialize the option structure gengetopt_args_info (default 1) */
  int check_required; /**< @brief whether to check that all required options were provided (default 1) */
  int check_ambiguity; /**< @brief whether to check for options already specified in the option structure gengetopt_args_info (default 0) */
  int print_errors; /**< @brief whether getopt_long should print an error message for a bad option (default 1) */
} ;

/** @brief the purpose string of the program */
extern const char *gengetopt_args_info_purpose;
/** @brief the usage string of the program */
extern const char *gengetopt_args_info_usage;
/** @brief the description string of the program */
extern const char *gengetopt_args_info_description;
/** @brief all the lines making the help output */
extern const char *gengetopt_args_info_help[];

/**
 * The command line parser
 * @param argc the number of command line options
 * @param argv the command line options
 * @param args_info the structure where option information will be stored
 * @return 0 if everything went fine, NON 0 if an error took place
 */
int cmdline_parser (int argc, char **argv,
  struct gengetopt_args_info *args_info);

/**
 * The command line parser (version with additional parameters - deprecated)
 * @param argc the number of command line options
 * @param argv the command line options
 * @param args_info the structure where option information will be stored
 * @param override whether to override possibly already present options
 * @param initialize whether to initialize the option structure my_args_info
 * @param check_required whether to check that all required options were provided
 * @return 0 if everything went fine, NON 0 if an error took place
 * @deprecated use cmdline_parser_ext() instead
 */
int cmdline_parser2 (int argc, char **argv,
  struct gengetopt_args_info *args_info,
  int override, int initialize, int check_required);

/**
 * The command line parser (version with additional parameters)
 * @param argc the number of command line options
 * @param argv the command line options
 * @param args_info the structure where option information will be stored
 * @param params additional parameters for the parser
 * @return 0 if everything went fine, NON 0 if an error took place
 */
int cmdline_parser_ext (int argc, char **argv,
  struct gengetopt_args_info *args_info,
  struct cmdline_parser_params *params);

/**
 * Save the contents of the option struct into an already open FILE stream.
 * @param outfile the stream where to dump options
 * @param args_info the option struct to dump
 * @return 0 if everything went fine, NON 0 if an error took place
 */
int cmdline_parser_dump(FILE *outfile,
  struct gengetopt_args_info *args_info);

/**
 * Save the contents of the option struct into a (text) file.
 * This file can be read by the config file parser (if generated by gengetopt)
 * @param filename the file where to save
 * @param args_info the option struct to save
 * @return 0 if everything went fine, NON 0 if an error took place
 */
int cmdline_parser_file_save(const char *filename,
  struct gengetopt_args_info *args_info);

/**
 * Print the help
 */
void cmdline_parser_print_help(void);
/**
 * Print the version
 */
void cmdline_parser_print_version(void);

/**
 * Initializes all the fields a cmdline_parser_params structure 
 * to their default values
 * @param params the structure to initialize
 */
void cmdline_parser_params_init(struct cmdline_parser_params *params);

/**
 * Allocates dynamically a cmdline_parser_params structure and initializes
 * all its fields to their default values
 * @return the created and initialized cmdline_parser_params structure
 */
struct cmdline_parser_params *cmdline_parser_params_create(void);

/**
 * Initializes the passed gengetopt_args_info structure's fields
 * (also set default values for options that have a default)
 * @param args_info the structure to initialize
 */
void cmdline_parser_init (struct gengetopt_args_info *args_info);
/**
 * Deallocates the string fields of the gengetopt_args_info structure
 * (but does not deallocate the structure itself)
 * @param args_info the structure to deallocate
 */
void cmdline_parser_free (struct gengetopt_args_info *args_info);

/**
 * The config file parser (deprecated version)
 * @param filename the name of the config file
 * @param args_info the structure where option information will be stored
 * @param override whether to override possibly already present options
 * @param initialize whether to initialize the option structure my_args_info
 * @param check_required whether to check that all required options were provided
 * @return 0 if everything went fine, NON 0 if an error took place
 * @deprecated use cmdline_parser_config_file() instead
 */
int cmdline_parser_configfile (const char *filename,
  struct gengetopt_args_info *args_info,
  int override, int initialize, int check_required);

/**
 * The config file parser
 * @param filename the name of the config file
 * @param args_info the structure where option information will be stored
 * @param params additional parameters for the parser
 * @return 0 if everything went fine, NON 0 if an error took place
 */
int cmdline_parser_config_file (const char *filename,
  struct gengetopt_args_info *args_info,
  struct cmdline_parser_params *params);

/**
 * Checks that all the required options were specified
 * @param args_info the structure to check
 * @param prog_name the name of the program that will be used to print
 *   possible errors
 * @return
 */
int cmdline_parser_required (struct gengetopt_args_info *args_info,
  const char *prog_name);


#ifdef __cplusplus
}
#endif /* __cplusplus */
#endif /* CMDLINE_H */
