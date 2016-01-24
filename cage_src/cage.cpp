/**
 * Copyright 2015 Adam Bloniarz
 * 
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 * 
 *     http://www.apache.org/licenses/LICENSE-2.0
 * 
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
#include <iostream>
#include <vector>
#include <sys/time.h>
#include <time.h>
#include <thread>
#include <math.h>
#include <tclap/CmdLine.h>
#include <utility>
#include "dataloader.h"
#include "pelt.h"

using namespace std;
using cage::CageParams;
using cage::Dataloader;
using cage::pelt;

CageParams parse_parameters(int argc, char** argv) {

  CageParams parsed_params;

  try {
    TCLAP::CmdLine cmd("CAGe - Changepoint detection for efficient variant calling", ' ', "0.1");

    TCLAP::UnlabeledValueArg<string> contig_arg("contig", "contig name", true, "", "contig", cmd);
    TCLAP::UnlabeledValueArg<int> start_arg("start", "start position", true, 0, "start", cmd);
    TCLAP::UnlabeledValueArg<int> end_arg("end", "end position", true, 0, "end", cmd);
    TCLAP::UnlabeledValueArg<int> step_arg("step", "step size", true, 200, "stepsize", cmd);
    TCLAP::UnlabeledValueArg<double> beta_arg("beta", "beta parameter for PELT", true, 0, "beta", cmd);
    TCLAP::UnlabeledValueArg<string> output_file_arg("output_file", "File to output the changepoints determined by CAGe", true, "", "cage_output_file", cmd);
    TCLAP::SwitchArg verbose_switch_arg("v", "verbose", "print verbose output of CAGe", cmd, false);
    TCLAP::ValueArg<string> output_VCF_arg("o", "output_vcf", "File to output variants called when running CAGe", true, "CAGe++.vcf", "VCF_output_file");
    TCLAP::ValueArg<string> input_SNP_db_arg("s", "input_SNP_db", "Filename of sqlite3 SNP database", true, "", "SNP_input_db");
    cmd.xorAdd(output_VCF_arg, input_SNP_db_arg);

    cmd.parse(argc, argv);

    parsed_params.dataloader_params.contig = contig_arg.getValue();
    parsed_params.dataloader_params.start = start_arg.getValue();
    parsed_params.dataloader_params.end = end_arg.getValue();
    parsed_params.dataloader_params.step = step_arg.getValue();
    parsed_params.beta = beta_arg.getValue();
    parsed_params.output_file = output_file_arg.getValue();
    parsed_params.dataloader_params.verbose = verbose_switch_arg.getValue();

    if (output_VCF_arg.isSet()) {
      parsed_params.dataloader_params.call_variants = true;
      parsed_params.dataloader_params.output_VCF = output_VCF_arg.getValue();
    } else {
      parsed_params.dataloader_params.call_variants = false;
      parsed_params.dataloader_params.input_SNP_db = input_SNP_db_arg.getValue();
    }
  } catch (TCLAP::ArgException &e) {
    cerr << "Error: " << e.error() << " for arg " << e.argId() << std::endl;
    throw runtime_error("Could not parse command line arguments");
  }

  if (parsed_params.dataloader_params.start >= parsed_params.dataloader_params.end) {
    throw std::runtime_error("Start position must be less than end position");
  }

  if (parsed_params.beta < 0) {
    throw std::runtime_error("Error: beta must be nonnegative");
  }

  return parsed_params;
}	

double wall_time ()
{
  struct timeval t;
  gettimeofday (&t, NULL);
  return 1.*t.tv_sec + 1.e-6*t.tv_usec;
}


int main(int argc, char** argv) {

  // To speed up calls to getline
  ios_base::sync_with_stdio(false);

  try {
    // Parse command line parameters
    CageParams params = parse_parameters(argc, argv);

    double begin_t = wall_time();

    string ref;
    // Get the reference (first line of stdin)
    cout << "Reading the complete reference ... ";
    getline(cin, ref);
    cout << "Done" << endl;
    Dataloader data_loader(params.dataloader_params, move(ref)); 	

    // Run pelt
    pelt(data_loader, params);
    double end_t = wall_time();

    cout << "Time to compute cage: " << end_t - begin_t << " seconds\n";

  } catch (const std::exception& e) {
    cerr << e.what() << endl; 
    return 1;
  }

  return 0;
}
