#include <iostream>
#include <limits>
#include <iomanip>
#include <string>
#include <cstring>
#include <fstream>
#include <sstream>
#include <chrono>
#include <cmath>
#include <vector>
#include <tuple>
#include <algorithm>

using namespace std::chrono;

int rounded(double val);

int main (int argc, char* argv[])
{
  //int coverage_thr_ancestor = 0;
  //int std_dev_cov_ancestor1 = 0;
  //int std_dev_cov_ancestor2 = 0;
  //int max_cov_a = 0;
  //int min_cov_a = 0;
  //int coverage_thr_endpoint = 0;
  //int std_dev_cov_endpoint1 = 0;
  //int std_dev_cov_endpoint2 = 0;
  //int max_cov_e = 0;
  //int min_cov_e = 0;

  if (argc != 4)
  {
    std::cout << "Use: " << argv[0] << "[common_lines_evolved.dat] [output_path] [copy_number]" << std::endl; 
    return 1;
  }

  std::ofstream file_reeds_p(std::string{argv[2]}+"file_reeds_p"+std::string{argv[3]}+".dat", std::ios_base::binary);
  if (!file_reeds_p.is_open())
  {
    std::cerr << "File file_reeds.dat not opened. Exit" << std::endl;
    return 1;
  }

  std::ofstream file_sel_p(std::string{argv[2]}+"file_no_mutations_p"+std::string{argv[3]}+".dat", std::ios_base::binary);
  if (!file_sel_p.is_open())
  {
    std::cerr << "File file_no_mutations.dat not opened. Exit" << std::endl;
    return 1;
  }


  std::ifstream in_file_endpoint(std::string{argv[1]}, std::ios_base::binary);
  if (!in_file_endpoint.is_open())
  {
    std::cerr << "File endpoint.dat not opened. Exit" << std::endl;
    return 1;
  }
  /////////////////////////////////////////////////////////////////////////////////////////////////////

  auto start_timer = high_resolution_clock::now(); //reference base 
  std::string A_s = "A"; 
  std::string G_s = "G";
  std::string C_s = "C";
  std::string T_s = "T";

  //endpoint
  std::string endpoint_line;
  std::string chr_e;
  std::string base_num_e;
  std::string ref_e;
  std::string A_count_e;
  std::string G_count_e;
  std::string C_count_e;
  std::string T_count_e;
  std::string a_count_e;
  std::string g_count_e;
  std::string c_count_e;
  std::string t_count_e;
  int endpoint_coverage = 0;
  int endpoint_forward = 0;
  int endpoint_reverse = 0;
  double endpoint_mut_freq = 0.0;
  std::string endpoint_mutated_in;
  size_t lines_endpoint = 0;
 
  size_t selected_bases_p = 0;
  size_t n_bases_CN_mut = 0;

  double CN_end = std::stod(argv[3]);
  //coverage_thr_ancestor = rounded(std::stod(argv[4]));
  //std_dev_cov_ancestor1 = rounded(std::stod(argv[5]));
  //std_dev_cov_ancestor2 = rounded(std::stod(argv[5]));
  //min_cov_a = coverage_thr_ancestor - std_dev_cov_ancestor1;
  //max_cov_a = coverage_thr_ancestor + std_dev_cov_ancestor2;

  //coverage_thr_endpoint = rounded(std::stod(argv[6]));
  //std_dev_cov_endpoint1 = rounded(std::stod(argv[7]));
  //std_dev_cov_endpoint2 = rounded(std::stod(argv[7]));
  //min_cov_e = coverage_thr_endpoint - std_dev_cov_endpoint1;
  //max_cov_e = coverage_thr_endpoint + std_dev_cov_endpoint2;
  ///////////////////////////////////////////////////////////////////////

  std::cout << "CN = " << CN_end << std::endl;

  while (std::getline(in_file_endpoint, endpoint_line))
  {
    //endpoint
    int mutated_reads = 0;
    lines_endpoint += 1;
    std::stringstream data_endpoint(endpoint_line);
    data_endpoint >> chr_e >> base_num_e >> ref_e >> A_count_e >> G_count_e >> C_count_e >> T_count_e >>
                                                     a_count_e >> g_count_e >> c_count_e >> t_count_e; 

    if (ref_e == A_s || ref_e == G_s || ref_e == C_s || ref_e == T_s)
    { 
      selected_bases_p += 1;
      //endpoint coverage
      endpoint_forward = std::stoi(A_count_e) + std::stoi(G_count_e) + std::stoi(C_count_e) + std::stoi(T_count_e);
      endpoint_reverse = std::stoi(a_count_e) + std::stoi(g_count_e) + std::stoi(c_count_e) + std::stoi(t_count_e);
      endpoint_coverage = endpoint_forward + endpoint_reverse;

      if (ref_e == A_s)
      {
        mutated_reads = std::max(std::stoi(T_count_e)+std::stoi(t_count_e),
                                 std::stoi(C_count_e)+std::stoi(c_count_e)); 
        mutated_reads = std::max(mutated_reads, std::stoi(G_count_e)+std::stoi(g_count_e));
        endpoint_mut_freq = static_cast<double>(mutated_reads)/static_cast<double>(endpoint_coverage);
        if (mutated_reads != 0)
        {
          if (mutated_reads == (std::stoi(G_count_e)+std::stoi(g_count_e))) { endpoint_mutated_in = "G";}
          if (mutated_reads == (std::stoi(C_count_e)+std::stoi(c_count_e))) { endpoint_mutated_in = "C";}
          if (mutated_reads == (std::stoi(T_count_e)+std::stoi(t_count_e))) { endpoint_mutated_in = "T";}
        }
        else
        {
          endpoint_mutated_in = "NN";
        }
      } 

      if (ref_e == G_s)
      {
        mutated_reads = std::max(std::stoi(T_count_e)+std::stoi(t_count_e),
                                 std::stoi(C_count_e)+std::stoi(c_count_e)); 
        mutated_reads = std::max(mutated_reads, std::stoi(A_count_e)+std::stoi(a_count_e));
        endpoint_mut_freq = static_cast<double>(mutated_reads)/static_cast<double>(endpoint_coverage);
        if (mutated_reads != 0)
        {
          if (mutated_reads == (std::stoi(A_count_e)+std::stoi(a_count_e))) { endpoint_mutated_in = "A";}
          if (mutated_reads == (std::stoi(C_count_e)+std::stoi(c_count_e))) { endpoint_mutated_in = "C";}
          if (mutated_reads == (std::stoi(T_count_e)+std::stoi(t_count_e))) { endpoint_mutated_in = "T";}
        }
        else
        {
          endpoint_mutated_in = "NN";
        }
      }

      if (ref_e == C_s)
      {
        mutated_reads = std::max(std::stoi(T_count_e)+std::stoi(t_count_e),
                                 std::stoi(A_count_e)+std::stoi(a_count_e)); 
        mutated_reads = std::max(mutated_reads, std::stoi(G_count_e)+std::stoi(g_count_e));
        endpoint_mut_freq = static_cast<double>(mutated_reads)/static_cast<double>(endpoint_coverage);
        if (mutated_reads != 0)
        {
          if (mutated_reads == (std::stoi(A_count_e)+std::stoi(a_count_e))) { endpoint_mutated_in = "A";}
          if (mutated_reads == (std::stoi(G_count_e)+std::stoi(g_count_e))) { endpoint_mutated_in = "G";}
          if (mutated_reads == (std::stoi(T_count_e)+std::stoi(t_count_e))) { endpoint_mutated_in = "T";}
        }
        else
        {
          endpoint_mutated_in = "NN";
        }
      }

      if (ref_e == T_s)
      {
        mutated_reads = std::max(std::stoi(A_count_e)+std::stoi(a_count_e),
                                 std::stoi(C_count_e)+std::stoi(c_count_e)); 
        mutated_reads = std::max(mutated_reads, std::stoi(G_count_e)+std::stoi(g_count_e));
        endpoint_mut_freq = static_cast<double>(mutated_reads)/static_cast<double>(endpoint_coverage);
        if (mutated_reads != 0)
        {
          if (mutated_reads == (std::stoi(A_count_e)+std::stoi(a_count_e))) { endpoint_mutated_in = "A";}
          if (mutated_reads == (std::stoi(G_count_e)+std::stoi(g_count_e))) { endpoint_mutated_in = "G";}
          if (mutated_reads == (std::stoi(C_count_e)+std::stoi(c_count_e))) { endpoint_mutated_in = "C";}
        }
        else
        {
          endpoint_mutated_in = "NN";
        }
      }
       
      if (endpoint_mut_freq > 0.0)
      { 
        n_bases_CN_mut += 1; 

        std::string t_chr_e;
        for( char c : chr_e ) if (std::isalnum(c)) t_chr_e += c;

        file_reeds_p << t_chr_e << "\t" << base_num_e << "\t" << ref_e << "\t" << endpoint_mutated_in << "\t" << endpoint_coverage << "\t" << mutated_reads << "\t" << CN_end << std::endl;
      }
    }

    if (lines_endpoint % 10000000 == 0)
    {
      std::cout << lines_endpoint << "\t" << chr_e << std::endl; 
    }
  }
 
  in_file_endpoint.close();

  size_t no_mutations_p = selected_bases_p - n_bases_CN_mut;

  std::cout << "Selected bases: " << selected_bases_p << std::endl;
  std::cout << "Not mutated number of bases: " << no_mutations_p << std::endl;

  file_sel_p << no_mutations_p << std::endl;

  file_reeds_p.close();
  file_sel_p.close();

  auto stop_timer = high_resolution_clock::now();
  auto measure_time = duration_cast<minutes>(stop_timer - start_timer);

  std::cout << "Execution time:  " << measure_time.count() << "  min." << std::endl;

  return 0;
}
