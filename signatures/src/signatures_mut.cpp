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

std::pair<bool, std::string> same_reeds(std::string A_count, std::string G_count,
                                        std::string C_count, std::string T_count,
                                        std::string a_count, std::string g_count,
                                        std::string c_count, std::string t_count,
                                        int total_coverage);

std::string multiple_occur(std::string A_count, std::string G_count, std::string C_count, std::string T_count,
                           std::string a_count, std::string g_count, std::string c_count, std::string t_count, 
                           std::string A_s, std::string G_s, std::string C_s, std::string T_s);
///////////////////////////////////////////////////////////////////////////////////////////////////////

int main (int argc, char* argv[])
{

  if (argc != 6)
  {
    std::cout << "Use: " << argv[0] << " [lines_mutations.dat] [lines_control.dat] [lines_evolved.dat] [output_path] [copy_number]" << std::endl; 
    return 1;
  }

  std::ofstream file_sign_p(std::string{argv[4]}+"signatures_p"+std::string{argv[5]}+".dat", std::ios_base::binary);
  if (!file_sign_p.is_open())
  {
    std::cerr << "File file_signatures.dat not opened. Exit" << std::endl;
    return 1;
  }

  /////////////////////////////////////////////////////////////////////////////////////////////////////

  auto start_timer = high_resolution_clock::now(); //reference base 
  std::string A_s = "A"; 
  std::string G_s = "G";
  std::string C_s = "C";
  std::string T_s = "T";

  //file mutations
  std::string mutation_line;
  std::string chr_m;
  std::string base_num_m;
  std::string anc_m;
  std::string mut_m;
  std::string cov_endpoint_m;
  std::string mut_read_m;
  std::string CN_m;

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

  //ancestor
  std::string ancestor_line;
  std::string chr_a;
  std::string base_num_a;
  std::string ref_a;
  std::string A_count_a;
  std::string G_count_a;
  std::string C_count_a;
  std::string T_count_a;
  std::string a_count_a;
  std::string g_count_a;
  std::string c_count_a;
  std::string t_count_a;

  std::ifstream in_file_mutations(std::string{argv[1]}, std::ios_base::binary);
  if (!in_file_mutations.is_open())
  {
    std::cerr << "File mutations.dat not opened. Exit" << std::endl;
    return 1;
  }

  std::ifstream in_file_ancestor(std::string{argv[2]}, std::ios_base::binary);
  if (!in_file_ancestor.is_open())
  {
    std::cerr << "File ancestor.dat not opened. Exit" << std::endl;
    return 1;
  }
 
  std::ifstream in_file_endpoint(std::string{argv[3]}, std::ios_base::binary);
  if (!in_file_endpoint.is_open())
  {
    std::cerr << "File endpoint.dat not opened. Exit" << std::endl;
    return 1;
  }

  ///////////////////////////////////////////////////////////////////////

  size_t lines_mut = 0;
  auto oldpos_end = in_file_endpoint.tellg();
  auto oldpos_anc = in_file_ancestor.tellg();

  while (std::getline(in_file_mutations, mutation_line))
  {
    //mutations
    lines_mut += 1;
    std::stringstream data_mutation(mutation_line);
    data_mutation >> chr_m >> base_num_m >> anc_m >> mut_m >> cov_endpoint_m >> mut_read_m >> CN_m; 

    std::string prev_e;
    std::string curr_e;
    std::string next_e;

    std::string ref_p;
    std::string ref_c;
    std::string ref_n;

    std::string prev_a;
    std::string curr_a;
    std::string next_a;

    //std::cout << "start prev endpoint" << std::endl;
    //endpoint
    while (std::getline(in_file_endpoint, endpoint_line))
    {
      //prev
      std::stringstream data_endpoint(endpoint_line);
      data_endpoint >> chr_e >> base_num_e >> ref_e >> A_count_e >> G_count_e >> C_count_e >> T_count_e >>
                                                       a_count_e >> g_count_e >> c_count_e >> t_count_e; 
      int endpoint_forward = std::stoi(A_count_e) + std::stoi(G_count_e) + std::stoi(C_count_e) + std::stoi(T_count_e);
      int endpoint_reverse = std::stoi(a_count_e) + std::stoi(g_count_e) + std::stoi(c_count_e) + std::stoi(t_count_e);
      int endpoint_coverage = endpoint_forward + endpoint_reverse;

      //std::cout << chr_e << "\t" << base_num_e << "\t" << ref_e << "\t" << A_count_e << "\t" << G_count_e << "\t" << C_count_e << "\t" << T_count_e << "\t" << a_count_e << "\t" << g_count_e << "\t" << c_count_e << "\t" << t_count_e << std::endl;
      if ((std::stoi(chr_e) == std::stoi(chr_m)) && (std::stoi(base_num_e) == (std::stoi(base_num_m)-1)))
      {
        oldpos_end = in_file_endpoint.tellg();
        ref_p = ref_e; 
        std::pair<bool, std::string> same_reed = same_reeds(A_count_e, G_count_e, C_count_e, T_count_e, a_count_e, g_count_e, c_count_e, t_count_e, endpoint_coverage);
        if (same_reed.first == true)
        {
          prev_e = same_reed.second; 
          break;
        }
        else
        {
          prev_e = multiple_occur(A_count_e, G_count_e, C_count_e, T_count_e, a_count_e, g_count_e, c_count_e, t_count_e, A_s, G_s, C_s, T_s);
          break;
        }
        break;
      }

      if ((std::stoi(chr_e) == std::stoi(chr_m)) && (std::stoi(base_num_e) > (std::stoi(base_num_m)-1)))
      {
        //oldpos_end = in_file_endpoint.tellg();
        break;
      }
      
      if (std::stoi(chr_e) > std::stoi(chr_m))
      {
        //oldpos_end = in_file_endpoint.tellg();
        break;
      }
    } 

    //std::cout << "end prev endpoint" << std::endl;

    //in_file_endpoint.seekg(oldpos_end);
    //in_file_endpoint.clear();
    //in_file_endpoint.seekg(0, std::ios::beg);

    //oldpos_end = in_file_endpoint.tellg();

    //std::cout << "start curr endpoint" << std::endl;
    while (std::getline(in_file_endpoint, endpoint_line))
    {
      //curr
      std::stringstream data_endpoint(endpoint_line);
      data_endpoint >> chr_e >> base_num_e >> ref_e >> A_count_e >> G_count_e >> C_count_e >> T_count_e >>
                                                       a_count_e >> g_count_e >> c_count_e >> t_count_e; 

      //std::cout << chr_e << "\t" << base_num_e << "\t" << ref_e << "\t" << A_count_e << "\t" << G_count_e << "\t" << C_count_e << "\t" << T_count_e << "\t" << a_count_e << "\t" << g_count_e << "\t" << c_count_e << "\t" << t_count_e << std::endl;
      if ((std::stoi(chr_e) == std::stoi(chr_m)) && (std::stoi(base_num_e) == std::stoi(base_num_m)))
      {
        curr_e = mut_m; 
        ref_c = ref_e; 
        break;
      }

      if ((std::stoi(chr_e) == std::stoi(chr_m)) && (std::stoi(base_num_e) > std::stoi(base_num_m)))
      {
        //oldpos_end = in_file_endpoint.tellg();
        break;
      }
      
      if (std::stoi(chr_e) > std::stoi(chr_m))
      {
        //oldpos_end = in_file_endpoint.tellg();
        break;
      }
    }
    //std::cout << "end curr endpoint" << std::endl;

    //in_file_endpoint.seekg(oldpos_end);
    //in_file_endpoint.clear();
    //in_file_endpoint.seekg(0, std::ios::beg);

    //oldpos_end = in_file_endpoint.tellg();

    //std::cout << "start next endpoint" << std::endl;
    while (std::getline(in_file_endpoint, endpoint_line))
    {
      //next
      std::stringstream data_endpoint(endpoint_line);
      data_endpoint >> chr_e >> base_num_e >> ref_e >> A_count_e >> G_count_e >> C_count_e >> T_count_e >>
                                                       a_count_e >> g_count_e >> c_count_e >> t_count_e; 
      int endpoint_forward = std::stoi(A_count_e) + std::stoi(G_count_e) + std::stoi(C_count_e) + std::stoi(T_count_e);
      int endpoint_reverse = std::stoi(a_count_e) + std::stoi(g_count_e) + std::stoi(c_count_e) + std::stoi(t_count_e);
      int endpoint_coverage = endpoint_forward + endpoint_reverse;

      //std::cout << chr_e << "\t" << base_num_e << "\t" << ref_e << "\t" << A_count_e << "\t" << G_count_e << "\t" << C_count_e << "\t" << T_count_e << "\t" << a_count_e << "\t" << g_count_e << "\t" << c_count_e << "\t" << t_count_e << std::endl;
      if ((std::stoi(chr_e) == std::stoi(chr_m)) && (std::stoi(base_num_e) == (std::stoi(base_num_m)+1)))
      {
        ref_n = ref_e; 
        std::pair<bool, std::string> same_reed = same_reeds(A_count_e, G_count_e, C_count_e, T_count_e, a_count_e, g_count_e, c_count_e, t_count_e, endpoint_coverage);
        if (same_reed.first == true)
        {
          next_e = same_reed.second; 
          break;
        }
        else
        {
          next_e = multiple_occur(A_count_e, G_count_e, C_count_e, T_count_e, a_count_e, g_count_e, c_count_e, t_count_e, A_s, G_s, C_s, T_s);
          break;
        }
        break;
      }

      if ((std::stoi(chr_e) == std::stoi(chr_m)) && (std::stoi(base_num_e) > (std::stoi(base_num_m)+1)))
      {
        //oldpos_end = in_file_endpoint.tellg();
        break;
      }
      
      if (std::stoi(chr_e) > std::stoi(chr_m))
      {
        //oldpos_end = in_file_endpoint.tellg();
        break;
      }
    }
    //std::cout << "end next endpoint" << std::endl;

    //in_file_endpoint.clear();
    //in_file_endpoint.seekg(0, std::ios::beg);


    //std::cout << "start prev ancestor" << std::endl;
    //ancestor
    while (std::getline(in_file_ancestor, ancestor_line))
    {
      //prev
      std::stringstream data_ancestor(ancestor_line);
      data_ancestor >> chr_a >> base_num_a >> ref_a >> A_count_a >> G_count_a >> C_count_a >> T_count_a >>
                                                       a_count_a >> g_count_a >> c_count_a >> t_count_a; 
      int ancestor_forward = std::stoi(A_count_a) + std::stoi(G_count_a) + std::stoi(C_count_a) + std::stoi(T_count_a);
      int ancestor_reverse = std::stoi(a_count_a) + std::stoi(g_count_a) + std::stoi(c_count_a) + std::stoi(t_count_a);
      int ancestor_coverage = ancestor_forward + ancestor_reverse;

      if ((std::stoi(chr_a) == std::stoi(chr_m)) && (std::stoi(base_num_a) == (std::stoi(base_num_m)-1)))
      {
        oldpos_anc = in_file_ancestor.tellg();
        std::pair<bool, std::string> same_reed = same_reeds(A_count_a, G_count_a, C_count_a, T_count_a, a_count_a, g_count_a, c_count_a, t_count_a, ancestor_coverage);
        if (same_reed.first == true)
        {
          prev_a = same_reed.second; 
          break;
        }
        else
        {
          prev_a = multiple_occur(A_count_a, G_count_a, C_count_a, T_count_a, a_count_a, g_count_a, c_count_a, t_count_a, A_s, G_s, C_s, T_s);
          break;
        }
        break;
      }

      if ((std::stoi(chr_a) == std::stoi(chr_m)) && (std::stoi(base_num_a) > (std::stoi(base_num_m)-1)))
      {
        //oldpos_anc = in_file_ancestor.tellg();
        break;
      }
      
      if (std::stoi(chr_a) > std::stoi(chr_m))
      {
        //oldpos_end = in_file_ancestor.tellg();
        break;
      }
    }
    //std::cout << "end prev ancestor" << std::endl;

    //in_file_ancestor.seekg(oldpos_anc);
    //in_file_ancestor.clear();
    //in_file_ancestor.seekg(0, std::ios::beg);

    //std::cout << "start curr ancestor" << std::endl;
    while (std::getline(in_file_ancestor, ancestor_line))
    {
      //curr
      std::stringstream data_ancestor(ancestor_line);
      data_ancestor >> chr_a >> base_num_a >> ref_a >> A_count_a >> G_count_a >> C_count_a >> T_count_a >>
                                                       a_count_a >> g_count_a >> c_count_a >> t_count_a; 
      if ((std::stoi(chr_a) == std::stoi(chr_m)) && (std::stoi(base_num_a) == std::stoi(base_num_m)))
      {
        curr_a = anc_m; 
        break;
      }

      if ((std::stoi(chr_a) == std::stoi(chr_m)) && (std::stoi(base_num_a) > std::stoi(base_num_m)))
      {
        //oldpos_end = in_file_endpoint.tellg();
        break;
      }
      
      if (std::stoi(chr_a) > std::stoi(chr_m))
      {
        //oldpos_end = in_file_endpoint.tellg();
        break;
      }
    }
    //std::cout << "end curr ancestor" << std::endl;

    //in_file_ancestor.seekg(oldpos_anc);
    //in_file_ancestor.clear();
    //in_file_ancestor.seekg(0, std::ios::beg);

    //std::cout << "start next ancestor" << std::endl;
    while (std::getline(in_file_ancestor, ancestor_line))
    {
      //next
      std::stringstream data_ancestor(ancestor_line);
      data_ancestor >> chr_a >> base_num_a >> ref_a >> A_count_a >> G_count_a >> C_count_a >> T_count_a >>
                                                       a_count_a >> g_count_a >> c_count_a >> t_count_a; 
      int ancestor_forward = std::stoi(A_count_a) + std::stoi(G_count_a) + std::stoi(C_count_a) + std::stoi(T_count_a);
      int ancestor_reverse = std::stoi(a_count_a) + std::stoi(g_count_a) + std::stoi(c_count_a) + std::stoi(t_count_a);
      int ancestor_coverage = ancestor_forward + ancestor_reverse;

      if ((std::stoi(chr_a) == std::stoi(chr_m)) && (std::stoi(base_num_a) == (std::stoi(base_num_m)+1)))
      {
        std::pair<bool, std::string> same_reed = same_reeds(A_count_a, G_count_a, C_count_a, T_count_a, a_count_a, g_count_a, c_count_a, t_count_a, ancestor_coverage);
        if (same_reed.first == true)
        {
          next_a = same_reed.second; 
          break;
        }
        else
        {
          next_a = multiple_occur(A_count_a, G_count_a, C_count_a, T_count_a, a_count_a, g_count_a, c_count_a, t_count_a, A_s, G_s, C_s, T_s);
          break;
        }
        break;
      }

      if ((std::stoi(chr_a) == std::stoi(chr_m)) && (std::stoi(base_num_a) > (std::stoi(base_num_m)+1)))
      {
        //oldpos_anc = in_file_ancestor.tellg();
        break;
      }
      
      if (std::stoi(chr_a) > std::stoi(chr_m))
      {
        //oldpos_end = in_file_ancestor.tellg();
        break;
      }
    }
    //std::cout << "end next ancestor" << std::endl;
    //in_file_ancestor.clear();
    //in_file_ancestor.seekg(0, std::ios::beg);

    if (!prev_e.empty() && !curr_e.empty() && !next_e.empty() && !prev_a.empty() && !curr_a.empty() && !next_a.empty())
    {
      file_sign_p << chr_m << "\t" << base_num_m << "\t" << prev_a << "\t" << curr_a << "\t" << next_a 
                                                 << "\t" << prev_e << "\t" << curr_e << "\t" << next_e 
                                                 << "\t" << ref_p << "\t" << ref_c << "\t" << ref_n << std::endl;
    }
    //else
    //{
    //  std::cout << chr_m << "\t" << base_num_m << "\t" << "Not found!" << std::endl;
    //}
     
    in_file_endpoint.seekg(oldpos_end);
    in_file_ancestor.seekg(oldpos_anc);

    if (lines_mut % 100000 == 0)
    {
      std::cout << lines_mut << "\t" << chr_m << std::endl; 
    }
  }

  in_file_mutations.close();
  in_file_endpoint.close();
  in_file_ancestor.close();
  file_sign_p.close();

  auto stop_timer = high_resolution_clock::now();
  auto measure_time = duration_cast<minutes>(stop_timer - start_timer);

  std::cout << "Execution time:  " << measure_time.count() << "  min." << std::endl;

  return 0;
}

std::pair<bool, std::string> same_reeds(std::string A_count, std::string G_count,
                                         std::string C_count, std::string T_count,
                                         std::string a_count, std::string g_count,
                                         std::string c_count, std::string t_count,
                                         int total_coverage)
{
  std::pair<bool, std::string> ret_val = std::make_pair(false, "NN");

  int ATOT = std::stoi(A_count) + std::stoi(a_count);
  int GTOT = std::stoi(G_count) + std::stoi(g_count);
  int CTOT = std::stoi(C_count) + std::stoi(c_count);
  int TTOT = std::stoi(T_count) + std::stoi(t_count);

  if (ATOT == total_coverage)
  {
    ret_val.first = true;
    ret_val.second = "A";
    return ret_val;
  }

  if (GTOT == total_coverage)
  {
    ret_val.first = true;
    ret_val.second = "G";
    return ret_val;
  }

  if (CTOT == total_coverage)
  {
    ret_val.first = true;
    ret_val.second = "C";
    return ret_val;
  }

  if (TTOT == total_coverage)
  {
    ret_val.first = true;
    ret_val.second = "T";
    return ret_val;
  }

  return ret_val;
}

std::string multiple_occur(std::string A_count, std::string G_count, std::string C_count, std::string T_count,
                           std::string a_count, std::string g_count, std::string c_count, std::string t_count,
                           std::string A_s, std::string G_s, std::string C_s, std::string T_s)
{
  std::string ret_occur;  
  int ATOT = std::stoi(A_count) + std::stoi(a_count);
  int GTOT = std::stoi(G_count) + std::stoi(g_count);
  int CTOT = std::stoi(C_count) + std::stoi(c_count);
  int TTOT = std::stoi(T_count) + std::stoi(t_count);

  std::vector<int> numbers;
  numbers.push_back(ATOT);
  numbers.push_back(GTOT);
  numbers.push_back(CTOT);
  numbers.push_back(TTOT);

  std::sort(std::begin(numbers), std::end(numbers));

  for (auto el : numbers)
  {
    if (el != 0)
    {
      if (el == ATOT)
      {
        ret_occur += A_s;
        continue;
      }
      else if (el == GTOT)
      {
        ret_occur += G_s;
        continue;
      }
      else if (el == CTOT)
      {
        ret_occur += C_s;
        continue;
      }
      else
      {
        ret_occur += T_s;
        continue;
      }
    }
  }

  return ret_occur;
}

