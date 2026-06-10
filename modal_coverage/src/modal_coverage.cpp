#include <iostream>
#include <iomanip>
#include <string>
#include <cstring>
#include <fstream>
#include <sstream>
#include <chrono>
#include <cmath>
#include <vector>
#include <set>
#include <algorithm>
#include <cassert>
#include <utility>

using namespace std::chrono;

void distribution(std::vector<int> &coverages, std::vector<std::pair<int,int>> &counts);
std::pair<int, int>  modal(std::vector<std::pair<int,int>> &counts);

int main (int argc, char* argv[])
{

  if (argc < 3)
  {
    std::cout << "Use: " << argv[0] << " <path_pseudo-pileup.dat> <path_output_file.dat> " << std::endl;
    return 1;
  }

  std::string seq_file = std::string{argv[1]};
  std::string output_path_file = std::string{argv[2]};
  //input file
  std::ifstream in_file(seq_file, std::ios_base::binary);
  if (!in_file.is_open())
  {
    std::cerr << "File infile.dat not opened. Exit" << std::endl;
    return 1;
  }

  //output file
  std::ofstream modal_out(output_path_file, std::ios_base::binary);
  if (!modal_out.is_open())
  {
    std::cerr << "File modal_out.dat not opened. Exit" << std::endl;
    return 1;
  }

  auto start_timer = high_resolution_clock::now();

  std::vector<int> coverages;
  std::vector<std::pair<int,int>> counts;

  //reference base
  std::string A_s = "A";
  std::string G_s = "G";
  std::string C_s = "C";
  std::string T_s = "T";

  //ancestor
  std::string line;
  std::string chr;
  std::string base_num;
  std::string ref;
  std::string A_count;
  std::string G_count;
  std::string C_count;
  std::string T_count;
  std::string a_count;
  std::string g_count;
  std::string c_count;
  std::string t_count;
  int sequencing_coverage = 0;
  size_t lines = 0;

  while (std::getline(in_file, line))
  {
    std::stringstream data_seq(line);
    lines += 1;
    data_seq >> chr >> base_num >> ref >> A_count >> G_count >> C_count >> T_count >> 
                                          a_count >> g_count >> c_count >> t_count; 

    sequencing_coverage = std::stoi(A_count) + std::stoi(G_count) + std::stoi(C_count) + std::stoi(T_count) +
                          std::stoi(a_count) + std::stoi(g_count) + std::stoi(c_count) + std::stoi(t_count);

    coverages.push_back(sequencing_coverage);

    if (lines % 10000000 == 0)
    {
      std::cout << lines << "\t" << chr << std::endl;
    }
  }

  distribution(coverages, counts);
 
  std::pair<int, int> modal_cov = modal(counts); 
  
  std::cout << modal_cov.first << "\t" << modal_cov.second << std::endl; 
 
  modal_out << modal_cov.first << std::endl;
  modal_out << modal_cov.second << std::endl;

  modal_out.close();

  auto stop_timer = high_resolution_clock::now();
  auto measure_time = duration_cast<minutes>(stop_timer - start_timer);

  std::cout << "Execution time:  " << measure_time.count() << "  min." << std::endl;

  return 0;
}

void distribution(std::vector<int> &coverages, std::vector<std::pair<int,int>> &counts)
{
  std::sort(coverages.begin(), coverages.end());

  std::set<int> temp_set;

  for (auto &el : coverages)
  {
    temp_set.insert(el);
  }
  
  for (const auto &el : temp_set)
  {
    int mycount = std::count(coverages.begin(), coverages.end(), el);
    auto p_count = std::make_pair(el, mycount);
    counts.push_back(p_count);
  }
}

std::pair<int, int> modal(std::vector<std::pair<int,int>> &counts)
{
  std::sort(counts.begin(), counts.end(),
            [] (const std::pair<int, int> &a,
                const std::pair<int, int> &b) {
                return a.second > b.second;});

  return counts[0];
}
