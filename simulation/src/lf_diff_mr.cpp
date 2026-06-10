#include <iostream>
#include <filesystem>
#include <fstream>
#include <string>
#include <sstream>
#include <chrono>
#include <vector>
#include <random>
#include "tree.h"
#include "random.h"
#include "XoshiroCpp.hpp"

using namespace std::chrono;

int main(int argc, char* argv[])
{
 
  if (argc != 3)
  {
    std::cout  << "Use: " << argv[0] << "  data.dat   path_output_dir" << std::endl;
    return 1;
  }

  std::ifstream data_file(std::string{argv[1]});
  if (!data_file.is_open())
  {
    std::cerr << "File data.dat not opened. Exit" << std::endl;
    return 1;
  }

  Random rnd(SEED "/Primes", SEED "/seed.in");
  for (size_t i = 0; i < 100000; ++i)
  {
    rnd.rannyu();
  } 

  std::string seq_name;
  size_t n_blocks = 0;
  size_t n_steps = 0; 
  size_t n_bases = 0; //10^8 or 10^9
  size_t coverage = 0;
  size_t ploidy = 1;
  size_t n_zero = 0;
  double death_probability = 0.0;
  double seq_error = 0.0;
  double K = 0.0;
  size_t min_survived_nodes_number = 0;
  size_t max_generations = 0;
  long double MUT_rate = 0.0;
  size_t clonal_mutations = 0;
  long double MUT_rate_higher = 0.0;
  int n_levels = 0;
  int layer_id = 0;
  data_file >> seq_name >> n_blocks >> n_steps >> n_bases >> 
               coverage >> ploidy >> n_zero >> death_probability >> 
               seq_error >> K >> min_survived_nodes_number >> 
               max_generations >> MUT_rate >> clonal_mutations >>
               MUT_rate_higher >> n_levels >> layer_id;

  std::cout << seq_name << std::endl;
  std::cout << "Number of simulation blocks: " << n_blocks << std::endl;
  std::cout << "Number of simulation steps: " << n_steps << std::endl;
  std::cout << "Number of bases: " << n_bases << std::endl;
  std::cout << "Coverage: " << coverage << std::endl;
  std::cout << "Ploidy: " << ploidy << std::endl;
  std::cout << "N0: " << n_zero << std::endl;
  std::cout << "Death probability: " << death_probability << std::endl;
  std::cout << "Sequencing error: " << seq_error << std::endl;
  std::cout << "K: " << K << std::endl;
  std::cout << "Min number of survived cells: " << min_survived_nodes_number << std::endl;
  std::cout << "Max number of generations: " << max_generations << std::endl;
  std::cout << "Mutation rate: " << MUT_rate << std::endl;
  std::cout << "Number of clonal mutations = " << clonal_mutations << std::endl;
  std::cout << "Hyper mutation rate = " << MUT_rate_higher << std::endl;
  std::cout << "Minimum number of Macro levels = " << n_levels << std::endl;
  std::cout << "Layer selected = " << layer_id << std::endl;

  auto start_timer = high_resolution_clock::now();

  constexpr std::uint64_t seed_mutations = 12345;
  XoshiroCpp::Xoshiro256PlusPlus rng(seed_mutations);
  auto seed_r = rng(); 
  static std::mt19937 engine(seed_r);

  std::vector<Node<size_t>> sequential_tree; 
  Tree<size_t> tree(sequential_tree, death_probability, seq_error, K, 
                    max_generations, n_bases, coverage, ploidy, n_zero);
  tree.warm_up_random();
  tree.lambda_clonal = clonal_mutations;

  for (size_t j = 0; j < n_steps; ++j)
  {
    std::cout << "Step = " << j << std::endl;

    std::string file_freqs = std::string{argv[2]}+std::string{"/"}+"file_frequencies_"+std::to_string(j)+".dat";
    std::string file_nodes = std::string{argv[2]}+std::string{"/"}+"file_n_nodes_"+std::to_string(j)+".dat";

    //if (std::filesystem::exists(file_freqs) && std::filesystem::exists(file_nodes))
    if (std::filesystem::exists(file_freqs))
    {
      std::cout << "File " << file_freqs << " exists in the output folder." << std::endl;
      std::cout << "Skip" << std::endl;
      tree.warm_up_random();
      tree.r_number(engine);
      continue;
    }
    else
    {
      std::ofstream file_out_freqs(file_freqs, std::ios_base::binary);
      std::ofstream file_out_nodes(file_nodes, std::ios_base::binary);
      tree.start_tree(MUT_rate, min_survived_nodes_number, max_generations);
      tree.create_macro_levels(n_levels);
      tree.set_mutations_diff_mr(engine, MUT_rate_higher, layer_id);
      //std::cout << tree;
      for (const auto& el : tree.freq_mutations)
      {
        file_out_freqs << el << std::endl;
      }
      
      file_out_nodes << tree.get_survived_nodes() << std::endl;

      file_out_nodes.close();
      file_out_freqs.close();
    }
  }
  
  auto stop_timer = high_resolution_clock::now();
  auto measure_time = duration_cast<minutes>(stop_timer - start_timer);
  std::cout << "Execution time:  " << measure_time.count() << "  min." << std::endl;
  return 0;
}
