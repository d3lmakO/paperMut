#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <iomanip>
#include <chrono>

std::string get_right_of_delim(std::string const& str, std::string const& delim)
{
  return str.substr(str.find(delim) + delim.size());
}

int main(int argc, char* argv[])
{
  const std::vector<std::string> chr_perm = {"I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX", "X", "XI", "XII", "XIII", "XIV", "XV", "XVI"};

  if (argc != 3)
  {
    std::cout << "Use: " << argv[0] << " [SRRXX_X_X_pileup.dat]  [SRRXX_X_X_mpileup.dat]" << std::endl; 
    return 1;
  }
  
  std::ifstream gz_in_buf(std::string{argv[1]}, std::ios_base::in | std::ios_base::binary);
  if (!gz_in_buf.is_open())
  {
    std::cerr << "File .dat not opened. Exit" << std::endl;
    return -1;
  }

  std::ofstream gz_out_buf(std::string{argv[2]}, std::ios_base::out | std::ios_base::binary);
  if (!gz_out_buf.is_open())
  { 
    std::cerr << "File .dat not opened. Exit" << std::endl;
    return -1;
  }
 
  std::string r_line;
  std::string chr_s;
  std::string base_num_s;
  std::string ref_s;
  std::string cov_s;
  std::string count_s;
  std::string mq_s;

  size_t count_lines = 0;
  ///////////////////////////
  //std::string delim = "\t";

  auto start_timer = std::chrono::high_resolution_clock::now();

  while (std::getline(gz_in_buf, r_line))
  {
    count_lines += 1;
    std::stringstream data_pileup(r_line);
    data_pileup >> chr_s >> base_num_s >> ref_s >> cov_s >> count_s >> mq_s;
    //
    //chr number selection
    //
    std::string chr_sel = get_right_of_delim(chr_s, "r");
    int chr_num = 0;
    if (std::find(chr_perm.begin(), chr_perm.end(), chr_sel) != chr_perm.end())
    {
      if (chr_sel == "I") {chr_num = 1;}
      if (chr_sel == "II") {chr_num = 2;}
      if (chr_sel == "III") {chr_num = 3;}
      if (chr_sel == "IV") {chr_num = 4;}
      if (chr_sel == "V") {chr_num = 5;}
      if (chr_sel == "VI") {chr_num = 6;}
      if (chr_sel == "VII") {chr_num = 7;}
      if (chr_sel == "VIII") {chr_num = 8;}
      if (chr_sel == "IX") {chr_num = 9;}
      if (chr_sel == "X") {chr_num = 10;}
      if (chr_sel == "XI") {chr_num = 11;}
      if (chr_sel == "XII") {chr_num = 12;}
      if (chr_sel == "XIII") {chr_num = 13;}
      if (chr_sel == "XIV") {chr_num = 14;}
      if (chr_sel == "XV") {chr_num = 15;}
      if (chr_sel == "XVI") {chr_num = 16;}
    }
    //
    if ((std::stoi(cov_s) != 0) && (ref_s == "A" || ref_s == "G"|| ref_s == "C" || ref_s == "T" || ref_s == "a" || ref_s == "g"|| ref_s == "c" || ref_s == "t") && (chr_num != 0)) //check if coverage is zero
    { 
      //parsing the last string: read base (. -> forward | , -> backward)  and the strand (forward AGCT | backward agct)

      //forward
      size_t A_count = 0; 
      size_t G_count = 0; 
      size_t C_count = 0; 
      size_t T_count = 0; 
      //backward
      size_t a_count = 0; 
      size_t g_count = 0; 
      size_t c_count = 0; 
      size_t t_count = 0; 

      size_t dot = 0;
      size_t comma = 0;
      //iteration over reads substring
      for (std::string::iterator i = count_s.begin(); i != count_s.end(); ++i)
      { 
        if (*i == '.' && *(i-1) != '^') dot++;
        if (*i == ',' && *(i-1) != '^') comma++;
        if (*i == 'A' && *(i-1) != '^') A_count++;
        if (*i == 'a' && *(i-1) != '^') a_count++;
        if (*i == 'G' && *(i-1) != '^') G_count++;
        if (*i == 'g' && *(i-1) != '^') g_count++;
        if (*i == 'C' && *(i-1) != '^') C_count++;
        if (*i == 'c' && *(i-1) != '^') c_count++;
        if (*i == 'T' && *(i-1) != '^') T_count++;
        if (*i == 't' && *(i-1) != '^') t_count++;
        if ((*i == '+' || *i == '-') && *(i-1) != '^') 
        {
          size_t first = std::distance(count_s.begin(), i+1);  
          size_t last = count_s.find_first_not_of("0123456789", first);
          std::advance(i, std::stoi(count_s.substr(first, last - first)) +
                          count_s.substr(first, last - first).length());
        }
      }

      //count reads in forward and reverse strains
      if (ref_s == "A" || ref_s == "a") { A_count += dot; a_count += comma; }
      if (ref_s == "G" || ref_s == "g") { G_count += dot; g_count += comma; } 
      if (ref_s == "C" || ref_s == "c") { C_count += dot; c_count += comma; }
      if (ref_s == "T" || ref_s == "t") { T_count += dot; t_count += comma; }

      //here we transform reference in uppercase nucleotides if necessary  
      if (ref_s == "a") { ref_s = "A"; }
      if (ref_s == "g") { ref_s = "G"; }
      if (ref_s == "c") { ref_s = "C"; }
      if (ref_s == "t") { ref_s = "T"; }

      //print on .dat file
      gz_out_buf << chr_num << "\t" 
          << base_num_s << "\t" <<  ref_s << "\t"
          << A_count << "\t" << G_count << "\t" << C_count << "\t" << T_count << "\t"
          << a_count << "\t" << g_count << "\t" << c_count << "\t" << t_count << "\n";
    }

    if (count_lines % 10000000 == 0)
    {
      std::cout << count_lines << "\n";
    }
  }

  auto stop_timer =  std::chrono::high_resolution_clock::now();
  auto measure_time =  std::chrono::duration_cast<std::chrono::minutes>(stop_timer - start_timer);

  std::cout << "Time taken to transform input file .gz in output file .dat:  " << measure_time.count() << "  min." << std::endl;

  gz_in_buf.close();

  gz_out_buf.close(); 

  return 0;
}
