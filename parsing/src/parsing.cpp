#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <iomanip>
#include <chrono>
#include <boost/iostreams/filtering_streambuf.hpp>
//#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/filter/gzip.hpp>
//#include </usr/include/boost/iostreams/filter/gzip.hpp>

namespace bi = boost::iostreams;

int main(int argc, char* argv[])
{
  constexpr size_t num_split = 7;

  if (argc != 3)
  {
    std::cout << "Use: " << argv[0] << " [CRCXXXX-XX-X-X.callable_CN-X_pseudompileup.gz]  [CRCXXXX-XX-X.[ancestor|endpoint]_CN-X_pseudopileup.dat]" << std::endl; 
    return 1;
  }
  
  std::ifstream gz_in_buf(std::string{argv[1]}, std::ios_base::in | std::ios_base::binary);
  if (!gz_in_buf.is_open())
  {
    std::cerr << "File .gz not opened. Exit" << std::endl;
    return -1;
  }

  std::string out_file = argv[2]; 
  std::ofstream gz_out_buf(out_file, std::ios_base::out | std::ios_base::binary);
  if (!gz_out_buf.is_open())
  { 
    std::cerr << "File .dat not opened. Exit" << std::endl;
    return -1;
  }

  //input buffer
  bi::filtering_streambuf<bi::input> in_buf;
  in_buf.push(bi::gzip_decompressor());
  in_buf.push(gz_in_buf); 
 
  std::istream instream_buf(&in_buf); 
  std::string r_line;
  size_t count_lines = 0;
  ///////////////////////////
  std::string delim = "\t";

  auto start_timer = std::chrono::high_resolution_clock::now();

  while (std::getline(instream_buf, r_line))
  {
    count_lines += 1;
    std::vector<std::string> v_string;
    v_string.reserve(num_split);
    size_t first = 0;
    while (first < r_line.size())
    {
      const auto second = r_line.find_first_of(delim, first);    
      if (first != second)
      { 
        v_string.emplace_back(r_line.substr(first, second - first));
      }

      if (second == std::string::npos) break;
      first = second + 1;
    }
    
    if (std::stoi(v_string[4]) != 0) //check if coverage is zero
    { 
      //v_string[0].substr(3,0) //chromosome, es: chr1-->1
      //v_string[2] //loci number 
      //v_string[3] //reference base
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
      for (std::string::iterator i = v_string[5].begin(); i != v_string[5].end(); ++i)
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
          size_t first = std::distance(v_string[5].begin(), i+1);  
          size_t last = v_string[5].find_first_not_of("0123456789", first);
          std::advance(i, std::stoi(v_string[5].substr(first, last - first)) +
                          v_string[5].substr(first, last - first).length());
        }
      }

      //count reads in forward and reverse strains
      if (v_string[3][0] == 'A') { A_count += dot; a_count += comma; }
      if (v_string[3][0] == 'G') { G_count += dot; g_count += comma; } 
      if (v_string[3][0] == 'C') { C_count += dot; c_count += comma; }
      if (v_string[3][0] == 'T') { T_count += dot; t_count += comma; }


      //print on .dat file
      gz_out_buf << (static_cast<size_t>(v_string[0][4]) != std::string::npos ? 
                     std::string{v_string[0][3]}+std::string{v_string[0][4]} : 
                     std::string{v_string[0][3]}) 
          << '\t' << v_string[2] << '\t' <<  v_string[3][0] << '\t'
          << A_count << '\t' << G_count << '\t' << C_count << '\t' << T_count << '\t'  
          << a_count << '\t' << g_count << '\t' << c_count << '\t' << t_count << '\n';
    }

    if (count_lines % 10000000 == 0)
    {
      std::cout << count_lines << '\n';
    }
  }

  auto stop_timer =  std::chrono::high_resolution_clock::now();
  auto measure_time =  std::chrono::duration_cast<std::chrono::minutes>(stop_timer - start_timer);

  std::cout << "Time taken to transform input file .gz in output file .dat:  " << measure_time.count() << "  min." << std::endl;

  gz_in_buf.close();

  gz_out_buf.close(); 

  return 0;
}
