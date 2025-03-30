#ifndef LAD_FILTERED_STREAM_BUF_H
#define LAD_FILTERED_STREAM_BUF_H

#include <iostream>
#include <sstream>
#include <streambuf>
#include <string>
#include <unordered_map>
#include <vector>

class LADFilteredStreamBuf : public std::streambuf {
public:
  LADFilteredStreamBuf(std::ostream &originalStream, const std::vector<std::string> &filterStrings);
  LADFilteredStreamBuf(std::ostream &originalStream);
  ~LADFilteredStreamBuf();

protected:
  int sync() override;
  int overflow(int c) override;

public:
  void addFilterString(const std::string &filterString);
  const std::unordered_map<std::string, int> &getFilterCounts() const; // New method to get filter counts
  void printFilterCounts() const {
    std::cout << "The following strings were filtered out so many times:" << std::endl;
    for (const auto &pair : filterCounts) {
      // Print the filtered strings and their counts
      std::cout << "Filtered string: " << pair.first << ", Count: " << pair.second << std::endl;
    }
  }

private:
  std::ostream &originalStream;                      // Reference to the original stream
  std::streambuf *originalBuffer;                    // Pointer to the original stream buffer
  std::string buffer;                                // Buffer to store the current line
  std::vector<std::string> filterStrings;            // Strings to filter out
  std::unordered_map<std::string, int> filterCounts; // Map to count filtered strings

  void flushBuffer();
};

#endif // LAD_FILTERED_STREAM_BUF_H
