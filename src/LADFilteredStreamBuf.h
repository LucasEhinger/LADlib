#ifndef LAD_FILTERED_STREAM_BUF_H
#define LAD_FILTERED_STREAM_BUF_H

#include <TError.h> // For SetErrorHandler and gError
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
  bool incRootErrorMessages = false; // Flag to include ROOT error messages

public:
  static LADFilteredStreamBuf *currentInstance; // Declare static member
  static void StaticCustomErrorHandler(int level, bool abort, const char *location, const char *message);
  void addFilterString(const std::string &filterString);
  void includeRootErrorMessages(bool include = true) {
    incRootErrorMessages = include;
  } // New method to include ROOT error messages
  const std::unordered_map<std::string, int> &getFilterCounts() const; // New method to get filter counts
  const bool getincludeRootErrorMessages() const { // Getter for includeRootErrorMessages
    return incRootErrorMessages;
  }
  const std::vector<std::string> &getFilterStrings() const {           // Getter for filter strings
    return filterStrings;
  }
  void incrementFilterCount(const std::string &filterString) { // Incrementer for counts
    ++filterCounts[filterString];
  }
  void printFilterCounts() const {
    std::cout << "--------------------------------------------" << std::endl;
    std::cout << "LADFilteredStreamBuf: Filtered strings count" << std::endl;
    for (const auto &pair : filterCounts) {
      std::cout << "Filtered string: \"" << pair.first << "\",    Count: " << pair.second << std::endl;
    }
    std::cout << "--------------------------------------------" << std::endl;
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
