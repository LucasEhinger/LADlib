#include "LADFilteredStreamBuf.h"
#include <Rtypes.h> // For Int_t and Bool_t
#include <TError.h> // For SetErrorHandler and gError
#include <TError.h> // Ensure gError is defined
#include <TROOT.h>  // Ensure ROOT headers are included properly
#include <algorithm>
// using namespace RO
// // Removed invalid namespace usage

LADFilteredStreamBuf *LADFilteredStreamBuf::currentInstance = nullptr; // Define static member

LADFilteredStreamBuf::LADFilteredStreamBuf(std::ostream &originalStream, const std::vector<std::string> &filterStrings)
    : originalStream(originalStream), filterStrings(filterStrings), originalBuffer(originalStream.rdbuf()) {
  currentInstance = this; // Set current instance
  // Initialize filter counts
  for (const auto &filter : filterStrings) {
    filterCounts[filter] = 0;
  }

  // Redirect the stream to this custom buffer
  originalStream.rdbuf(this);

  // SetErrorHandler(&LADFilteredStreamBuf::StaticCustomErrorHandler); // Use a static member function
}

LADFilteredStreamBuf::LADFilteredStreamBuf(std::ostream &originalStream)
    : originalStream(originalStream), originalBuffer(originalStream.rdbuf()) {
  // Redirect the stream to this custom buffer
  currentInstance = this; // Set current instance
  originalStream.rdbuf(this);

  // SetErrorHandler(&LADFilteredStreamBuf::StaticCustomErrorHandler); // Use a static member function
}

LADFilteredStreamBuf::~LADFilteredStreamBuf() {
  // Restore the original buffer when the object is destroyed
  originalStream.rdbuf(originalBuffer);
  // Print the filter counts
  printFilterCounts();
  currentInstance = nullptr; // Reset current instance
}

void LADFilteredStreamBuf::StaticCustomErrorHandler(int level, bool abort, const char *location, const char *message) {
  std::string msg(message);
  std::string loc(location);
  if (currentInstance &&
      currentInstance->getincludeRootErrorMessages()) { // Ensure includeRootErrorMessages is a member variable
    //     // If the message is not included, return early

    // Retrieve the current instance of LADFilteredStreamBuf
    if (LADFilteredStreamBuf::currentInstance) {
      for (const auto &filter : currentInstance->filterStrings) {
        if (msg.find(filter) != std::string::npos || loc.find(filter) != std::string::npos ) {
          // Increment the count for the matched filter
          currentInstance->filterCounts[filter]++;
          // Suppress messages containing any of the filter strings
          return;
        }
      }
    }
  }

  // Otherwise, print the message to the terminal
  fprintf(stderr, "ROOT [%d] %s: %s\n", level, location, message);

  if (abort) {
    std::abort();
  }
}

int LADFilteredStreamBuf::sync() {
  flushBuffer();
  return 0;
}

int LADFilteredStreamBuf::overflow(int c) {
  if (c == EOF) {
    flushBuffer();
    return c;
  }

  // Add the character to the buffer
  buffer += static_cast<char>(c);

  // If a newline is encountered, process the buffer
  if (c == '\n' || c == '\r') {
    flushBuffer();
  }

  return c;
}

void LADFilteredStreamBuf::flushBuffer() {
  if (!buffer.empty()) {
    // Check if the buffer contains any of the filter strings
    bool shouldFilter = false;
    for (const auto &filter : filterStrings) {
      if (buffer.find(filter) != std::string::npos) {
        shouldFilter = true;
        filterCounts[filter]++; // Increment the count for the matched filter
        break;
      }
    }

    // If no filter string is found, write the buffer to the original stream
    if (!shouldFilter) {
      originalBuffer->sputn(buffer.c_str(), buffer.size());
      originalBuffer->pubsync(); // Ensure the buffer is flushed to the original stream
    }

    // Clear the buffer
    buffer.clear();
  }
}

void LADFilteredStreamBuf::addFilterString(const std::string &filterString) {
  filterStrings.push_back(filterString);
  filterCounts[filterString] = 0; // Initialize the count for the new filter string
}

const std::unordered_map<std::string, int> &LADFilteredStreamBuf::getFilterCounts() const { return filterCounts; }
