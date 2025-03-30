#include "LADFilteredStreamBuf.h"
#include <algorithm>

LADFilteredStreamBuf::LADFilteredStreamBuf(std::ostream &originalStream, const std::vector<std::string> &filterStrings)
    : originalStream(originalStream), filterStrings(filterStrings), originalBuffer(originalStream.rdbuf()) {
    // Initialize filter counts
    for (const auto &filter : filterStrings) {
        filterCounts[filter] = 0;
    }

    // Redirect the stream to this custom buffer
    originalStream.rdbuf(this);
}

LADFilteredStreamBuf::LADFilteredStreamBuf(std::ostream &originalStream)
    : originalStream(originalStream), originalBuffer(originalStream.rdbuf()) {
    // Redirect the stream to this custom buffer
    originalStream.rdbuf(this);
}

LADFilteredStreamBuf::~LADFilteredStreamBuf() {
    // Restore the original buffer when the object is destroyed
    originalStream.rdbuf(originalBuffer);
    // Print the filter counts
    printFilterCounts();
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
    if (c == '\n') {
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
        }

        // Clear the buffer
        buffer.clear();
    }
}

void LADFilteredStreamBuf::addFilterString(const std::string &filterString) {
    filterStrings.push_back(filterString);
    filterCounts[filterString] = 0; // Initialize the count for the new filter string
}

const std::unordered_map<std::string, int>& LADFilteredStreamBuf::getFilterCounts() const {
    return filterCounts;
}
