#ifndef basic_exception_H
#define basic_exception_H

#include "basic_constants.hpp"
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

class error_message : public std::exception {
public:
  std::string filename;
  std::string name;
  int line;
  std::string message;

  error_message(const std::string &filename_, const std::string &name_, const int line_, const std::string &message_ = std::string()) : filename(filename_), name(name_), line(line_), message(message_) { std::cout << nani() << std::endl; };

  std::string nani() const {
    std::ostringstream ss;
    ss << colorReset << "----------------------------------------------------\n" << colorReset << "    FILE: " << red << filename << ":" << Green << line << "\n" << colorReset << "FUNCTION: " << Magenta << name << "\n" << colorReset << " message: " << Red << message << colorReset << "\n" << colorReset << "----------------------------------------------------\n";
    return ss.str();
  };

  void print() const { std::cout << nani() << std::endl; };
};

inline std::string message(const std::string &filename_, const std::string &name_, const int line_, const std::string &message_) {
  std::ostringstream ss;
  ss << colorReset << Green << " LINE: " << line_ << colorReset << Magenta << " FUNCTION: " << name_ << colorReset << blue << " FILE: " << filename_ << colorReset << red << " : " << message_ << colorReset;
  return ss.str();
};

#endif
