/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "oops/util/Duration.h"

#include<stdint.h>

#include <istream>
#include <locale>
#include <sstream>
#include <string>

namespace util {

// -----------------------------------------------------------------------------

Duration::Duration() : seconds_(0) {}

// -----------------------------------------------------------------------------

Duration::~Duration() {}

// -----------------------------------------------------------------------------

Duration::Duration(const int64_t s) : seconds_(s) {}

// -----------------------------------------------------------------------------

Duration::Duration(const std::string & s) {
  this->set(s);
}

// -----------------------------------------------------------------------------

void Duration::set(const std::string & str) {
  this ->seconds_ = 0;

  std::istringstream is(str);

  char c = is.get();
  if (!is.good()) {failBadFormat(str);}

  // strip off any initial '-' and/or 'P'

  bool negative = false;
  if (c == '-') {
    negative = true;
    c = is.get();
    if (!is.good()) {failBadFormat(str);}
  }

  if (c != 'P') {failBadFormat(str);}

  // We now expect either 'T', or an integer followed by 'D'

  c = is.get();
  if (!is.good()) {failBadFormat(str);}

  if (c != 'T') {
    is.unget();
    std::string days = eatDigits(is);
    if (days.length() == 0) {failBadFormat(str);}
    c = is.get();
    if (!is.good() || c != 'D') {failBadFormat(str);}
    std::istringstream mys(days);
    int64_t ndays;
    mys >> ndays;
    if (mys.fail()) {failBadFormat(str);}
    this ->seconds_ = 86400*ndays;
  } else {
    is.unget();
  }

  // Now, we can have the end of the string, or a 'T'

  c = is.get();
  if (!is.eof()) {
    if (c != 'T') {failBadFormat(str);}

    // Now we can have an integer followed by 'H', 'M' or 'S'

    bool hoursdone = false;
    bool minsdone = false;
    bool secsdone = false;

    for (int i = 0; i < 3; i++) {
      c = is.peek();
      if (is.eof()) break;
      std::string num = eatDigits(is);
      c = is.get();
      if (!is.good()) {failBadFormat(str);}

      if (c == 'H' && !hoursdone && !minsdone && !secsdone) {
        std::istringstream mys(num);
        int64_t hours;
        mys >> hours;
        if (mys.fail()) {failBadFormat(str);}
        this ->seconds_ += 3600*hours;
        hoursdone = true;
      } else if (c == 'M' && !minsdone && !secsdone) {
        std::istringstream mys(num);
        int64_t mins;
        mys >> mins;
        if (mys.fail()) {failBadFormat(str);}
        this ->seconds_ += 60*mins;
        minsdone = true;
      } else if (c == 'S' && !secsdone) {
        std::istringstream mys(num);
        int64_t secs;
        mys >> secs;
        if (mys.fail()) {failBadFormat(str);}
        this ->seconds_ += secs;
        secsdone = true;
      } else {
        failBadFormat(str);
      }
    }
  }

  c = is.peek();
  if (!is.eof()) {failBadFormat(str);}

  if (negative) this->negate();
}

// -----------------------------------------------------------------------------

void Duration::failBadFormat(const std::string& str) const {
  std::string message = "Badly formed duration string: ";
  message.append(str);
  throw eckit::BadValue(message);
}

// -----------------------------------------------------------------------------

std::string Duration::eatDigits(std::istream & is) {
  std::string str;
  char c;
  std::locale loc;
  while (isdigit(c = is.get(), loc)) {
    if (!is.good()) throw eckit::BadCast("Error converting string to Duration");
    str.push_back(c);
  }
  is.unget();
  return str;
}

// -----------------------------------------------------------------------------

int64_t Duration::toSeconds () const {return seconds_;}

// -----------------------------------------------------------------------------

std::string Duration::toString() const {
  std::ostringstream os;
  int64_t remainder = seconds_;

  if (remainder < 0) {
    os << "-P";
    remainder = -remainder;
  } else {
    os << "P";
  }

  int64_t days = remainder/86400;
  if (days != 0) {
    os << days << "D";
    remainder = remainder - 86400*days;
  }

  if (remainder != 0 || days == 0) {
    os << "T";
  }

  int64_t hours = remainder/3600;
  if (hours != 0) {
    os << hours << "H";
    remainder = remainder - 3600*hours;
  }

  int64_t minutes = remainder/60;
  if (minutes != 0) {
    os << minutes << "M";
    remainder = remainder - 60*minutes;
  }

  if (remainder != 0 || seconds_ == 0) {
    os << remainder << "S";
  }

  return os.str();
}

// -----------------------------------------------------------------------------

bool Duration::operator==(const Duration& other) const {
  return this->seconds_ == other.seconds_;
}

// -----------------------------------------------------------------------------

bool Duration::operator!=(const Duration& other) const {
  return this->seconds_ != other.seconds_;
}

// -----------------------------------------------------------------------------

bool Duration::operator<(const Duration& other) const {
  return this->seconds_ < other.seconds_;
}

// -----------------------------------------------------------------------------

bool Duration::operator<=(const Duration& other) const {
  return this->seconds_ <= other.seconds_;
}

// -----------------------------------------------------------------------------

bool Duration::operator>(const Duration& other) const {
  return this->seconds_ > other.seconds_;
}

// -----------------------------------------------------------------------------

bool Duration::operator>=(const Duration& other) const {
  return this->seconds_ >= other.seconds_;
}

// -----------------------------------------------------------------------------

int Duration::operator%(const Duration& other) const {
  return this->seconds_ % other.seconds_;
}

// -----------------------------------------------------------------------------

void Duration::operator+=(const Duration & other) {
  seconds_ += other.seconds_;
}

// -----------------------------------------------------------------------------

void Duration::operator-=(const Duration & other) {
  seconds_ -= other.seconds_;
}

// -----------------------------------------------------------------------------

Duration operator+ (const Duration & d1, const Duration & d2) {
  Duration sum(d1);
  sum += d2;
  return sum;
}

// -----------------------------------------------------------------------------

Duration operator- (const Duration & d1, const Duration & d2) {
  Duration sum(d1);
  sum -= d2;
  return sum;
}

}  // namespace util
