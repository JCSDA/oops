/*
 * (C) Copyright 2020 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "oops/util/parameters/PropertyJsonSchema.h"

#include <sstream>
#include <string>

#include "eckit/log/Channel.h"

namespace oops {

std::string toString(const PropertyJsonSchema &schema) {
  std::stringstream str;
  {
    eckit::Channel ch;
    ch.setStream(str);

    ch << '{';
    if (!schema.empty()) {
      ch << '\n';
      bool needsCommaAndNewline = false;
      {
        eckit::AutoIndent indent(ch);
        for (const auto &nameAndValue : schema) {
          if (needsCommaAndNewline)
            ch << ",\n";
          ch << '"' << nameAndValue.first << '"' << ": " << nameAndValue.second;
          needsCommaAndNewline = true;
        }
      }
      if (needsCommaAndNewline)
        ch << '\n';  // this was the last item, so comma is not needed, just a newline
    }
    ch << '}';
  }
  return str.str();
}

}  // namespace oops
