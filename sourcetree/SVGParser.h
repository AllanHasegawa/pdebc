/*
 Copyright 2012 Allan Yoshio Hasegawa

 Licensed under the Apache License, Version 2.0 (the "License");
 you may not use this file except in compliance with the License.
 You may obtain a copy of the License at

 http://www.apache.org/licenses/LICENSE-2.0

 Unless required by applicable law or agreed to in writing, software
 distributed under the License is distributed on an "AS IS" BASIS,
 WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 See the License for the specific language governing permissions and
 limitations under the License.
 -----------------------------------------------------------------------------
 */

#ifndef SVGPARSER_H_
#define SVGPARSER_H_

#include <exception>
#include <vector>
#include <array>
#include "Vec2.h"

class ErrorParsingSVG : public std::exception {

};

class SVGParser {
 public:
  static void GetDataPoints(const char* file_name,
                            std::vector<Vec2>& data_points);

 private:
  SVGParser();
  virtual ~SVGParser();

  static std::vector<std::string> StringSplit(const std::string& input,
                                              const std::string& str_regex);
};

#endif /* SVGPARSER_H_ */
