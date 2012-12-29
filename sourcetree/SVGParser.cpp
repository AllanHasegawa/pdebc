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

#include "SVGParser.h"

#include <regex>
#include <vector>
#include <array>
#include <string>

#include "external/tinyxml2/tinyxml2.h"
#include "Vec2.h"

SVGParser::SVGParser() {
  // TODO Auto-generated constructor stub

}

SVGParser::~SVGParser() {
  // TODO Auto-generated destructor stub
}

void SVGParser::GetDataPoints(const char* file_name,
                              std::vector<Vec2>& data_points) {
  using namespace tinyxml2;
  using namespace std;

  XMLDocument doc;
  if (doc.LoadFile(file_name) != XML_SUCCESS) {
    throw ErrorParsingSVG();
  }

  data_points.clear();

  XMLElement* g_element = doc.FirstChildElement("svg");

  XMLElement* g_paths = g_element->FirstChildElement("path");

  if (g_paths != nullptr) {
    const char* a_type = g_paths->Attribute("d");
    if (a_type != 0) {
      string stall(a_type);

      vector<string> tokens = StringSplit(stall, "\\s+");

      // Something is wrong if there is not enough lines
      if (tokens.size() < 4) {
        throw ErrorParsingSVG();
      }

      // first 3 tokens are useless
      // lets get the first of each line
      int last_i;
      for (int i = 3; i < tokens.size(); i += 3) {
        string s = tokens[i];
        auto m_tokens = StringSplit(s, ",");
        Vec2 cp;
        cp.x = atof(m_tokens[0].c_str());
        cp.y = atof(m_tokens[1].c_str());
        data_points.push_back(cp);
        last_i = i;
      }
      string s = tokens[last_i+1];
      auto m_tokens = StringSplit(s, ",");
      Vec2 cp;
      cp.x = atof(m_tokens[0].c_str());
      cp.y = atof(m_tokens[1].c_str());
      data_points.push_back(cp);
    }
  }
}

std::vector<std::string> SVGParser::StringSplit(const std::string& input,
                                                const std::string& str_regex) {
  using namespace std;
  regex ws_re(str_regex);
  vector<string> v;
  sregex_token_iterator i(input.begin(), input.end(), ws_re, -1);
  while (i != sregex_token_iterator()) {
    v.push_back(i->str());
    i++;
  }
  return v;
}
