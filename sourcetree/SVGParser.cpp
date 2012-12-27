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

#include "/home/aran/Projects/libs/pdebc/sourcetree/SVGParser.h"

#include <vector>
#include <array>

#include "external/tinyxml2/tinyxml2.h"

SVGParser::SVGParser() {
	// TODO Auto-generated constructor stub

}

SVGParser::~SVGParser() {
	// TODO Auto-generated destructor stub
}

void SVGParser::GetControlPoints(const char* file_name,
		std::vector<std::array<double, 2>>& control_points) {
	using namespace tinyxml2;


	XMLDocument doc;
	if (doc.LoadFile(file_name) != XML_SUCCESS) {
		throw ErrorParsingSVG();
	}

	control_points.clear();

	XMLElement* g_element = doc.FirstChildElement("svg")->FirstChildElement(
			"g");

	XMLElement* g_paths = g_element->FirstChildElement("path");

	while (g_paths != NULL) {
		const char* a_type = g_paths->Attribute("sodipodi:type");
		if (strcmp(a_type, "arc") == 0) {
			double cx;
			double cy;
			g_paths->QueryDoubleAttribute("cx", &cx);
			g_paths->QueryDoubleAttribute("cy", &cy);
			const char* a_cx = g_paths->Attribute("sodipodi:cx");
			const char* a_cy = g_paths->Attribute("sodipodi:cy");
			printf("Path: %s\n%s\n%s\n", a_type, a_cx, a_cy);
			std::array<double,2> cp;
			cp[0] = cx;
			cp[1] = cy;
			control_points.push_back(cp);
		}
		g_paths = g_paths->NextSiblingElement("path");
	}
}
