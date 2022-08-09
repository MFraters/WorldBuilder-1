/*
  Copyright (C) 2021 by the authors of the World Builder code.

  This file is part of the World Builder.

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published
   by the Free Software Foundation, either version 2 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#define DOCTEST_CONFIG_SUPER_FAST_ASSERTS

#include "doctest/doctest.h"

#include "world_builder/objects/contours.h"

using namespace WorldBuilder;
using doctest::Approx;
using doctest::Contains;
/*
TEST_CASE("Contours")
{
  std::cout << "Countours" << std::endl;
  std::vector<std::pair<double,std::vector<double>>> contours_1;
  contours_1.emplace_back(0., std::vector<double> {1,5,3,5,4,6,5,5.5,6.5,4.5,8,6});
  contours_1.emplace_back(10., std::vector<double> {4.5,3.5,5.5,3.5});
  contours_1.emplace_back(10., std::vector<double> {1,1,3,1,5,1,7.5,1,8.5,1});
  Objects::Contours a(contours_1);

  for (auto &i : a.connectivity)
    {
      for (unsigned int j = 0; j < i.size(); j++)
        {
          std::cout << j << ": ";
          for (unsigned int k : i[j])
            {
              std::cout << " " << k;
            }
          std::cout << std::endl;
        }
    }

}*/