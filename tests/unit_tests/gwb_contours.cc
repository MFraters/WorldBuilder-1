/*
  Copyright (C) 2023 by the authors of the World Builder code.

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

namespace WorldBuilder
{
  namespace Features
  {
    namespace FaultModels
    {
      namespace Composition
      {
        class Interface;
      }  // namespace Composition
      namespace Grains
      {
        class Interface;
      }  // namespace Grains
      namespace Temperature
      {
        class Interface;
      }  // namespace Temperature
    }  // namespace FaultModels
    namespace SubductingPlateModels
    {
      namespace Composition
      {
        class Interface;
      }  // namespace Composition
      namespace Grains
      {
        class Interface;
      }  // namespace Grains
      namespace Temperature
      {
        class Interface;
      }  // namespace Temperature
    }  // namespace SubductingPlateModels
  }  // namespace Features
} // namespace WorldBuilder
TEST_CASE("Contours")
{
  std::cout << "Countours" << std::endl;
  std::unique_ptr<CoordinateSystems::Interface> cartesian_system = CoordinateSystems::Interface::create("cartesian", nullptr);;
  std::cout << "ut flag 1" << std::endl;
  CoordinateSystem c = CoordinateSystem::cartesian;
  std::cout << "c=" << c << std::endl;
  std::vector<std::vector<Point<2> > > points = {{Point<2>(1,5,c),Point<2>(3,5,c),Point<2>(4,6,c),Point<2>(5,5.5,c),Point<2>(6.5,4.5,c),Point<2>(8,6,c)},
    {Point<2>(4.5,3.5,c),Point<2>(5.5,3.5,c)},
    {Point<2>(1,1,c),Point<2>(3,2,c),Point<2>(5,1,c),Point<2>(7.5,1,c),Point<2>(8.5,1,c)}
  };
  std::cout << "ut flag 2" << std::endl;
  std::vector<double> depths = {0.,1,3};
  std::vector<std::vector<double> > angle_contraints;
  std::vector<std::vector<double> > thicknesses = {{5.,5.,5.,5.,5.,5.},{2.5,2.5},{10.,10.,10.,10.,10.,10.}};
  std::cout << "ut flag 2.5" << std::endl;
  Objects::Contours<Features::SubductingPlateModels::Temperature::Interface,Features::SubductingPlateModels::Composition::Interface,Features::SubductingPlateModels::Grains::Interface> a(points, depths,thicknesses, 60., {}, {}, {}, {}, {});

  std::cout << "ut flag 3" << std::endl;
  const Point<3> position(4.75,3.,60-3.2,c);
  Objects::NaturalCoordinate natural_coordinate = Objects::NaturalCoordinate(position,*cartesian_system);

  std::cout << "ut flag 4" << std::endl;
  const Point<2> reference_point(4.5,-10.,c);

  const std::vector<std::vector<double> > interpolation_properties;

  std::cout << "ut flag 5" << std::endl;
  Objects::DistanceInterpolationData result = a.distance_interpolation_data(position,natural_coordinate,reference_point,cartesian_system,interpolation_properties,60.,false);

  std::cout << "signed_distance_from_bezier_surface:" << result.signed_distance_from_bezier_surface << ", distance_along_surface = " << result.distance_along_surface << std::endl;
}
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