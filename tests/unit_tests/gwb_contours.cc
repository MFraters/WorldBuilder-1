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

#include "world_builder/coordinate_system.h"
#include "world_builder/objects/bezier_curve.h"
#include "world_builder/objects/bezier_surface_patch.h"
#include "world_builder/objects/bezier_surface_patches.h"
#include "world_builder/point.h"
#define DOCTEST_CONFIG_SUPER_FAST_ASSERTS

#include "doctest/doctest.h"

#include "world_builder/objects/contours.h"
#include "world_builder/utilities.h"

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
  std::vector<std::vector<double> > top_truncations = {{0.,0.,0.,0.,0.,0.},{0.0,0.0},{0.,0.,0.,0.,0.,0.}};
  std::cout << "ut flag 2.5" << std::endl;
  Objects::Contours<Features::SubductingPlateModels::Temperature::Interface,Features::SubductingPlateModels::Composition::Interface,Features::SubductingPlateModels::Grains::Interface> a(points, depths,thicknesses,top_truncations, 60., {}, {}, {}, {}, {});

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

TEST_CASE("Contours: Connectivity")
{
  std::vector<Objects::BezierCurve> curves =
  {
    Objects::BezierCurve(
    {
      Point<2>(1.,5.,cartesian),
      Point<2>(3.,5.,cartesian),
      Point<2>(4.,6.,cartesian),
      Point<2>(5.,5.5,cartesian),
      Point<2>(6.5,4.5,cartesian),
      Point<2>(8.,5.5,cartesian),
    }),
    Objects::BezierCurve(
    {
      Point<2>(4.5,3.75,cartesian),
      Point<2>(5.5,3.75,cartesian),
    }),
    Objects::BezierCurve(
    {
      Point<2>(1.,2.,cartesian),
      Point<2>(3.,2.,cartesian),
      Point<2>(4.5,1.,cartesian),
      Point<2>(5.5,1.5,cartesian),
      Point<2>(6.5,3.5,cartesian),
      Point<2>(8,2.,cartesian),
    }),
  };
  const std::vector<std::vector<std::vector<unsigned int> > > computed_connectivity = Utilities::create_contour_connectivity(curves);
  const std::vector<std::vector<std::vector<unsigned int> > > reference_connectivity = {{{0},{0},{1},{1},{1},{1}},{{0, 1, 2},{3, 4, 5}}};

  REQUIRE(computed_connectivity.size() == reference_connectivity.size());
  for (size_t i = 0; i < computed_connectivity.size(); ++i)
    {
      REQUIRE(computed_connectivity[i].size() == reference_connectivity[i].size());
      for (size_t j = 0; j < computed_connectivity[i].size(); ++j)
        {
          REQUIRE(computed_connectivity[i][j].size() == reference_connectivity[i][j].size());
          for (size_t k = 0; k < computed_connectivity[i][j].size(); ++k)
            {
              CHECK(computed_connectivity[i][j][k] == reference_connectivity[i][j][k]);
            }
        }
    }
}

TEST_CASE("Bezier Curve: approximate length")
{
  using doctest::Approx;
  Objects::BezierCurve curve = Objects::BezierCurve({Point<2>(0.0,0.0,cartesian),Point<2>(0.0,1.0,cartesian)});
  CHECK(1.0 == Approx(curve.approximate_length(1)));
  CHECK(1.0 == Approx(curve.approximate_length(2)));
  CHECK(1.0 == Approx(curve.approximate_length(3)));
  CHECK(1.0 == Approx(curve.approximate_length(10)));


  curve = Objects::BezierCurve({Point<2>(0.0,0.0,cartesian),Point<2>(1.0,0.0,cartesian)});
  CHECK(1.0 == Approx(curve.approximate_length(1)));
  CHECK(1.0 == Approx(curve.approximate_length(2)));
  CHECK(1.0 == Approx(curve.approximate_length(3)));
  CHECK(1.0 == Approx(curve.approximate_length(10)));

  curve = Objects::BezierCurve({Point<2>(0.0,0.0,cartesian),Point<2>(1.0,1.0,cartesian)});
  CHECK(std::sqrt(2) == Approx(curve.approximate_length(1)));
  CHECK(std::sqrt(2) == Approx(curve.approximate_length(2)));
  CHECK(std::sqrt(2) == Approx(curve.approximate_length(3)));
  CHECK(std::sqrt(2) == Approx(curve.approximate_length(4)));
  CHECK(std::sqrt(2) == Approx(curve.approximate_length(5)));
  CHECK(std::sqrt(2) == Approx(curve.approximate_length(6)));
  CHECK(std::sqrt(2) == Approx(curve.approximate_length(7)));
  CHECK(std::sqrt(2) == Approx(curve.approximate_length(8)));
  CHECK(std::sqrt(2) == Approx(curve.approximate_length(9)));
  CHECK(std::sqrt(2) == Approx(curve.approximate_length(10)));

  curve = Objects::BezierCurve({Point<2>(0.0,0.0,cartesian),
                                Point<2>(0.0,1.0,cartesian),
                                Point<2>(0.0,2.0,cartesian)
                               });

  size_t counter = 0;
  for (auto control_point_i : curve.get_control_points())
    {
      counter++;
    }
  CHECK(2. == Approx(curve.approximate_length(1)));
  CHECK(2. == Approx(curve.approximate_length(2)));
  CHECK(2. == Approx(curve.approximate_length(3)));
  CHECK(2. == Approx(curve.approximate_length(4)));
  CHECK(2. == Approx(curve.approximate_length(5)));
  CHECK(2. == Approx(curve.approximate_length(6)));
  CHECK(2. == Approx(curve.approximate_length(7)));
  CHECK(2. == Approx(curve.approximate_length(8)));
  CHECK(2. == Approx(curve.approximate_length(9)));
  CHECK(2. == Approx(curve.approximate_length(10)));


  curve = Objects::BezierCurve({Point<2>(0.0,0.0,cartesian),
                                Point<2>(1.0,1.0,cartesian),
                                Point<2>(2.0,2.0,cartesian)
                               });
  CHECK(2.*std::sqrt(2) == Approx(curve.approximate_length(1)));
  CHECK(2.*std::sqrt(2) == Approx(curve.approximate_length(2)));
  CHECK(2.*std::sqrt(2) == Approx(curve.approximate_length(3)));
  CHECK(2.*std::sqrt(2) == Approx(curve.approximate_length(4)));
  CHECK(2.*std::sqrt(2) == Approx(curve.approximate_length(5)));
  CHECK(2.*std::sqrt(2) == Approx(curve.approximate_length(6)));
  CHECK(2.*std::sqrt(2) == Approx(curve.approximate_length(7)));
  CHECK(2.*std::sqrt(2) == Approx(curve.approximate_length(8)));
  CHECK(2.*std::sqrt(2) == Approx(curve.approximate_length(9)));
  CHECK(2.*std::sqrt(2) == Approx(curve.approximate_length(10)));


  curve = Objects::BezierCurve({Point<2>(0.0,0.0,cartesian),
                                Point<2>(1.0,1.0,cartesian),
                                Point<2>(2.0,3.0,cartesian)
                               });
  CHECK(3.6502815399 == Approx(curve.approximate_length(1)));
  CHECK(3.6514377005 == Approx(curve.approximate_length(2)));
  CHECK(3.8413917324 == Approx(curve.approximate_length(10)));
  CHECK(3.8469724468 == Approx(curve.approximate_length(100)));
  CHECK(3.847156369 == Approx(curve.approximate_length(200)));
  CHECK(3.847156369 == Approx(curve.approximate_length(1000)));
  CHECK(3.847156369 == Approx(curve.approximate_length(2000)));
  CHECK(3.847156369 == Approx(curve.approximate_length(10000)));
  CHECK(3.847156369 == Approx(curve.approximate_length(20000)));
}



TEST_CASE("Bezier Curve: length map")
{
  using doctest::Approx;
  {
    Objects::BezierCurve curve = Objects::BezierCurve({Point<2>(0.0,0.0,cartesian),Point<2>(0.0,10.0,cartesian)});

    const std::vector<std::vector<double> > parameter_to_length_map = curve.compute_parameter_to_length_map(5);
    const double arc_length = parameter_to_length_map[parameter_to_length_map.size()-1][parameter_to_length_map[0].size()-1];
    CHECK(10. == Approx(arc_length));
    CHECK(parameter_to_length_map.size() == 1);
    CHECK(parameter_to_length_map[0].size() == 6);
    CHECK(parameter_to_length_map[0][0] == Approx(0.00)); // 0
    CHECK(parameter_to_length_map[0][1] == Approx(1.04)); // 0.2
    CHECK(parameter_to_length_map[0][2] == Approx(3.52)); // 0.4
    CHECK(parameter_to_length_map[0][3] == Approx(6.48)); // 0.6
    CHECK(parameter_to_length_map[0][4] == Approx(8.96)); // 0.8
    CHECK(parameter_to_length_map[0][5] == Approx(10.0)); // 1.0

    CHECK(curve.distance_to_t(parameter_to_length_map,-1) == Approx(-0.1));
    CHECK(curve.distance_to_t(parameter_to_length_map,0.0) == Approx(0.0));
    CHECK(curve.distance_to_t(parameter_to_length_map,1) == Approx(0.1923076923));
    CHECK(curve.distance_to_t(parameter_to_length_map,1.04) == Approx(0.2));
    CHECK(curve.distance_to_t(parameter_to_length_map,2) == Approx(0.2774193548));
    CHECK(curve.distance_to_t(parameter_to_length_map,3.52) == Approx(0.4));
    CHECK(curve.distance_to_t(parameter_to_length_map,4) == Approx(0.4324324324));
    CHECK(curve.distance_to_t(parameter_to_length_map,5) == Approx(0.5));
    CHECK(curve.distance_to_t(parameter_to_length_map,6) == Approx(0.5675675676));
    CHECK(curve.distance_to_t(parameter_to_length_map,6.48) == Approx(0.6));
    CHECK(curve.distance_to_t(parameter_to_length_map,8) == Approx(0.7225806452));
    CHECK(curve.distance_to_t(parameter_to_length_map,8.96) == Approx(0.8));
    CHECK(curve.distance_to_t(parameter_to_length_map,9) == Approx(0.8076923077));
    CHECK(curve.distance_to_t(parameter_to_length_map,10.0) == Approx(1.0));
    CHECK(curve.distance_to_t(parameter_to_length_map,11) == Approx(1.1));
  }
  {
    Objects::BezierCurve curve = Objects::BezierCurve({Point<2>(0.0,0.0,cartesian),Point<2>(0.0,10.0,cartesian),Point<2>(0.0,20.0,cartesian)});

    const std::vector<std::vector<double> > parameter_to_length_map = curve.compute_parameter_to_length_map(5);
    const double arc_length = parameter_to_length_map[parameter_to_length_map.size()-1][parameter_to_length_map[0].size()-1];
    CHECK(20. == Approx(arc_length));
    CHECK(parameter_to_length_map.size() == 2);
    CHECK(parameter_to_length_map[0].size() == 6);
    CHECK(parameter_to_length_map[0][0] == Approx(0.000)); // 0
    CHECK(parameter_to_length_map[0][1] == Approx(1.616)); // 0.2
    CHECK(parameter_to_length_map[0][2] == Approx(3.808)); // 0.4
    CHECK(parameter_to_length_map[0][3] == Approx(6.192)); // 0.6
    CHECK(parameter_to_length_map[0][4] == Approx(8.384)); // 0.8
    CHECK(parameter_to_length_map[0][5] == Approx(10.0)); // 1.0
    CHECK(parameter_to_length_map[1][1] == Approx(11.616)); // 1.2
    CHECK(parameter_to_length_map[1][2] == Approx(13.808)); // 1.4
    CHECK(parameter_to_length_map[1][3] == Approx(16.192)); // 1.6
    CHECK(parameter_to_length_map[1][4] == Approx(18.384)); // 1.8
    CHECK(parameter_to_length_map[1][5] == Approx(20.0)); // 2.0

    CHECK(curve.distance_to_t(parameter_to_length_map,-1) == Approx(-0.05));
    CHECK(curve.distance_to_t(parameter_to_length_map,0.0) == Approx(0.0));
    CHECK(curve.distance_to_t(parameter_to_length_map,1) == Approx(0.1237623762));
    CHECK(curve.distance_to_t(parameter_to_length_map,1.616) == Approx(0.2));
    CHECK(curve.distance_to_t(parameter_to_length_map,2) == Approx(0.2350364964));
    CHECK(curve.distance_to_t(parameter_to_length_map,3.808) == Approx(0.4));
    CHECK(curve.distance_to_t(parameter_to_length_map,4) == Approx(0.4161073826));
    CHECK(curve.distance_to_t(parameter_to_length_map,5) == Approx(0.5));
    CHECK(curve.distance_to_t(parameter_to_length_map,6) == Approx(0.5838926174));
    CHECK(curve.distance_to_t(parameter_to_length_map,6.192) == Approx(0.6));
    CHECK(curve.distance_to_t(parameter_to_length_map,8) == Approx(0.7649635036));
    CHECK(curve.distance_to_t(parameter_to_length_map,8.384) == Approx(0.8));
    CHECK(curve.distance_to_t(parameter_to_length_map,9) == Approx(0.8762376238));
    CHECK(curve.distance_to_t(parameter_to_length_map,10.0) == Approx(1.0));

    CHECK(curve.distance_to_t(parameter_to_length_map,11) == Approx(1.1893939394));
    CHECK(curve.distance_to_t(parameter_to_length_map,11.616) == Approx(1.2));
    CHECK(curve.distance_to_t(parameter_to_length_map,12) == Approx(1.2350364964));
    CHECK(curve.distance_to_t(parameter_to_length_map,13.808) == Approx(1.4));
    CHECK(curve.distance_to_t(parameter_to_length_map,14) == Approx(1.4161073826));
    CHECK(curve.distance_to_t(parameter_to_length_map,15) == Approx(1.5));
    CHECK(curve.distance_to_t(parameter_to_length_map,16) == Approx(1.5838926174));
    CHECK(curve.distance_to_t(parameter_to_length_map,16.192) == Approx(1.6));
    CHECK(curve.distance_to_t(parameter_to_length_map,18) == Approx(1.7649635036));
    CHECK(curve.distance_to_t(parameter_to_length_map,18.384) == Approx(1.8));
    CHECK(curve.distance_to_t(parameter_to_length_map,19) == Approx(1.8762376238));
    CHECK(curve.distance_to_t(parameter_to_length_map,20.0) == Approx(2.0));
    CHECK(curve.distance_to_t(parameter_to_length_map,21) == Approx(2.05));
  }
}


TEST_CASE("Bezier Patch")
{
  using doctest::Approx;
  {
    auto control_points = std::array<std::array<Point<3>,4>,4> {{{{Point<3>(0.,0,0,cartesian),Point<3>(0.,1,0,cartesian),Point<3>(0.,2,0,cartesian),Point<3>(0.,3,0,cartesian)}},
        {Point<3>(0.,0,1,cartesian),Point<3>(0.,1,1.,cartesian),Point<3>(0.,2,1.,cartesian),Point<3>(0.,3,1.,cartesian)},
        {Point<3>(0.,0,2,cartesian),Point<3>(0.,1,2.,cartesian),Point<3>(0.,2,2.,cartesian),Point<3>(0.,3,2.,cartesian)},
        {Point<3>(0.,0,3,cartesian),Point<3>(0.,1,3.,cartesian),Point<3>(0.,2,3.,cartesian),Point<3>(0.,3,3.,cartesian)}}};
    
    auto patch = Objects::BezierSurfacePatch(control_points);

    CHECK(patch.point_on_surface(0, 0) == Point<3>(0,0,0,cartesian));
    CHECK(patch.point_on_surface(0.5, 0) == Point<3>(0.5,0,0,cartesian));


  }

}

TEST_CASE("Bezier Curve: create_patches_from_contours")
{
  using doctest::Approx;
  {
    std::vector<Objects::BezierCurve> curves =
    {
      Objects::BezierCurve({Point<2>(0.0,0.0,cartesian),Point<2>(0.0,10.0,cartesian),Point<2>(0.0,20.0,cartesian)}),
      Objects::BezierCurve({Point<2>(10.0,0.0,cartesian),Point<2>(10.0,10.0,cartesian),Point<2>(10.0,20.0,cartesian)})
    };
    std::vector<double> depths = {0,10};

    auto patches = Objects::BezierSurfacePatches(curves, depths, {{0,45,90.},{0,450,900}}, {}, {}, {});



  }

}