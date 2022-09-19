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

#include "world_builder/utilities.h"
#include "world_builder/coordinate_system.h"
#include "world_builder/objects/bezier_curve.h"
#include "world_builder/objects/closest_point_on_curve.h"

using namespace WorldBuilder;
using doctest::Approx;
using doctest::Contains;

/**
 * Compare the given two std::vector<double> entries with an epsilon (using Catch::Approx)
 */
inline void compare_vectors_approx(
  const std::vector<double> &computed,
  const std::vector<double> &expected)
{
  CHECK(computed.size() == expected.size());
  for (unsigned int i=0; i< computed.size(); ++i)
    {
      INFO("vector index i=" << i << ": ");
      CHECK(computed[i] == Approx(expected[i]));
    }
}

TEST_CASE("Cubic roots solver")
{
  std::cout << "Cubic roots solver:" << std::endl;
  double a = 1.;
  double b = -3.;
  double c = -144.;
  double d = 432.;
  auto cubic_roots = Objects::BezierCurve::solve_cubic_equation_real(a,b,c,d);
  for (size_t i = 0; i < cubic_roots.size(); ++i)
    {
      CHECK(std::abs(a*cubic_roots[i]*cubic_roots[i]*cubic_roots[i]
                     +b*cubic_roots[i]*cubic_roots[i]
                     +c*cubic_roots[i] + d) < 1e-10);
    }
  compare_vectors_approx(cubic_roots, {12,3,-12});

  a = 2.;
  b = -3.;
  c = -144.;
  d = 432.;
  cubic_roots = Objects::BezierCurve::solve_cubic_equation_real(a,b,c,d);
  for (size_t i = 0; i < cubic_roots.size(); ++i)
    {
      CHECK(std::abs(a*cubic_roots[i]*cubic_roots[i]*cubic_roots[i]
                     +b*cubic_roots[i]*cubic_roots[i]
                     +c*cubic_roots[i] + d) < 1e-10);
    }
  compare_vectors_approx(cubic_roots, {7.30783,3.25969,-9.06752});

  // see https://www.youtube.com/watch?v=N-KXStupwsc at 3:18
  a = 1.;
  b = 0.;
  c = -15.;
  d = -126.;
  cubic_roots = Objects::BezierCurve::solve_cubic_equation_real(a,b,c,d);
  for (size_t i = 0; i < cubic_roots.size(); ++i)
    {
      CHECK(std::abs(a*cubic_roots[i]*cubic_roots[i]*cubic_roots[i]
                     +b*cubic_roots[i]*cubic_roots[i]
                     +c*cubic_roots[i] + d) < 1e-10);
    }
  compare_vectors_approx(cubic_roots, {6.});


  // see https://www.youtube.com/watch?v=N-KXStupwsc at 8:12
  a = 0.;
  b = 2.;
  c = 8.;
  d = -6.;
  cubic_roots = Objects::BezierCurve::solve_cubic_equation_real(a,b,c,d);

  for (size_t i = 0; i < cubic_roots.size(); ++i)
    {
      CHECK(std::abs(a*cubic_roots[i]*cubic_roots[i]*cubic_roots[i]
                     +b*cubic_roots[i]*cubic_roots[i]
                     +c*cubic_roots[i] + d) < 1e-10);
    }
  compare_vectors_approx(cubic_roots, {-2+sqrt(7),-2-sqrt(7)});

  a = 1.;
  b = 2.;
  c = 8.;
  d = 0.;
  cubic_roots = Objects::BezierCurve::solve_cubic_equation_real(a,b,c,d);

  for (size_t i = 0; i < cubic_roots.size(); ++i)
    {
      CHECK(std::abs(a*cubic_roots[i]*cubic_roots[i]*cubic_roots[i]
                     +b*cubic_roots[i]*cubic_roots[i]
                     +c*cubic_roots[i] + d) < 1e-10);
    }
  compare_vectors_approx(cubic_roots, {0.});

  // https://www.wolframalpha.com/widgets/view.jsp?id=578d50248844454e46e24e9ed230843d
  a = 1.;
  b = 5.;
  c = 2.;
  d = 0.;
  cubic_roots = Objects::BezierCurve::solve_cubic_equation_real(a,b,c,d);

  for (size_t i = 0; i < cubic_roots.size(); ++i)
    {
      CHECK(std::abs(a*cubic_roots[i]*cubic_roots[i]*cubic_roots[i]
                     +b*cubic_roots[i]*cubic_roots[i]
                     +c*cubic_roots[i] + d) < 1e-10);
    }
  compare_vectors_approx(cubic_roots, {0.,sqrt(17.)/2.-5./2.,-5./2.-sqrt(17)/2.});

  // see https://www.youtube.com/watch?v=N-KXStupwsc at 19:35
  a = 5.;
  b = 0.;
  c = 0.;
  d = 0.;
  cubic_roots = Objects::BezierCurve::solve_cubic_equation_real(a,b,c,d);

  for (size_t i = 0; i < cubic_roots.size(); ++i)
    {
      CHECK(std::abs(a*cubic_roots[i]*cubic_roots[i]*cubic_roots[i]
                     +b*cubic_roots[i]*cubic_roots[i]
                     +c*cubic_roots[i] + d) < 1e-10);
    }
  compare_vectors_approx(cubic_roots, {0.});

  // see https://www.youtube.com/watch?v=N-KXStupwsc at 24:22
  a = 1.;
  b = 0.;
  c = -6.;
  d = -40.;
  cubic_roots = Objects::BezierCurve::solve_cubic_equation_real(a,b,c,d);
  for (size_t i = 0; i < cubic_roots.size(); ++i)
    {
      CHECK(std::abs(a*cubic_roots[i]*cubic_roots[i]*cubic_roots[i]
                     +b*cubic_roots[i]*cubic_roots[i]
                     +c*cubic_roots[i] + d) < 1e-10);
    }
  compare_vectors_approx(cubic_roots, {4.});

  // see https://www.youtube.com/watch?v=N-KXStupwsc at 27:20
  a = 1.;
  b = 0.;
  c = -6.;
  d = -4.;
  cubic_roots = Objects::BezierCurve::solve_cubic_equation_real(a,b,c,d);
  for (size_t i = 0; i < cubic_roots.size(); ++i)
    {
      CHECK(std::abs(a*cubic_roots[i]*cubic_roots[i]*cubic_roots[i]
                     +b*cubic_roots[i]*cubic_roots[i]
                     +c*cubic_roots[i] + d) < 1e-10);
    }
  compare_vectors_approx(cubic_roots, {1+sqrt(3.),1-sqrt(3.),-2.});

}


TEST_CASE("Bezier curves")
{
  std::cout << "Bezier curves" << std::endl;
  std::vector<std::pair<double,std::vector<double>>> contours_1;
  std::vector<Point<2> > points = {Point<2>(0,4,CoordinateSystem::cartesian),
                                   Point<2>(1,3,CoordinateSystem::cartesian),
                                   Point<2>(2,1,CoordinateSystem::cartesian),
                                   Point<2>(3,2,CoordinateSystem::cartesian),
                                   Point<2>(4,3,CoordinateSystem::cartesian),
                                   Point<2>(5,3,CoordinateSystem::cartesian),
                                   Point<2>(5.,2,CoordinateSystem::cartesian)
                                  };

  std::vector<double> angles(points.size(),NAN);
  Objects::BezierCurve b1(points, angles);
  for (size_t i = 0; i < b1.control_points.size(); ++i)
    {
      std::cout << i << ": " << b1.control_points[i] << ", angle = " << b1.angles[i] << "(" << b1.angles[i] *180./M_PI << "), lengths = " << b1.lengths[i] << std::endl;
    }
  std::cout << "last angle = " << b1.angles[b1.control_points.size()]<< "(" << b1.angles[b1.control_points.size()] *180./M_PI << ")" << std::endl;;

  // now compute closest point on line
  Point<2> test_point(0.1,4,CoordinateSystem::cartesian);
  Objects::ClosestPointOnCurve result = b1.closest_point_on_curve_segment(test_point);
  std::cout << "closest point on curve to " << test_point << " = " << result.point << ", distance = " << result.distance<< std::endl << std::endl;

  test_point = Point<2>(0.2,4,CoordinateSystem::cartesian);
  result = b1.closest_point_on_curve_segment(test_point);
  std::cout << "closest point on curve to " << test_point << " = " << result.point << ", distance = " << result.distance << std::endl << std::endl;

  test_point = Point<2>(0.5,3.75,CoordinateSystem::cartesian);
  result = b1.closest_point_on_curve_segment(test_point);
  std::cout << "closest point on curve to " << test_point << " = " << result.point << ", distance = " << result.distance << std::endl << std::endl;

  test_point = Point<2>(1,3.1,CoordinateSystem::cartesian);
  result = b1.closest_point_on_curve_segment(test_point);
  std::cout << "closest point on curve to " << test_point << " = " << result.point << ", distance = " << result.distance << std::endl << std::endl;

  test_point = Point<2>(1,1,CoordinateSystem::cartesian);
  result = b1.closest_point_on_curve_segment(test_point);
  std::cout << "closest point on curve to " << test_point << " = " << result.point << ", distance = " << result.distance << std::endl << std::endl;

  test_point = Point<2>(2,2,CoordinateSystem::cartesian);
  result = b1.closest_point_on_curve_segment(test_point);
  std::cout << "closest point on curve to " << test_point << " = " << result.point << ", distance = " << result.distance << std::endl << std::endl;

  test_point = Point<2>(6,2.5,CoordinateSystem::cartesian);
  result = b1.closest_point_on_curve_segment(test_point);
  std::cout << "closest point on curve to " << test_point << " = " << result.point << ", distance = " << result.distance << std::endl << std::endl;

  test_point = Point<2>(4,2.5,CoordinateSystem::cartesian);
  result = b1.closest_point_on_curve_segment(test_point);
  std::cout << "closest point on curve to " << test_point << " = " << result.point << ", distance = " << result.distance << std::endl << std::endl;

  test_point = Point<2>(6,1.5,CoordinateSystem::cartesian);
  result = b1.closest_point_on_curve_segment(test_point);
  std::cout << "closest point on curve to " << test_point << " = " << result.point << ", distance = " << result.distance << std::endl << std::endl;

  test_point = Point<2>(6,1.85,CoordinateSystem::cartesian);
  result = b1.closest_point_on_curve_segment(test_point);
  std::cout << "closest point on curve to " << test_point << " = " << result.point << ", distance = " << result.distance << std::endl << std::endl;

  test_point = Point<2>(6,1.8,CoordinateSystem::cartesian);
  result = b1.closest_point_on_curve_segment(test_point);
  std::cout << "closest point on curve to " << test_point << " = " << result.point << ", distance = " << result.distance << std::endl << std::endl;
}
/*
TEST_CASE("Circle lines")
{
  std::cout << "Circe lines" << std::endl;
  std::vector<std::pair<double,std::vector<double>>> contours_1;
  std::vector<Point<2> > points = {Point<2>(0,4,CoordinateSystem::cartesian),
                                   Point<2>(1,3,CoordinateSystem::cartesian),
                                   Point<2>(2,1,CoordinateSystem::cartesian),
                                   Point<2>(3,2,CoordinateSystem::cartesian),
                                   Point<2>(4,3,CoordinateSystem::cartesian),
                                   Point<2>(5,3,CoordinateSystem::cartesian)
                                  };
  Utilities::CircleLine c1(points);

  for (size_t i = 0; i < c1.circle_centers.size(); ++i)
    {
      std::cout << i << ": " << c1.circle_centers[i] << ", angle = " << c1.angles[i] << std::endl;
    }
  std::cout << "last angle = " << c1.angles[c1.circle_centers.size()];
//for (auto &i : a.connectivity){
//  for (unsigned int j = 0; j < i.size(); j++){
//    std::cout << j << ": ";
//    for (unsigned int k : i[j]){
//      std::cout << " " << k;
//    }
//    std::cout << std::endl;
//  }
//}

}*/