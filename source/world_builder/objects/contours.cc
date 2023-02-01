/*
  Copyright (C) 2018 - 2023 by the authors of the World Builder code.

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
#include "world_builder/objects/contours.h"
#include "world_builder/assert.h"
#include "world_builder/utilities.h"

#include <iostream>


namespace WorldBuilder
{
  namespace Objects
  {
    Contours::Contours()
    {}

    Contours::Contours(const std::vector<std::vector<Point<2> > > &points_,
                       const std::vector<double> depths_,
                       const std::vector<std::vector<double> > &thicknesses_,
                       const double start_radius_,
                       const std::vector<std::vector<double> > &angle_contraints_,
                       const std::vector<std::vector<Point<2> > > &directions_)
      :
      points(points_),
      depths(depths_),
      angle_contraints(angle_contraints_),
      thicknesses(thicknesses_),
      directions(directions_),
      start_radius(start_radius_)
    {
      bool bool_cartesian = points[0][0].get_coordinate_system() == cartesian;
      const size_t n_curves = points.size();
      std::vector<std::vector<double> > angle_contraints_horizontal(points.size());
      for (size_t i = 0; i < points.size(); ++i){
        angle_contraints_horizontal[i].resize(points[i].size(),NaN::DQNAN);
      }

      // if angle contrainst is empty fill it with the right amount of entries and set the values to signaling nan.
      if (angle_contraints.size() == 0)
        {
          angle_contraints = angle_contraints_horizontal;
        }

      WBAssertThrow(angle_contraints.size() == points.size(),
                    "Error: incorrect number of curves in angle contstrains: " << angle_contraints.size() << ". Expected: " << points.size());
      for (size_t curve_i = 0; curve_i < angle_contraints.size(); ++curve_i)
        WBAssertThrow(angle_contraints[curve_i].size() == points[curve_i].size(),
                      "Error: incorrect number of points in curve " << curve_i << " in angle contstrains: " << angle_contraints[curve_i].size()
                      << ". Expected: " << points[curve_i].size());

      for (size_t curve_i = 0; curve_i < n_curves; ++curve_i)
        {
          contour_curves.emplace_back(BezierCurve(points[curve_i],angle_contraints_horizontal[curve_i]));
          std::cout << "     angles (" << contour_curves[curve_i].angles.size() << ") " << curve_i << ":" << std::endl;
          for (size_t angle_i = 0; angle_i < contour_curves[curve_i].angles.size(); ++angle_i)
            {
              std::cout << contour_curves[curve_i].angles[angle_i] << " (" << contour_curves[curve_i].angles[angle_i]*180/Consts::PI << "), ";
            }
          std::cout  << std::endl;
        }

      max_thickness_on_curve.resize(n_curves);
      for (size_t curve_i = 0; curve_i < n_curves; ++curve_i)
        {
          std::cout << "curve " << curve_i << ": thicknesses[curve_i].size() = " << thicknesses[curve_i].size() << std::endl;
          for (size_t section_i = 0; section_i < thicknesses[curve_i].size(); ++section_i)
            {
              std::cout << "thicknesses[" << curve_i << "][" << section_i << "]" << thicknesses[curve_i][section_i]
                        << ",  max_thickness_on_curve[curve_i] = " <<  max_thickness_on_curve[curve_i] << std::endl;
              if (thicknesses[curve_i][section_i] > max_thickness_on_curve[curve_i])
                {
                  max_thickness_on_curve[curve_i] = thicknesses[curve_i][section_i];
                }
            }
        }


      // if directions is empty, set it to the right amount of entries (the last curve doesn't have a direction)
      // and the the values to of the closest point below.
      // If there is only one set it to the value. Otherwise it is assumed to have the correct number of curves and points
      // (which is checked afterwards).
      if (directions.size() == 0)
        {
          directions.resize(angle_contraints.size()-1, {});
          for (size_t curve_i = 0; curve_i < directions.size(); ++curve_i)
            {
              for (size_t point_i = 0; point_i < angle_contraints[curve_i].size(); ++point_i)
                {
                  // find the closest point below and point to it
                  directions[curve_i].emplace_back(contour_curves[curve_i+1].closest_point_on_curve_segment(points[0][point_i]).point);
                }
            }
        }
      else if (directions.size() == 1 && directions[0].size() == 0)
        {
          directions.resize(angle_contraints.size()-1, {});
          for (size_t curve_i = 0; curve_i < directions.size(); ++curve_i)
            directions[curve_i].resize(angle_contraints[curve_i].size(),directions[0][0]);
        }


      WBAssertThrow(directions.size() == points.size()-1,
                    "Error: incorrect number of curves in directions: " << directions.size() << ". Expected: " << points.size()-1);
      for (size_t curve_i = 0; curve_i < directions.size(); ++curve_i)
        {
          WBAssertThrow(directions[curve_i].size() == points[curve_i].size(),
                        "Error: incorrect number of points in curve " << curve_i << " in directions: " << directions[curve_i].size()
                        << ". Expected: " << points[curve_i].size());
        }

      /*std::cout << "============================================= construct connectivity ========================================" << std::endl << std::endl;
      WBAssertThrow(points.size() >= 2, "The countours object needs at least 2 conours.");
      // first create the connectivity
      // The connectivity is created by connecting each point with one or many points in a countour below. One poiont below can have many connections to a point above.
      // Connect each point with the closest point on the contour below.
      // Assumptions: Countours always get deeper and the points are ordered in the same direction.
      //connectivity.resize(points.size()-1);
      connectivity.resize(n_curves);
      for (unsigned int contour_i = 0; contour_i < n_curves-1; contour_i++)
        {
          std::cout << " ---------------- " << contour_i << "--------------------- " << std::endl;

          connectivity[contour_i].resize(points[contour_i].size());
          // start with connecting the edge to prevent an if statement within the loop
          connectivity[contour_i][0].emplace_back(0);

          // now start connecting the rest
          double depth = depths[contour_i];
          double depth_next = depths[contour_i+1];
          unsigned int countour_entry_i = 0;
          unsigned int next_countour_entry_i = 0;
          std::cout << "start search for contour " << contour_i << ", contours[contour_i].second.size() = " << points[contour_i].size() << ",  contours[contour_i+1].second.size() = " <<  (points[contour_i+1].size())-1 << std::endl;
          while (countour_entry_i < points[contour_i].size() && next_countour_entry_i < points[contour_i+1].size())
            {
              std::cout << "countour_entry_i: " << countour_entry_i << ", next_countour_entry_i: " << next_countour_entry_i << std::endl;
              // compute distances
              // between the current top point and the next down point.
              Point<2> current_point = points[contour_i][countour_entry_i];
              Point<2> next_point = points[contour_i][countour_entry_i+1];
              //double x = contours[contour_i].second[2*countour_entry_i];
              //double y = contours[contour_i].second[2*countour_entry_i+1];
              //double x_next = contours[contour_i+1].second[2*(next_countour_entry_i+1)];
              //double y_next = contours[contour_i+1].second[2*(next_countour_entry_i+1)+1];

              std::cout << "x:y = " << x << ":" << y << ", next: " << x_next << ":" << y_next << std::endl;

              double distance_1 = std::sqrt((depth-depth_next)*(depth-depth_next) + (x-x_next)*(x-x_next) + (y-y_next)*(y-y_next));
              double distance_2 = std::numeric_limits<double>::infinity();
              double distance_3 = std::numeric_limits<double>::infinity();
              bool countour_next_entry = false;
              bool next_countour_next_entry = false;
              if (contours[contour_i].second.size() > 2*(countour_entry_i+1))
                {
                  // between the next top point and the current down point.
                  x = contours[contour_i].second[2*(countour_entry_i+1)];
                  y = contours[contour_i].second[2*(countour_entry_i+1)+1];
                  x_next = contours[contour_i+1].second[2*next_countour_entry_i];
                  y_next = contours[contour_i+1].second[2*next_countour_entry_i+1];

                  distance_2 = std::sqrt((depth-depth_next)*(depth-depth_next) + (x-x_next)*(x-x_next) + (y-y_next)*(y-y_next));

                  std::cout << "x:y = " << x << ":" << y << ", next: " << x_next << ":" << y_next << std::endl;

                  if (distance_1 > distance_2)
                    countour_next_entry = true;
                  // between the next top point and the next down point.
                  x = contours[contour_i].second[2*(countour_entry_i+1)];
                  y = contours[contour_i].second[2*(countour_entry_i+1)+1];
                  x_next = contours[contour_i+1].second[2*(next_countour_entry_i+1)];
                  y_next = contours[contour_i+1].second[2*(next_countour_entry_i+1)+1];
                  distance_3 = std::sqrt((depth-depth_next)*(depth-depth_next) + (x-x_next)*(x-x_next) + (y-y_next)*(y-y_next));
                  std::cout << "x:y = " << x << ":" << y << ", next: " << x_next << ":" << y_next << std::endl;
                }
              else
                {
                  next_countour_next_entry = true;
                  std::cout <<  "distance 1 = " << distance_1 << ", distance 2 = " << distance_2 << "distance 3 = " << distance_3 << std::endl;
                }
              // I don't think that if distance_1 is shorter than distance_2, that distance
              if (distance_1> distance_3 && distance_2> distance_3)
                {
                  countour_next_entry = true;
                  next_countour_next_entry = true;
                }
              if (distance_2> distance_3)
                next_countour_next_entry = true;

              if (countour_next_entry)
                countour_entry_i += 1;
              if (next_countour_next_entry)
                next_countour_entry_i += 1;
              if (next_countour_entry_i > contours[contour_i+1].second.size()/2.-1)
                {
                  break;
                }
              connectivity[contour_i][countour_entry_i].emplace_back(next_countour_entry_i);
            }
        }*/

      std::cout << "============================================= compute vertical angle at each point ========================================" << std::endl << std::endl;
      // first deterimine the angles for curve 0
      if (points.size() > 2)
        {
          std::cout << "ct flag 50" << std::endl;
          // first set the angles at the top
          for (size_t point_i = 0; point_i < angle_contraints[0].size(); ++point_i)
            {
              std::cout << "ct flag 51" << std::endl;
              if (std::isnan(angle_contraints[0][point_i]))
                {
                  std::cout << "ct flag 52: points[0][point_i].coordiante_system()=" << points[0][point_i].get_coordinate_system() << std::endl;
                  // The first angle at each point is determined by the angle of the closest point below and the closest point below that
                  ClosestPointOnCurve closest_point_below = contour_curves[1].closest_point_on_curve_segment(points[0][point_i]);
                  std::cout << "ct flag 52.5: points[0][point_i].coordiante_system()=" << points[0][point_i].get_coordinate_system() << ", closest_point_below.point.get_coordinate_system()= " << closest_point_below.point.get_coordinate_system() << ", closest_point_below= " << closest_point_below.point<< std::endl;

                  ClosestPointOnCurve closest_point_below_below = contour_curves[2].closest_point_on_curve_segment(closest_point_below.point);

                  std::cout << "ct flag 53" << std::endl;
                  // first we need to compute the normal to the plane defined by the vector going from the point down and to closest_point_below_below.
                  const Point<2> &basis_point = points[0][point_i];
                  const Point<2> &point_b = closest_point_below.point;
                  const Point<2> &point_bb = closest_point_below_below.point;

                  std::cout << "ct flag 54" << std::endl;
                  const Point<3> basis_point_cartesian = bool_cartesian ? Point<3>(basis_point[0],basis_point[1],start_radius-depths[0],cartesian) : Utilities::spherical_to_cartesian_coordinates(Point<3>(start_radius-depths[0],basis_point[0],basis_point[1],spherical).get_array());
                  const Point<3> point_b_cartesian = bool_cartesian ? Point<3>(point_b[0],point_b[1],start_radius-depths[1],cartesian) : Utilities::spherical_to_cartesian_coordinates(Point<3>(start_radius-depths[1],point_b[0],point_b[1],spherical).get_array());
                  const Point<3> point_bb_cartesian = bool_cartesian ? Point<3>(point_bb[0],point_bb[1],start_radius-depths[2],cartesian) : Utilities::spherical_to_cartesian_coordinates(Point<3>(start_radius-depths[2],point_bb[0],point_bb[1],spherical).get_array());
                  const Point<3> basis_point_bottom_cartesian = bool_cartesian ? Point<3>(basis_point[0],basis_point[1],0,cartesian) : Utilities::spherical_to_cartesian_coordinates(Point<3>(0,basis_point[0],basis_point[1],spherical).get_array());
                  const Point<3> basis_point_bottom_cartesian_normalized = basis_point_bottom_cartesian/basis_point_bottom_cartesian.norm();

                  std::cout << "ct flag 55" << std::endl;
                  Point<3> normal_to_plane = Utilities::cross_product(point_bb_cartesian-basis_point_cartesian,basis_point_bottom_cartesian-basis_point_cartesian);
                  normal_to_plane = normal_to_plane/normal_to_plane.norm();

                  std::cout << "ct flag 56" << std::endl;
                  // next compute the x and y axis
                  const Point<3> y_axis = (basis_point_cartesian-basis_point_bottom_cartesian)/(basis_point_cartesian-basis_point_bottom_cartesian).norm();
                  const Point<3> x_axis = Utilities::cross_product(y_axis,normal_to_plane)/Utilities::cross_product(y_axis,normal_to_plane).norm();

                  std::cout << "ct flag 57" << std::endl;
                  // compute basis_point and point_bb in 2d cross_section
                  Point<2> basis_point_2d(x_axis * (basis_point_cartesian-basis_point_bottom_cartesian),
                                          y_axis * (basis_point_cartesian-basis_point_bottom_cartesian),
                                          cartesian);

                  std::cout << "ct flag 58" << std::endl;
                  Point<2> point_b_2d(x_axis * (point_b_cartesian-basis_point_bottom_cartesian),
                                      y_axis * (point_b_cartesian-basis_point_bottom_cartesian),
                                      cartesian);

                  std::cout << "ct flag 59" << std::endl;
                  Point<2> point_bb_2d(x_axis * (point_bb_cartesian-basis_point_bottom_cartesian),
                                       y_axis * (point_bb_cartesian-basis_point_bottom_cartesian),
                                       cartesian);

                  std::cout << "ct flag 60" << std::endl;
                  const double top_b_angle = atan2((point_b_2d-basis_point_2d)[1],(point_b_2d-basis_point_2d)[0])+M_PI;
                  const double top_bb_angle = atan2((point_bb_2d-basis_point_2d)[1],(point_bb_2d-basis_point_2d)[0])+M_PI;
                  const double diff_angle = top_b_angle-top_bb_angle;

                  std::cout << "ct flag 61" << std::endl;
                  angle_contraints[0][point_i] = top_b_angle+0.5*diff_angle;



                  std::cout  << ", basis_point: " << basis_point << ", point_b: " << point_b  << ", point_bb: " << point_bb
                             << ", x_axis = " << x_axis << ", y_axis = " << y_axis
                             << ", basis_point_2d = " << basis_point_2d << ", point_b_2d = " << point_b_2d  << ", point_bb_2d = " << point_bb_2d
                             << ", top_b_angle = " << top_b_angle << "(" << (top_b_angle) *180/M_PI << ")"
                             << ", top_bb_angle = " << top_bb_angle << "(" << (top_bb_angle) *180/M_PI << ")"
                             << ", diff_angle = " << diff_angle << "(" << (diff_angle) *180/M_PI << ")"
                             << ", top_angle = " << angle_contraints[0][point_i] << "(" << (angle_contraints[0][point_i]) *180/M_PI << ")" << std::endl;
                }
            }

          std::cout << "-----------------------------" << std::endl;
          // next loop over all curves up to the last one and set the angles
          for (size_t curve_i = 1; curve_i < angle_contraints.size()-1; ++curve_i)
            {
              for (size_t point_i = 0; point_i < angle_contraints[curve_i].size(); ++point_i)
                {
                  if (std::isnan(angle_contraints[curve_i][point_i]))
                    {
                      // The first angle at each point is determined by the angle of the closest point below and the closest point below that
                      ClosestPointOnCurve closest_point_above = contour_curves[curve_i-1].closest_point_on_curve_segment(points[curve_i][point_i]);
                      ClosestPointOnCurve closest_point_below = contour_curves[curve_i+1].closest_point_on_curve_segment(points[curve_i][point_i]);

                      // first we need to compute the normal to the plane defined by the vector going from the point down and to closest_point_below_below.
                      const Point<2> &basis_point = points[curve_i][point_i];
                      const Point<2> &point_a = closest_point_above.point;
                      const Point<2> &point_b = closest_point_below.point;

                      const Point<3> basis_point_cartesian = bool_cartesian ? Point<3>(basis_point[0],basis_point[1],start_radius-depths[curve_i],cartesian) : Utilities::spherical_to_cartesian_coordinates(Point<3>(start_radius-depths[curve_i],basis_point[0],basis_point[1],spherical).get_array());
                      const Point<3> point_a_cartesian = bool_cartesian ? Point<3>(point_a[0],point_a[1],start_radius-depths[curve_i-1],cartesian) : Utilities::spherical_to_cartesian_coordinates(Point<3>(start_radius-depths[curve_i-1],point_a[0],point_a[1],spherical).get_array());
                      const Point<3> point_b_cartesian = bool_cartesian ? Point<3>(point_b[0],point_b[1],start_radius-depths[curve_i+1],cartesian) : Utilities::spherical_to_cartesian_coordinates(Point<3>(start_radius-depths[curve_i+1],point_b[0],point_b[1],spherical).get_array());
                      const Point<3> basis_point_bottom_cartesian = bool_cartesian ? Point<3>(basis_point[0],basis_point[1],0,cartesian) : Utilities::spherical_to_cartesian_coordinates(Point<3>(0,basis_point[0],basis_point[1],spherical).get_array());
                      const Point<3> basis_point_bottom_cartesian_normalized = basis_point_bottom_cartesian/basis_point_bottom_cartesian.norm();

                      Point<3> normal_to_plane = Utilities::cross_product(point_b_cartesian-basis_point_cartesian,basis_point_bottom_cartesian-basis_point_cartesian);
                      normal_to_plane = normal_to_plane/normal_to_plane.norm();

                      // next compute the x and y axis
                      const Point<3> y_axis = (basis_point_cartesian-basis_point_bottom_cartesian)/(basis_point_cartesian-basis_point_bottom_cartesian).norm();
                      const Point<3> x_axis = Utilities::cross_product(y_axis,normal_to_plane)/Utilities::cross_product(y_axis,normal_to_plane).norm();

                      Point<2> basis_point_2d(x_axis * (basis_point_cartesian-basis_point_bottom_cartesian),
                                              y_axis * (basis_point_cartesian-basis_point_bottom_cartesian),
                                              cartesian);

                      Point<2> point_a_2d(x_axis * (point_a_cartesian-basis_point_bottom_cartesian),
                                          y_axis * (point_a_cartesian-basis_point_bottom_cartesian),
                                          cartesian);

                      Point<2> point_b_2d(x_axis * (point_b_cartesian-basis_point_bottom_cartesian),
                                          y_axis * (point_b_cartesian-basis_point_bottom_cartesian),
                                          cartesian);

                      //angle_contraints[curve_i][point_i] = atan2((point_a_2d-point_b_2d)[1],(point_a_2d-point_b_2d)[0]);
                      const double top_a_angle = atan2((point_a_2d-basis_point_2d)[1],(point_a_2d-basis_point_2d)[0]);//+M_PI;
                      const double top_b_angle = atan2((point_b_2d-basis_point_2d)[1],(point_b_2d-basis_point_2d)[0])+M_PI;
                      const double diff_angle = top_a_angle-top_b_angle;

                      angle_contraints[curve_i][point_i] = (top_a_angle+top_b_angle)*0.5;


                      std::cout  << ", basis_point: " << basis_point << ", point_a: " << point_a  << ", point_b: " << point_b
                                 << ", x_axis = " << x_axis << ", y_axis = " << y_axis
                                 << ", basis_point_2d = " << basis_point_2d << ", point_a_2d = " << point_a_2d  << ", point_b_2d = " << point_b_2d
                                 << ", top_a_angle = " << top_a_angle << "(" << (top_a_angle) *180/M_PI << ")"
                                 << ", top_b_angle = " << top_b_angle << "(" << (top_b_angle) *180/M_PI << ")"
                                 << ", diff_angle = " << diff_angle << "(" << (diff_angle) *180/M_PI << ")"
                                 << ", angle = " << angle_contraints[curve_i][point_i] << "(" << (angle_contraints[curve_i][point_i]) *180/M_PI << ")" << std::endl;
                    }
                }
            }

          std::cout << "-----------------------------" << std::endl;
          // lastly loop over the points of the last cruve
          const size_t last_index = angle_contraints.size()-1;
          for (size_t point_i = 0; point_i < angle_contraints[last_index].size(); ++point_i)
            {
              if (std::isnan(angle_contraints[last_index][point_i]))
                {
                  // The first angle at each point is determined by the angle of the closest point below and the closest point below that
                  ClosestPointOnCurve closest_point_above = contour_curves[last_index-1].closest_point_on_curve_segment(points[last_index][point_i]);
                  ClosestPointOnCurve closest_point_above_above = contour_curves[last_index-2].closest_point_on_curve_segment(closest_point_above.point);

                  // first we need to compute the normal to the plane defined by the vector going from the point down and to closest_point_below_below.
                  const Point<2> &basis_point = points[last_index][point_i];
                  const Point<2> &point_a = closest_point_above.point;
                  const Point<2> &point_aa = closest_point_above_above.point;

                  const Point<3> basis_point_cartesian = bool_cartesian ? Point<3>(basis_point[0],basis_point[1],start_radius-depths[last_index],cartesian) : Utilities::spherical_to_cartesian_coordinates(Point<3>(start_radius-depths[last_index],basis_point[0],basis_point[1],spherical).get_array());
                  const Point<3> point_a_cartesian = bool_cartesian ? Point<3>(point_a[0],point_a[1],start_radius-depths[last_index-1],cartesian) : Utilities::spherical_to_cartesian_coordinates(Point<3>(start_radius-depths[last_index-1],point_a[0],point_a[1],spherical).get_array());
                  const Point<3> point_aa_cartesian = bool_cartesian ? Point<3>(point_aa[0],point_aa[1],start_radius-depths[last_index-2],cartesian) : Utilities::spherical_to_cartesian_coordinates(Point<3>(start_radius-depths[last_index-2],point_aa[0],point_aa[1],spherical).get_array());
                  const Point<3> basis_point_bottom_cartesian = bool_cartesian ? Point<3>(basis_point[0],basis_point[1],0,cartesian) : Utilities::spherical_to_cartesian_coordinates(Point<3>(0,basis_point[0],basis_point[1],spherical).get_array());
                  const Point<3> basis_point_bottom_cartesian_normalized = basis_point_bottom_cartesian/basis_point_bottom_cartesian.norm();

                  Point<3> normal_to_plane = Utilities::cross_product(point_aa_cartesian-basis_point_cartesian,basis_point_bottom_cartesian-basis_point_cartesian);
                  normal_to_plane = normal_to_plane/normal_to_plane.norm();

                  // next compute the x and y axis
                  const Point<3> y_axis = (basis_point_cartesian-basis_point_bottom_cartesian)/(basis_point_cartesian-basis_point_bottom_cartesian).norm();
                  const Point<3> x_axis = Utilities::cross_product(y_axis,normal_to_plane)/Utilities::cross_product(y_axis,normal_to_plane).norm();

                  // compute basis_point and point_aa in 2d cross_section
                  Point<2> basis_point_2d(x_axis * (basis_point_cartesian-basis_point_cartesian),
                                          y_axis * (basis_point_cartesian-basis_point_cartesian),
                                          cartesian);

                  Point<2> point_a_2d(x_axis * (point_a_cartesian-basis_point_cartesian),
                                      y_axis * (point_a_cartesian-basis_point_cartesian),
                                      cartesian);

                  Point<2> point_aa_2d(x_axis * (point_aa_cartesian-basis_point_cartesian),
                                       y_axis * (point_aa_cartesian-basis_point_cartesian),
                                       cartesian);

                  const double top_a_angle = -atan2((point_a_2d-basis_point_2d)[1],(point_a_2d-basis_point_2d)[0])+M_PI;
                  const double top_aa_angle = -atan2((point_aa_2d-basis_point_2d)[1],(point_aa_2d-basis_point_2d)[0])+M_PI;
                  //const double top_a_angle = atan2((basis_point_2d-point_aa_2d)[1],(basis_point_2d-point_aa_2d)[0]);
                  //const double top_aa_angle = atan2((point_a_2d-point_aa_2d)[1],(point_a_2d-point_aa_2d)[0]);
                  const double diff_angle = top_a_angle-top_aa_angle;

                  angle_contraints[last_index][point_i] = top_a_angle+0.5*diff_angle;



                  std::cout  << ", basis_point: " << basis_point << ", point_a: " << point_a  << ", point_aa: " << point_aa
                             << ", x_axis = " << x_axis << ", y_axis = " << y_axis
                             << ", basis_point_2d = " << basis_point_2d << ", point_a_2d-b = " << point_a_2d-basis_point_2d  << ", point_aa_2d-b = " << point_aa_2d-basis_point_2d
                             << ", top_a_angle = " << top_a_angle << "(" << (top_a_angle) *180/M_PI << ")"
                             << ", top_aa_angle = " << top_aa_angle << "(" << (top_aa_angle) *180/M_PI << ")"
                             << ", diff_angle = " << diff_angle << "(" << (diff_angle) *180/M_PI << ")"
                             << ", top_angle = " << angle_contraints[last_index][point_i] << "(" << (angle_contraints[last_index][point_i]) *180/M_PI << ")" << std::endl;
                }
            }
        }
      else
        {
          WBAssertThrow(false, "you provided only two lines, this has not been implemented yet.");
        }


      std::cout << "============================================= construct distance along surface ========================================" << std::endl << std::endl;
      distance_along_surface.resize(n_curves);
      // the distance on curve_i == 0 is zero, so start with curve 1
      distance_along_surface[0].resize(points[0].size(),0.);
      for (size_t curve_i = 1; curve_i < n_curves; ++curve_i)
        {
          std::cout << "ct flag 10" << std::endl;
          distance_along_surface[curve_i].resize(points[curve_i].size(),0.);
          for (size_t point_i = 0; point_i < points[curve_i].size(); ++point_i)
            {
              std::cout << "ct flag 11" << std::endl;

              const Point<3> basis_point_3D_cartesian = bool_cartesian
                                                        ?
                                                        Point<3>(points[curve_i][point_i][0],points[curve_i][point_i][1],start_radius-depths[curve_i],cartesian)
                                                        :
                                                        Utilities::spherical_to_cartesian_coordinates(Point<3>(start_radius-depths[curve_i],points[curve_i][point_i][0],points[curve_i][point_i][1],spherical).get_array());

              std::cout << "ct flag 12" << std::endl;
              const Point<3> basis_point_3D_cartesian_bottom = bool_cartesian
                                                               ?
                                                               Point<3>(points[curve_i][point_i][0],points[curve_i][point_i][1],0.,cartesian)
                                                               :
                                                               Utilities::spherical_to_cartesian_coordinates(Point<3>(0.,points[curve_i][point_i][0],points[curve_i][point_i][1],spherical).get_array());
              // loop over the points above and find the point for which the direction (3D) makes the smallest angle a direct line between the two points.
              // if they are the same (with an error of 1 degree?), use the smallest distance;

              std::cout << "ct flag 13" << std::endl;
              unsigned int lowest_angle = 800;
              size_t lowest_angle_point = std::numeric_limits<size_t>::max();
              double lowest_angle_distance = std::numeric_limits<double>::infinity();


              std::cout << "ct flag 14" << std::endl;
              for (size_t prev_point_i = 0; prev_point_i < points[curve_i-1].size(); ++prev_point_i)
                {

                  std::cout << "ct flag 15" << std::endl;
                  const Point<3> basis_point_3D_cartesian_up = bool_cartesian
                                                               ?
                                                               Point<3>(points[curve_i-1][prev_point_i][0],points[curve_i-1][prev_point_i][1],start_radius-depths[curve_i-1],cartesian)
                                                               :
                                                               Utilities::spherical_to_cartesian_coordinates(Point<3>(start_radius-depths[curve_i-1],points[curve_i-1][prev_point_i][0],points[curve_i-1][prev_point_i][1],spherical).get_array());

                  std::cout << "ct flag 16" << std::endl;
                  const Point<3> direction_3d_cartesian = (bool_cartesian
                                                           ?
                                                           Point<3>(directions[curve_i-1][prev_point_i][0],directions[curve_i-1][prev_point_i][1],start_radius-depths[curve_i], cartesian)
                                                           :
                                                           Utilities::spherical_to_cartesian_coordinates(Point<3>(start_radius-depths[curve_i],directions[curve_i-1][prev_point_i][0],directions[curve_i-1][prev_point_i][1],spherical).get_array()))
                                                          - basis_point_3D_cartesian_up;

                  std::cout << "ct flag 16" << std::endl;
                  const Point<3> top_bottom_3d_cartesian = (bool_cartesian
                                                            ?
                                                            Point<3>(directions[curve_i-1][prev_point_i][0],directions[curve_i-1][prev_point_i][1],start_radius-depths[curve_i], cartesian)
                                                            :
                                                            Utilities::spherical_to_cartesian_coordinates(Point<3>(start_radius-depths[curve_i],directions[curve_i-1][prev_point_i][0],directions[curve_i-1][prev_point_i][1],spherical).get_array()))
                                                           - basis_point_3D_cartesian_up;

                  std::cout << "ct flag 17" << std::endl;
                  const unsigned int angle = (unsigned int) std::round(std::abs(std::acos((direction_3d_cartesian*top_bottom_3d_cartesian)/(direction_3d_cartesian.norm()*top_bottom_3d_cartesian.norm()))*180/M_PI));
                  if (angle <= lowest_angle)
                    {
                      // compute the distance between the previous point and this point through a Bezier curve


                      Point<3> normal_to_plane = Utilities::cross_product(basis_point_3D_cartesian-basis_point_3D_cartesian_up,basis_point_3D_cartesian_bottom-basis_point_3D_cartesian_up);
                      normal_to_plane = normal_to_plane/normal_to_plane.norm();

                      // next compute the x and y axis
                      const Point<3> y_axis = (basis_point_3D_cartesian_up-basis_point_3D_cartesian_bottom)/(basis_point_3D_cartesian_up-basis_point_3D_cartesian_bottom).norm();
                      const Point<3> x_axis = Utilities::cross_product(y_axis,normal_to_plane)/Utilities::cross_product(y_axis,normal_to_plane).norm();

                      Point<2> basis_point_2D_up(x_axis * (basis_point_3D_cartesian_up-basis_point_3D_cartesian_up),
                                                 y_axis * (basis_point_3D_cartesian_up-basis_point_3D_cartesian_up),
                                                 cartesian);

                      Point<2> basis_point_2D(x_axis * (basis_point_3D_cartesian-basis_point_3D_cartesian_up),
                                              y_axis * (basis_point_3D_cartesian-basis_point_3D_cartesian_up),
                                              cartesian);

                      const std::vector<Point<2> > local_points = {basis_point_2D_up,basis_point_2D};
                      //const std::vector<double> local_angle_contraints = {angle_contraints[curve_i-1][prev_point_i],angle_contraints[curve_i][point_i]};

                      const double curve_length = (basis_point_2D_up-basis_point_2D).norm();//BezierCurve(local_points,local_angle_contraints).lengths[0]; // todo: improve
                      // replace the old value with the new values if the angle is lower
                      // or the angle is equal to the lowest angle bt the distance is lower.
                      if (angle < lowest_angle || curve_length < lowest_angle_distance)
                        {

                          // use this previous point
                          lowest_angle = angle;
                          lowest_angle_point = prev_point_i;
                          lowest_angle_distance = curve_length;
                        }
                    }
                }

              distance_along_surface[curve_i][point_i] = distance_along_surface[curve_i-1][lowest_angle_point]+lowest_angle_distance;
              /*
              // find the closest points on the previous curve.
              std::cout << "flag 1: curve_i: " << curve_i << ", point_i: " << point_i << ", point: " << points[curve_i][point_i] << std::endl;
              ClosestPointOnCurve closest_point_above = contour_curves[curve_i-1].closest_point_on_curve_segment(points[curve_i][point_i]);
              std::cout << "flag 2: closest_point_above.index: " << closest_point_above.index
                        << ", distance_along_surface.size(): " << distance_along_surface.size()
                        << ", distance_along_surface[curve_i].size(): " << distance_along_surface[curve_i].size()
                        << ", distance_along_surface[curve_i-1].size(): " << distance_along_surface[curve_i-1].size() << std::endl;
              std::cout << "flag 3: " << distance_along_surface[curve_i][point_i] << ", closest_point_above.point = " << closest_point_above.point << std::endl;

              //const double angle_above = angle_contraints[curve_i-1][]
              // find the closest point on the next curve
              //if (curve_i < contour_curves.size())
              //  {
              //    std::cout << "flag 1: curve_i: " << curve_i << ", point_i: " << point_i << ", point: " << points[curve_i][point_i] << std::endl;
              //    ClosestPointOnCurve closest_point_below = contour_curves[curve_i+1].closest_point_on_curve_segment(points[curve_i][point_i]);
              //    std::cout << "flag 2: closest_point_below.index: " << closest_point_below.index
              //              << ", distance_along_surface.size(): " << distance_along_surface.size()
              //              << ", distance_along_surface[curve_i].size(): " << distance_along_surface[curve_i].size()
              //              << ", distance_along_surface[curve_i-1].size(): " << distance_along_surface[curve_i-1].size() << std::endl;
              //    std::cout << "flag 3: " << distance_along_surface[curve_i][point_i] << ", closest_point_below.point = " << closest_point_below.point << std::endl;
              //  }
              //else
              //  {
              //    // This is the last curve. Compute the angle
              //  }


              // now find the index point on the curve below with the lowest angle (2d from above) with the line between the point and the closest
              // todo: if there is no curve below, skip.
              //for (size_t below_point_i = 0; below_point_i < points[curve_i+1].size(); ++below_point_i)
              //  {
              //  }

              // compute the 2D cross section and the location of the point in that 3D cross section
              Point<3> closest_point_above_3D = Point<3>(bool_cartesian ? closest_point_above.point[0] : start_radius-depths[curve_i-1],closest_point_above.point[1],bool_cartesian ? start_radius-depths[curve_i-1] : closest_point_above.point[2], cartesian);
              Point<3> point_3D = Point<3>(bool_cartesian ? points[curve_i][point_i][0] : start_radius-depths[curve_i-1],points[curve_i][point_i][1],bool_cartesian ? start_radius-depths[curve_i-1] : points[curve_i][point_i][2], cartesian);


              Point<3> closest_point_above_bottom = Point<3>(bool_cartesian ? closest_point_above.point[0] : 0.,closest_point_above.point[1],bool_cartesian ? 0 : closest_point_above.point[2], cartesian);
              //closest_point_on_line_bottom[bool_cartesian ? 2 : 0] = 0;
              std::cout << "compute compute_cross_section_axes" << std::endl;
              std::pair<Point<3>,Point<3> > axis = compute_cross_section_axes(closest_point_above_3D,point_3D,closest_point_above_bottom);

              // begin segment is {0,0} by definition?
              // now that we have the x and y axes computed, convert the 3d check point into a 2d one.
              //Point<2> begin_segment(0., 0., cartesian);
              //Point<2> end_segment(x_axis * points[curve_i][point_i], y_axis * points[curve_i][point_i], cartesian);
              std::vector<Point<2> > segment_points =
              {
                Point<2>(0., 0., cartesian),
                Point<2>(axis.first *point_3D, axis.second *point_3D, cartesian)
              };

              std::vector<double> segment_angle_constrains =
              {
                contour_curves[curve_i-1].angles[closest_point_above.index]+contour_curves[curve_i-1].angles[closest_point_above.index+1] *closest_point_above.interpolation_fraction,
                contour_curves[curve_i].angles[point_i]
              };
              std::cout << "segment_angle_constrains: " << segment_angle_constrains[0] << ":" << segment_angle_constrains[1] << std::endl;

              BezierCurve segment_curve = BezierCurve(segment_points,segment_angle_constrains);

              std::cout << "flag 3.2: segment length = " << segment_curve.lengths[0] <<std::endl;

              if (closest_point_above.index >= distance_along_surface[curve_i-1].size())
                {
                  std::cout << "flag 4a:";
                  distance_along_surface[curve_i][point_i] = distance_along_surface[curve_i-1][closest_point_above.index]+std::abs(closest_point_above.distance);
                }
              else
                {
                  std::cout << "flag 4b: " << distance_along_surface[curve_i-1][closest_point_above.index] << std::endl;
                  std::cout << "flag 5: " << distance_along_surface[curve_i-1][closest_point_above.index+1] << std::endl;
                  std::cout << "flag 6: " << closest_point_above.interpolation_fraction << std::endl;
                  distance_along_surface[curve_i][point_i] = distance_along_surface[curve_i-1][closest_point_above.index]+distance_along_surface[curve_i-1][closest_point_above.index+1]*closest_point_above.interpolation_fraction+std::abs(closest_point_above.distance);
                }*/
            }

        }


      std::cout << "angles: " << std::endl;
      for (size_t curve_i = 0; curve_i < angle_contraints.size(); ++curve_i)
        {
          for (size_t point_i = 0; point_i < angle_contraints[curve_i].size(); ++point_i)
            {
              std::cout << angle_contraints[curve_i][point_i]*180/M_PI  << " ";
            }
          std::cout << std::endl;
        }
      std::cout << std::endl;



      std::cout << "directions: " << std::endl;
      for (size_t curve_i = 0; curve_i < directions.size(); ++curve_i)
        {
          for (size_t point_i = 0; point_i < directions[curve_i].size(); ++point_i)
            {
              std::cout << directions[curve_i][point_i] << ", ";
            }
          std::cout << std::endl;
        }
      std::cout << std::endl;


      std::cout << "along surface lengths: " << std::endl;
      for (size_t curve_i = 0; curve_i < distance_along_surface.size(); ++curve_i)
        {
          for (size_t point_i = 0; point_i < distance_along_surface[curve_i].size(); ++point_i)
            {
              std::cout << distance_along_surface[curve_i][point_i] << ", ";
            }
          std::cout << std::endl;
        }
      std::cout << std::endl;

      //for (size_t curve_i = 0; curve_i < n_curves; ++curve_i)
      //  {
      //    double max_thickness_on_line = 0.;
      //    for (size_t section_i; section_i < thicknesses[curve_i].size(); ++section_i){
      //      if(thicknesses[curve_i][section_i] > max_thickness_on_line){
      //        max_thickness_on_line = thicknesses[curve_i][section_i];
      //      }
      //    }
      //    depth_ranges.emplace_back(depths[curve_i]-max_thickness_on_line,depths[curve_i+1]+max_thickness_on_line);
      //  }
      std::cout << "============================================= finished construction ========================================" << std::endl << std::endl;
    }

    std::pair<Point<3>,Point<3> >
    Contours::compute_cross_section_axes(Point<3> origin, Point<3> x_direction,Point<3> y_direction) const
    {
      Point<3> y_axis = origin - y_direction;
      Point<3> x_axis = origin - x_direction;
      // corner case:
      // The origin and the y direction are the same.
      // This is preventable so we assert this situation.
      WBAssert(y_axis.norm() > 2e-14, "The y_axis of the contour is zero, this should never happen. y_axis: " << y_axis << ", y_axis.norm() = " << y_axis.norm() << ".");

      if (x_axis.norm() < 2e-14)
        {
          // The origin and x_direction point are at the same loation. No easy fix, so return a Point with NaN
          // so that higher up function can properly deal with it. This basically means that the distance
          // between the two curves is the depth (90 degree dipping surface) or just zero.
          x_axis = Point<3>(NaN::DSNAN,NaN::DSNAN,NaN::DSNAN,cartesian);
        }

      x_axis = x_axis / x_axis.norm();
      y_axis = y_axis / y_axis.norm();

      return {x_axis,y_axis};
    }

    DistanceInterpolationData
    Contours::distance_interpolation_data(const Point<3> &check_point_cartesian,
                                          const Objects::NaturalCoordinate &check_point_natural,
                                          const Point<2> &reference_point,
                                          const std::unique_ptr<CoordinateSystems::Interface> &coordinate_system,
                                          const std::vector<std::vector<double> > &interpolation_properties,
                                          const double /*start_radius*/,
                                          const bool only_positive) const
    {
      // do some preparations
      const CoordinateSystem natural_coordinate_system = coordinate_system->natural_coordinate_system();
      const bool bool_cartesian = natural_coordinate_system == cartesian;

      const std::array<double,3> &check_point_surface_2d_array = check_point_natural.get_coordinates();
      const Point<3> check_point_surface(bool_cartesian ? check_point_surface_2d_array[0] : start_radius,
                                         check_point_surface_2d_array[1],
                                         bool_cartesian ? start_radius : check_point_surface_2d_array[2],
                                         natural_coordinate_system);
      const Point<2> check_point_surface_2d(check_point_natural.get_surface_coordinates(),
                                            natural_coordinate_system);


      const double check_point_depth = start_radius-check_point_natural.get_depth_coordinate();
      const size_t n_contour_curves = contour_curves.size();

      double min_distance = std::numeric_limits<double>::infinity();

      DistanceInterpolationData return_distance_interpolation_data;

      for (size_t curve_i = 0; curve_i < n_contour_curves-1; ++curve_i)
        {
          // first check whether the checkpoint is within the range of the curve
          std::cout << "-------> Curve " << curve_i << ", check_point_depth = " << check_point_depth << ", depths[curve_i] = " << depths[curve_i] << ", max_thickness_on_curve[curve_i] = " << max_thickness_on_curve[curve_i] << "(" << depths[curve_i]-max_thickness_on_curve[curve_i]  << ")"
                    << ", depths[curve_i+1] = " << depths[curve_i+1] << ", max_thickness_on_curve[curve_i+1] = " << max_thickness_on_curve[curve_i+1] << "(" << depths[curve_i+1]+max_thickness_on_curve[curve_i+1]  << ")" << std::endl;
          if (check_point_depth > depths[curve_i]-max_thickness_on_curve[curve_i] && check_point_depth < depths[curve_i+1]+max_thickness_on_curve[curve_i+1])
            {
              // potientially in range. Compute the distance for the point to the curve above and below.
              ClosestPointOnCurve closest_point_above = contour_curves[curve_i].closest_point_on_curve_segment(check_point_surface_2d);
              std::cout << "==> closest_point_above "
                        << "distance: " << closest_point_above.distance
                        << ", index: " << closest_point_above.index
                        << ", parametric_fraction: " << closest_point_above.parametric_fraction
                        << ", interpolation_fraction: " << closest_point_above.interpolation_fraction
                        << ", point: " << closest_point_above.point
                        << std::endl;
              ClosestPointOnCurve closest_point_below = contour_curves[curve_i+1].closest_point_on_curve_segment(check_point_surface_2d);
              std::cout << "==> closest_point_below "
                        << "distance: " << closest_point_below.distance
                        << ", index: " << closest_point_below.index
                        << ", parametric_fraction: " << closest_point_below.parametric_fraction
                        << ", interpolation_fraction: " << closest_point_below.interpolation_fraction
                        << ", point: " << closest_point_below.point
                        << std::endl;

              // get angle above and below and make a bezier curve out of it
              std::vector<double> local_angle_contraints =
              {
                // above ()
                (contour_curves[curve_i].angles[closest_point_above.index+1]-contour_curves[curve_i].angles[closest_point_above.index]) *closest_point_above.interpolation_fraction+contour_curves[curve_i].angles[closest_point_above.index],
                (contour_curves[curve_i+1].angles[closest_point_below.index+1]-contour_curves[curve_i+1].angles[closest_point_below.index]) *closest_point_below.interpolation_fraction+contour_curves[curve_i+1].angles[closest_point_below.index]
              };

              // get 3D location of closest point above and below
              Point<3> closest_point_above_3D_local =  bool_cartesian
                                                       ?
                                                       Point<3>(closest_point_above.point[0],closest_point_above.point[1],start_radius-depths[curve_i], cartesian)
                                                       :
                                                       Utilities::spherical_to_cartesian_coordinates(Point<3>(start_radius-depths[curve_i],closest_point_above.point[0],closest_point_above.point[1],spherical).get_array());

              Point<3> closest_point_below_3D_local =  bool_cartesian
                                                       ?
                                                       Point<3>(closest_point_below.point[0],closest_point_below.point[1],start_radius-depths[curve_i+1], cartesian)
                                                       :
                                                       Utilities::spherical_to_cartesian_coordinates(Point<3>(start_radius-depths[curve_i+1],closest_point_below.point[0],closest_point_below.point[1],spherical).get_array());

              // get 2D plane going through closest point above and below
              // Note: may or may not go through the actual point.
              //       Todo: test other options which will actually go through the point?

              // convert the 3D top and bottom point to the 2D plane coordinates (top will probably be 0,0)

              std::cout << "interpl. angle "
                        << ", above: " << local_angle_contraints[0] << "(" << contour_curves[curve_i].angles[closest_point_above.index] << ":" << contour_curves[curve_i].angles[closest_point_above.index+1] << ", " << closest_point_above.interpolation_fraction << ")"
                        << ", below: " << local_angle_contraints[1] << "(" << contour_curves[curve_i+1].angles[closest_point_below.index] << ":" << contour_curves[curve_i+1].angles[closest_point_below.index+1] << ", " << closest_point_below.interpolation_fraction<< ")"
                        << ", closest_point_above_3D_local: " << closest_point_above_3D_local << ", closest_point_below_3D_local: " << closest_point_below_3D_local
                        << std::endl;

              // todo: test this out with large d, might need to make a 2d cross section out of it...
              // on the other hand, the main goal is to make the surface smooth, which is exactly what this is doing
              // It therefore might not matter how the points are connected smoothly, as long as it is smooth and within bounderies,
              // with it should be.
              // todo: other option is to let spline go through the point, but on the other hand that assumes that the
              // closest point to the surface is directly above it, which may not always (or even most of the time) be true.
              const std::vector<Point<2> > local_points = {closest_point_above.point,closest_point_below.point};
              ClosestPointOnCurve closest_point_local = BezierCurve(local_points,local_angle_contraints).closest_point_on_curve_segment(check_point_surface_2d);
              std::cout << "==> closest_point_local "
                        << "distance: " << closest_point_local.distance
                        << ", index: " << closest_point_local.index
                        << ", parametric_fraction: " << closest_point_local.parametric_fraction
                        << ", interpolation_fraction: " << closest_point_local.interpolation_fraction
                        << ", point: " << closest_point_local.point
                        << std::endl;

              // so now we know where the "closest point on the surface" is. Compute the distance between it and the check point.
              //const double closest_point_on_surface_depth = depths[curve_i+1]*closest_point_local.interpolation_fraction+depths[curve_i]; // correct
              const double closest_point_on_surface_depth = depths[curve_i+1]*closest_point_local.parametric_fraction+depths[curve_i]; //approx
              const Point<3> closest_point_on_surface = Point<3>(closest_point_local.point[0],closest_point_local.point[1],start_radius-closest_point_on_surface_depth,cartesian);
              const double distance = (closest_point_on_surface-check_point_cartesian).norm();
              // deterimine the sign: up is negative, down is positive (kind of represents a local depth system)
              const double sign = closest_point_on_surface_depth > check_point_depth ? -1. : 1.;
              const double final_distance = only_positive && sign < 0. ? std::numeric_limits<double>::infinity() : sign*distance;
              std::cout << "interpolated depth:" << closest_point_on_surface_depth << ", closest_point_on_bezier_surface: " << closest_point_on_surface << ", distance = " << distance << ", sign = " << sign << std::endl;
              if (std::abs(final_distance < min_distance))
                {
                  // compute the distance along the path up. Approximate for now by just interpolating
                  // the distance values on the nearest points on the curve
                  return_distance_interpolation_data.distance_along_surface = 100e3;//closest_point_local.distance + (distance_along_surface[curve_i][closest_point_local.index]+distance_along_surface[curve_i][closest_point_local.index+1]*closest_point_local.interpolation_fraction);
                  return_distance_interpolation_data.signed_distance_from_bezier_surface = final_distance;

                  return_distance_interpolation_data.curve_above_index = curve_i;
                  return_distance_interpolation_data.curve_above_section_index = closest_point_above.index;
                  return_distance_interpolation_data.curve_above_interplation_fraction = closest_point_above.interpolation_fraction;
                  return_distance_interpolation_data.curve_below_index = curve_i+1;
                  return_distance_interpolation_data.curve_below_section_index = closest_point_below.index;
                  return_distance_interpolation_data.curve_below_interplation_fraction = closest_point_below.interpolation_fraction;

                  return_distance_interpolation_data.curve_local_interpolation_fraction = closest_point_local.interpolation_fraction;
                  //return_distance_interpolation_data.
                }

            }
        }

      return return_distance_interpolation_data;

    }

  } // namespace Objects
} // namespace WorldBuilder

