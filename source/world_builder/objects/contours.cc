/*
  Copyright (C) 2018 - 2021 by the authors of the World Builder code.

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
    //Contours::Contours()
    //{}

    Contours::Contours(std::vector<std::vector<Point<2> > > points,
                       std::vector<double> depths,
                       std::vector<std::vector<double> > angle_contraints,
                       std::vector<std::vector<double> > thicknesses)
      :
      points(points),
      depths(depths),
      angle_contraints(angle_contraints),
      thicknesses(thicknesses)
    {
      for (size_t curve_i = 0; curve_i < points.size(); ++curve_i)
        {
          contour_curves.emplace_back(BezierCurve(points[curve_i],angle_contraints[curve_i]));
        }
    }

    /*
        Contours::Contours(std::vector<std::pair<double,std::vector<double>>> &contours)
        {
          WBAssertThrow(contours.size() >= 2, "The countours object needs at least 2 conours.");
          // first create the connectivity
          // The connectivity is created by connecting each point with one or many points in a countour below. One poiont below can have many connections to a point above.
          // Connect each point with the closest point on the contour below.
          // Assumptions: Countours always get deeper and the points are ordered in the same direction.
          connectivity.resize(contours.size()-1);
          connectivity.resize(contours.size());
          for (unsigned int contour_i = 0; contour_i < contours.size()-1; contour_i++)
            {
              std::cout << " ---------------- " << contour_i << "--------------------- " << std::endl;

              connectivity[contour_i].resize(contours[contour_i].second.size()/2);
              // start with connecting the edge to prevent an if statement within the loop
              connectivity[contour_i][0].emplace_back(0);

              // now start connecting the rest
              double depth = contours[contour_i].first;
              double depth_next = contours[contour_i+1].first;
              unsigned int countour_entry_i = 0;
              unsigned int next_countour_entry_i = 0;
              std::cout << "start search for contour " << contour_i << ", contours[contour_i].second.size() = " << contours[contour_i].second.size()/2. << ",  contours[contour_i+1].second.size() = " <<  (contours[contour_i+1].second.size()/2.)-1 <<std::endl;
              while (countour_entry_i < contours[contour_i].second.size()/2. && next_countour_entry_i < contours[contour_i+1].second.size()/2.)
                {
                  std::cout << "countour_entry_i: " << countour_entry_i << ", next_countour_entry_i: " << next_countour_entry_i << std::endl;
                  // compute distances
                  // between the current top point and the next down point.
                  double x = contours[contour_i].second[2*countour_entry_i];
                  double y = contours[contour_i].second[2*countour_entry_i+1];
                  double x_next = contours[contour_i+1].second[2*(next_countour_entry_i+1)];
                  double y_next = contours[contour_i+1].second[2*(next_countour_entry_i+1)+1];

                  std::cout << "x:y = " << x << ":" << y << ", next: " << x_next << ":" << y_next << std::endl;

                  double distance_1 = std::sqrt((depth-depth_next)*(depth-depth_next) + (x-x_next)*(x-x_next) + (y-y_next)*(y-y_next));
                  double distance_2 = std::numeric_limits<double>::infinity();
                  double distance_3 = std::numeric_limits<double>::infinity();
                  bool countour_next_entry = false;
                  bool next_countour_next_entry = false;
                  if (contours[contour_i].second.size() > 2*(countour_entry_i+1))
                    {
                      // between the next top point and the current down ponit.
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
            }

          // second compute the dips and distances along the slab to the surface.

        }

        bool Contours::in_triangle(const std::array<std::array<double,3>,3> &points,
                                   const Point<2> check_point,
                                   double &interpolate_value) const
        {
          double factor = 20.;
          // based on https://stackoverflow.com/questions/2049582/how-to-determine-if-a-point-is-in-a-2d-triangle
          // compute s, t and area
          const double s_no_area = -(points[0][1]*points[2][0] - points[0][0]*points[2][1] + (points[2][1] - points[0][1])*check_point[0] + (points[0][0] - points[2][0])*check_point[1]);
          const double t_no_area = -(points[0][0]*points[1][1] - points[0][1]*points[1][0] + (points[0][1] - points[1][1])*check_point[0] + (points[1][0] - points[0][0])*check_point[1]);
          const double two_times_area = -(-points[1][1]*points[2][0] + points[0][1]*(-points[1][0] + points[2][0]) + points[0][0]*(points[1][1] - points[2][1]) + points[1][0]*points[2][1]);

          if (s_no_area >= -factor*std::numeric_limits<double>::epsilon() && t_no_area >= -factor*std::numeric_limits<double>::epsilon() && s_no_area+t_no_area-two_times_area<=two_times_area*factor*std::numeric_limits<double>::epsilon())
            {
              // point is in this triangle
              const double one_over_two_times_area = 1./two_times_area;
              const double s = one_over_two_times_area*s_no_area;
              const double t = one_over_two_times_area*t_no_area;
              interpolate_value = points[0][2]*(1-s-t)+points[1][2]*s+points[2][2]*t;
              return true;
            }
          return false;
        }

        double Contours::local_value(const Point<2> check_point) const
        {
          if (constant_value)
            {
              // just min and max are the same since it is constant. Just return min.
              return minimum;
            }
          // first find the closest centeroids
          const KDTree::IndexDistances index_distances = tree.find_closest_points(check_point);

          Point<2> other_point = check_point;
          KDTree::IndexDistances index_distances_other;
          bool spherical = false;
          if (check_point.get_coordinate_system() == CoordinateSystem::spherical)
            {
              spherical = true;
              other_point[0] += check_point[0] < 0 ? 2.0 * WorldBuilder::Utilities::Consts::PI : -2.0 * WorldBuilder::Utilities::Consts::PI;
              index_distances_other = tree.find_closest_points(other_point);
            }
          // try triangle of the closest centroid
          double interpolated_value = 0;

          if (in_triangle(triangles[tree.get_nodes()[index_distances.min_index].index],check_point,interpolated_value))
            {
              return interpolated_value;
            }
          else if (spherical && in_triangle(triangles[tree.get_nodes()[index_distances_other.min_index].index],other_point,interpolated_value))
            {
              return interpolated_value;
            }
          else
            {
              // if not found go to closets nodes
              // Todo: could remove the cosest node, because it was already tested. Could also sort based no distance.
              for (auto &index_distance: index_distances.vector)
                {
                  if (in_triangle(triangles[tree.get_nodes()[index_distance.index].index],check_point,interpolated_value))
                    {
                      return interpolated_value;
                    }
                  else if (spherical && in_triangle(triangles[tree.get_nodes()[index_distance.index].index],other_point,interpolated_value))
                    {
                      // This is probably non-optimal, but it seems to work better than expected
                      return interpolated_value;
                    }
                }

              // if still not found, go through all nodes
              // Todo: Although this shouldonly very rearly happen, could remove already tested nodes.
              for (const auto &nodes: tree.get_nodes())
                {
                  if (in_triangle(triangles[nodes.index],check_point,interpolated_value))
                    {
                      return interpolated_value;
                    }
                  else  if (spherical && in_triangle(triangles[nodes.index],other_point,interpolated_value))
                    {
                      return interpolated_value;
                    }
                }
              WBAssertThrow(false, "Internal error: The requested point was not in any triangle. "
                            << "This could be due to rounding errors if the difference between the check point and triangle points are small, "
                            << "or you are requesting a point ouside the bounderies defined by the additional points. The check point was "
                            << check_point[0] <<  ":" << check_point[1] << ".");
            }
          return 0;
        }*/
  } // namespace Objects
} // namespace WorldBuilder

