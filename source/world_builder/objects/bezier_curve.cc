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
#include "world_builder/objects/bezier_curve.h"
#include "world_builder/objects/cubic_spline.h"
#include "world_builder/assert.h"
#include "world_builder/nan.h"
#include "world_builder/utilities.h"

#include <sstream>
#include <cstddef>
#include <iomanip>

using namespace WorldBuilder;

namespace WorldBuilder
{
  namespace Objects
  {
    BezierCurve::BezierCurve(const std::vector<Point<2> > &p, const std::vector<double> &angle_constrains_input)
    {
      points = p;
      const size_t n_points = p.size();
      control_points.resize(n_points-1, {p[0],p[0]});
      lengths.resize(n_points-1,NaN::DSNAN);
      angles.resize(n_points,NaN::DSNAN);
      std::vector<double> angle_constrains = angle_constrains_input;
      angle_constrains.resize(n_points,NaN::DSNAN);

      // if no angle is provided, compute the angle as the average angle bewteen the previous and next point.
      // The first angle points at the second point and the last angle points at the second to last point.
      // The check points are set at a distance of 1/10th the line length from the point in the direction of the angle.
      if (std::isnan(angle_constrains[0]))
        {
          Point<2> P1P2 = points[1]-points[0];
          angles[0] = atan2(P1P2[0],P1P2[1]);
        }
      else
        {
          angles[0] = angle_constrains[0];
        }

      for (size_t p_i = 1; p_i < n_points-1; ++p_i)
        {
          bool clockwise = true;

          // first determine the angle
          if (std::isnan(angle_constrains[p_i]))
            {
              // get the average angle
              const Point<2> P1P2 = points[p_i-1]-points[p_i];
              const Point<2> P3P2 = points[p_i+1]-points[p_i];

              double diff_angle = atan2(P1P2[0]*P3P2[1] - P1P2[1]*P3P2[0] , P1P2[0]*P3P2[0] + P1P2[1]*P3P2[1]);
              const double angle_p1p2 = atan2(P1P2[1],P1P2[0]);
              const double angle_p3p1 = atan2(P3P2[1],P3P2[0]);
              const double average_angle = (angle_p1p2 + angle_p3p1)*0.5;
              angles[p_i] = average_angle;

              // if direction changes between clockwise and anti-clockwise
              // adjust direction by 90 degrees.
              //  if (diff_angle < 0 && !clockwise)
              //    {
              //      clockwise=true;
              //      angles[p_i] -= Consts::PI*0.5;
              //    }
              //  else if (diff_angle > 0 && clockwise)
              //    {
              //      clockwise=false;
              //      angles[p_i] -= Consts::PI*0.5;
              //    }
              angles[p_i] -= Consts::PI*0.5;
            }
          else
            {
              angles[p_i] = angle_constrains[p_i];
            }

        }

      if (std::isnan(angle_constrains[n_points-1]))
        {
          Point<2> P1P2 = points[n_points-1]-points[n_points-2];
          angles[n_points-1] = atan2(P1P2[0],P1P2[1]);
        }
      else
        {
          angles[n_points-1] = angle_constrains[n_points-1];
        }


      // next determine the location of the control points
      // the location of the control point is 1/10th p1p2 distance in the direction of the angle.
      // make sure the angle is pointing away from the next point, e.g.
      // the check point is on the other side of the of the line p1p2 compared to p3.
      const double fraction_of_length = 0.2;
      {
        const Point<2> &p1 = points[0];
        const Point<2> &p2 = points[1];
        const Point<2> &p3 = points[2];
        const double length = (points[0]-points[1]).norm(); // can be squared
        control_points[0][0][0] = cos(angles[0])*length*fraction_of_length+p1[0];
        control_points[0][0][1] = sin(angles[0])*length*fraction_of_length+p1[1];
        control_points[0][1][0] = cos(angles[1])*length*fraction_of_length+p2[0];
        control_points[0][1][1] = sin(angles[1])*length*fraction_of_length+p2[1];
        {
          const int side_of_line_1 =  (p1[0] - p2[0]) * (control_points[0][1][1] - p1[1])
                                      - (p1[1] - p2[1]) * (control_points[0][1][0] - p1[0])
                                      < 0 ? -1.0 : 1.0;
          const int side_of_line_2 =  (p1[0] - p2[0]) * (p3[1] - p1[1])
                                      - (p1[1] - p2[1]) * (p3[0] - p1[0])
                                      < 0 ? -1.0 : 1.0;
          if (side_of_line_1 == side_of_line_2)
            {
              // use a 180 degree rotated angle to create this control_point
              control_points[0][1][0] = cos(angles[1]+Consts::PI)*length*fraction_of_length+p2[0];
              control_points[0][1][1] = sin(angles[1]+Consts::PI)*length*fraction_of_length+p2[1];
            }
        }
        std::cout << "p1= " << points[0] << ", cp1= " << control_points[0][0] << ", cp2= " << control_points[0][1] << ", p2=" << points[1]
                  << ", angle1= " << angles[0] << ":" << angles[0]*180./Consts::PI << ", angle2= " << angles[0+1] << ":" << angles[1]*180./Consts::PI << std::endl;
      }

      for (size_t p_i = 1; p_i < n_points-1; ++p_i)
        {
          const Point<2> &p1 = points[p_i];
          const Point<2> &p2 = points[p_i+1];
          const Point<2> &p3 = points[p_i+2];
          const double length = (points[p_i]-points[p_i+1]).norm(); // can be squared
          control_points[p_i][0][0] = cos(angles[p_i])*length*fraction_of_length+p1[0];
          control_points[p_i][0][1] = sin(angles[p_i])*length*fraction_of_length+p1[1];

          {
            const int side_of_line_1 =  (p1[0] - p2[0]) * (control_points[p_i-1][1][1] - p1[1])
                                        - (p1[1] - p2[1]) * (control_points[p_i-1][1][0] - p1[0])
                                        < 0 ? -1.0 : 1.0;
            const int side_of_line_2 =  (p1[0] - p2[0]) * (control_points[p_i][0][1] - p1[1])
                                        - (p1[1] - p2[1]) * (control_points[p_i][0][0] - p1[0])
                                        < 0 ? -1.0 : 1.0;
            if (side_of_line_1 == side_of_line_2)
              {
                // use a 180 degree rotated angle to create this control_point
                control_points[p_i][0][0] = cos(angles[p_i]+Consts::PI)*length*fraction_of_length+p1[0];
                control_points[p_i][0][1] = sin(angles[p_i]+Consts::PI)*length*fraction_of_length+p1[1];
                std::cout << "switching a" << std::endl;
              }
          }

          control_points[p_i][1][0] = cos(angles[p_i+1])*length*fraction_of_length+points[p_i+1][0];
          control_points[p_i][1][1] = sin(angles[p_i+1])*length*fraction_of_length+points[p_i+1][1];

          {
            const int side_of_line_1 =  (p1[0] - p2[0]) * (control_points[p_i][1][1] - p1[1])
                                        - (p1[1] - p2[1]) * (control_points[p_i][1][0] - p1[0])
                                        < 0 ? -1.0 : 1.0;
            const int side_of_line_2 =  (p1[0] - p2[0]) * (p3[1] - p1[1])
                                        - (p1[1] - p2[1]) * (p3[0] - p1[0])
                                        < 0 ? -1.0 : 1.0;
            if (side_of_line_1 == side_of_line_2)
              {
                // use a 180 degree rotated angle to create this control_point
                control_points[p_i][1][0] = cos(angles[p_i+1]+Consts::PI)*length*fraction_of_length+p2[0];
                control_points[p_i][1][1] = sin(angles[p_i+1]+Consts::PI)*length*fraction_of_length+p2[1];
                std::cout << "switching b" << std::endl;
              }
          }

          std::cout << "p1= " << points[p_i] << ", cp1= " << control_points[p_i][0] //<< ", df=" << (points[p_i]-control_points[p_i][0]).norm()
                    << ", cp2= " << control_points[p_i][1] << ", p2=" << points[p_i+1] //<< ", df=" << (points[p_i+1]-control_points[p_i][1]).norm()
                    << ", angle1= " << angles[p_i] << ":" << angles[p_i]*180./Consts::PI << ", angle2= " << angles[p_i+1] << ":" << angles[p_i+1]*180./Consts::PI << std::endl;
        }



      /*points = p;
      const size_t n_points = p.size();
      // resize all vectors
      control_points.resize(n_points-1, p[0]);
      lengths.resize(n_points-1,NaN::DSNAN);
      angles.resize(n_points,NaN::DSNAN);
      std::vector<double> angle_constrains = angle_constrains_input;
      angle_constrains.resize(n_points,NaN::DSNAN);

      // if no angle provided, compute angles as an average angle between the previous and next point +90 degrees
      // when done, compute the angle for the first and last point using the second and second-to-last point respectively
      std::cout << "n_points = "<< n_points << std::endl;
      for (size_t p_i = 1; p_i < n_points-1; ++p_i)
        {
          std::cout << "p_i=" << p_i << ", n_points = " << n_points << std::endl;
          bool clockwise = true;
          if (std::isnan(angle_constrains[p_i]))
            {
              const Point<2> P1P2 = points[p_i-1]-points[p_i];
              const Point<2> P3P2 = points[p_i+1]-points[p_i];

              double diff_angle = atan2(P1P2[0]*P3P2[1] - P1P2[1]*P3P2[0] , P1P2[0]*P3P2[0] + P1P2[1]*P3P2[1]);
              const double angle_p1p2 = atan2(P1P2[1],P1P2[0]);
              const double angle_p3p1 = atan2(P3P2[1],P3P2[0]);
              const double average_angle = (angle_p1p2 + angle_p3p1)*0.5;
              angle_constrains[p_i] = average_angle;
              //if(diff_angle < 0 && !clockwise){
              //  clockwise=true;
              //  angle_constrains[p_i] -= Consts::PI*0.5;
              //  std::cout << "  switch 1" << std::endl;
              //} else if(diff_angle > 0 && clockwise) {
              //  clockwise=false;
              //  angle_constrains[p_i] -= Consts::PI*0.5;
              //  std::cout << "  switch 2" << std::endl;
              //}
              angle_constrains[p_i] = average_angle-Consts::PI*0.5*(p_i==2 ? 2 : 1);
              std::cout << p_i << ": angle_p1p2=" << angle_p1p2 << ", angle_p3p1=" << angle_p3p1 << ", aa=" << average_angle << ":" << average_angle*180./Consts::PI
                        << ", diff_angle = " << diff_angle << ":" << diff_angle*180./Consts::PI << (diff_angle > 0 ? " --> anti-clockwise" : " --> clockwise") << std::endl;
            }
        }

      std::cout << "======================" << std::endl;
      // now compute the first and last angle
      {
        Point<2> P1P2 = points[1]-points[0];
        double angle_p1p2 = atan2(P1P2[0],P1P2[1]);
        double average_angle = (angle_p1p2*9. + (angle_constrains[1] < -Consts::PI ? 2*Consts::PI+angle_constrains[1] : angle_constrains[1]))*0.1;//(angle_p1p2 +  (angle_constrains[1] < -Consts::PI ? 2*Consts::PI+angle_constrains[1] : angle_constrains[1]))*0.5;//(angle_p1p2 +angle_p1p2 +angle_p1p2 + angle_constrains[1])*0.25;


        const Point<2> P3P2 = points[2]-points[1];
        double diff_angle = atan2(P1P2[0]*P3P2[1] - P1P2[1]*P3P2[0] , P1P2[0]*P3P2[0] + P1P2[1]*P3P2[1]);

        angle_constrains[0] = average_angle;//-Consts::PI*0.5;//angle_constrains[1]*.5;//-Consts::PI*0.5;
        std::cout << "a angle_p1p2=" << angle_p1p2 << ":" << angle_p1p2*180./Consts::PI<< ", angle_constrains[1] = " << angle_constrains[1] << ":" << angle_constrains[1]*180./Consts::PI
                  << ", aa=" << average_angle << ":" << average_angle*180./Consts::PI
                  << ", aaa = " << (angle_p1p2 + angle_p1p2 + angle_p1p2 + (angle_constrains[1] < -Consts::PI ? 2*Consts::PI+angle_constrains[1] : angle_constrains[1]))*0.25*180./Consts::PI
                  << ", angle_constrains[1] = " << angle_constrains[1]  << ", ac=" << (angle_constrains[1] < -Consts::PI ? 2*Consts::PI+angle_constrains[1] : angle_constrains[1])<< std::endl;

        P1P2 = points[points.size()-2]-points[points.size()-1];
        angle_p1p2 = atan2(P1P2[1],P1P2[0]);
        average_angle = (angle_p1p2 +angle_p1p2 +angle_p1p2 + angle_constrains[angle_constrains.size()-2])*0.25;
        angle_constrains[angle_constrains.size()-1] = (average_angle);//-Consts::PI*0.5);
        std::cout << "b angle_p1p2=" << angle_p1p2 << ":" << angle_p1p2*180./Consts::PI<< ", angle_constrains[angle_constrains.size()-2] = " << angle_constrains[angle_constrains.size()-2] << ":" << angle_constrains[angle_constrains.size()-2]*180./Consts::PI
                  << ", aa=" << average_angle << ":" << average_angle*180./Consts::PI << ", diff_angle = " << diff_angle << ":" << diff_angle*180./Consts::PI << std::endl;
      }
      for (size_t p_i = 0; p_i < angle_constrains.size(); p_i++)
        {
          std::cout  << p_i << " ac = " << angle_constrains[p_i] << ":" << angle_constrains[p_i]*180./Consts::PI << std::endl;
        }*/


      /////////////////////////////////////////////////////////////////////////
      /*points = p;
      // first compute the factors for a monotome spline
      const size_t n_points = p.size();
      std::vector<double> x(n_points);
      std::vector<double> y(n_points);
      for (size_t i = 0; i < n_points; ++i)
        {
          x[i] = p[i][0];
          y[i] = p[i][1];
        }
      Utilities::interpolation x_spline;
      x_spline.set_points(x);
      Utilities::interpolation y_spline;
      y_spline.set_points(y);

      //// resize all vectors
      control_points.resize(n_points-1, p[0]);
      lengths.resize(n_points-1,NaN::DSNAN);
      angles.resize(n_points,NaN::DSNAN);
      std::vector<double> angle_constrains = angle_constrains_input;
      angle_constrains.resize(n_points,NaN::DSNAN);
      if (std::isnan(angle_constrains[0]))
        {
          const double &value_x = x_spline.m[0][3];
          const double &value_y = y_spline.m[0][3];
          const double value_x_h = x_spline.m[0][0]*0.5*0.5*0.5 + x_spline.m[0][1]*0.5*0.5 + x_spline.m[0][2]*0.5 + x_spline.m[0][3];
          const double value_y_h = y_spline.m[0][0]*0.5*0.5*0.5 + y_spline.m[0][1]*0.5*0.5 + y_spline.m[0][2]*0.5 + y_spline.m[0][3];
          angles[0] = std::atan2(value_y_h-value_y,value_x_h-value_x);
          std::cout << std::endl << "angle[0]:" << angles[0]*180./M_PI << ", value:" << value_x << ":" << value_y << ", value h:" << value_x_h << ":" << value_y_h << std::endl;
        }
      else
        {
          // NOTE: start angle assumes slabs or faults going down, which means they should provide a negative angle to get expected behavoir
          angles[0] = angle_constrains[0];
        }

      for (size_t i = 0; i < n_points-2; ++i)
        {
          //first compute next angle:
          if (std::isnan(angle_constrains[i+1]))
            {
              const double &value_x = x_spline.m[i+1][3];
              const double &value_y = y_spline.m[i+1][3];
              const double value_x_h = x_spline.m[i][0]*0.5*0.5*0.5 + x_spline.m[i][1]*0.5*0.5 + x_spline.m[i][2]*0.5 + x_spline.m[i][3];
              const double value_y_h = y_spline.m[i][0]*0.5*0.5*0.5 + y_spline.m[i][1]*0.5*0.5 + y_spline.m[i][2]*0.5 + y_spline.m[i][3];
              angles[i+1] = std::atan2(value_y-value_y_h,value_x-value_x_h);
              std::cout << "angle["<< i+1 << "]:" << angles[i+1]*180./M_PI << ", value:" << value_x << ":" << value_y << ", value h:" << value_x_h << ":" << value_y_h  << std::endl;
            }
          else
            {
              // NOTE: start angle doesn't assumes slabs or faults going down, which means they should provide a negative angle to get expected behavoir
              angles[i+1] = angle_constrains[i+1];
            }

          // now find the control point: where the two lines intersect:
          // p0.x
          // if angles are the same, the control point is either on the line or at infinity. Put it at P[i] for now
          if (std::abs(fmod((fmod(angles[i]-angles[i+1],180.) + 180.), 180.)) < std::numeric_limits<double>::epsilon()*10.)
            {
              control_points[i] = p[i];
              std::cout << "flag 1: "<< control_points[i] << std::endl;
            }
          else
            {
              std::cout << "flag 2" << std::endl;
              const double &x0 = p[i][0];
              const double &y0 = p[i][1];
              const double &x1 = p[i+1][0];
              const double &y1 = p[i+1][1];
              if (std::abs(fmod((fmod(angles[i], 180.) + 180.), 180.) - 90.) < std::numeric_limits<double>::epsilon()*10.)
                {
                  std::cout << "flag 3" << std::endl;
                  // vertical line at x = x0
                  control_points[i][0] = x0;
                  control_points[i][1] = std::tan(angles[i+1]) * (x0-x1) + y1;
                }
              else if (std::abs(fmod((fmod(angles[i+1], 180.) + 180.), 180.) - 90.) < std::numeric_limits<double>::epsilon()*10.)
                {
                  std::cout << "flag 4" << std::endl;
                  // vertical line at x = x0
                  control_points[i][0] = x1;
                  control_points[i][1] = std::tan(angles[i]) * (x1-x0) + y0;
                }
              std::cout << "flag 5" << std::endl;
              std::cout << ", flag 2.3: angles i: " << angles[i] << " (" << angles[i]*180./M_PI  << "), angles i+1: " << angles[i+1] << " (" << angles[i+1]*180./M_PI  << ")";
              // tan(angle) = opposite/adjecent
              const double m0 = std::tan(angles[i]); // Line 0: y = m0 (x - x0) + y0
              const double m1 = std::tan(angles[i+1]); // Line 1: y = m1 (x - x1) + y1
              // m0 (x - x0) + y0 = m1 (x - x1) + y1
              // m0*x - m0*x0 + y0 = m1*x - m1*x1 + y1
              // m0*x - m1*x =  -m1*x1 + y1 + m0*x0 - y0
              // (m0*-m1)x =  -m1*x1 + y1 + m0*x0 - y0
              // x =  (-m1*x1 + y1 + m0*x0 - y0)/(m0*-m1)
              // x =  (m0*x0 - m1*x1 + y1 - y0)/(m0*-m1)
              const double control_x = ((m0 * x0 - m1 * x1) - (y0 - y1)) / (m0 - m1);

              control_points[i][0] = control_x;
              control_points[i][1] = m0 * (control_x - x0) + y0;

              // compute length of segment
              lengths[i] = arc_length(p[i],control_points[i],p[i+1]);
              WBAssert(!std::isnan(lengths[i]),"");
            }
        }

      if (std::isnan(angle_constrains[n_points-1]))
        {
          const double value_x = x_spline.m[n_points-2][0] + x_spline.m[n_points-2][1] + x_spline.m[n_points-2][2] + x_spline.m[n_points-2][3];
          const double value_y = y_spline.m[n_points-2][0] + y_spline.m[n_points-2][1] + y_spline.m[n_points-2][2] + y_spline.m[n_points-2][3];
          // have the angle be pointing the the previous halfpoint instead of the next one
          const double value_x_h = x_spline.m[n_points-2][0]*0.5*0.5*0.5 + x_spline.m[n_points-2][1]*0.5*0.5 + x_spline.m[n_points-2][2]*0.5 + x_spline.m[n_points-2][3];
          const double value_y_h = y_spline.m[n_points-2][0]*0.5*0.5*0.5 + y_spline.m[n_points-2][1]*0.5*0.5 + y_spline.m[n_points-2][2]*0.5 + y_spline.m[n_points-2][3];
          angles[n_points-1] = std::atan2(value_y_h-value_y,value_x_h-value_x);
        }
      else
        {
          // NOTE: start angle assumes slabs or faults going down, which means they should provide a negative angle to get expected behavoir
          angles[n_points-1] = angle_constrains[n_points-1];
        }
      size_t i = n_points-2;

      if (std::abs(fmod((fmod(angles[i]-angles[i+1],180.) + 180.), 180.)) < std::numeric_limits<double>::epsilon()*10.)
        {
          control_points[i] = p[i];
        }
      else
        {
          const double &x0 = p[i][0];
          const double &y0 = p[i][1];
          const double &x1 = p[i+1][0];
          const double &y1 = p[i+1][1];
          if (std::abs(fmod((fmod(angles[i], 180.) + 180.), 180.) - 90.) < std::numeric_limits<double>::epsilon()*10.)
            {
              // vertical line at x = x0
              control_points[i][0] = x0;
              control_points[i][1] = std::tan(angles[i+1]) * (x0-x1) + y1;
            }
          else if (std::abs(fmod((fmod(angles[i+1], 180.) + 180.), 180.) - 90.) < std::numeric_limits<double>::epsilon()*10.)
            {
              // vertical line at x = x0
              control_points[i][0] = x1;
              control_points[i][1] = std::tan(angles[i]) * (x1-x0) + y0;
            }

          const double m0 = std::tan(angles[i]); // Line 0: y = m0 (x - x0) + y0
          const double m1 = std::tan(angles[i+1]); // Line 1: y = m1 (x - x1) + y1
          const double control_x = ((m0 * x0 - m1 * x1) - (y0 - y1)) / (m0 - m1);

          control_points[i][0] = control_x;
          control_points[i][1] = m0 * (control_x - x0) + y0;
        }

      // now that all the control points have been generated, check whether they make sense: The control point of line 1 need to be on a different side of the line than for line 2 if the angle between line 1 and 2 is largen than 90 degrees.
      std::cout << std::endl;
      for (size_t cp_i = 0; cp_i < control_points.size()-1; ++cp_i)
        {
          Point<2> p1 = points[cp_i];
          Point<2> p2 = points[cp_i+1];
          Point<2> p3 = points[cp_i+2];

          // check if angle between lines is more than 90 degrees
          Point<2> vector_1 = p2-p1;
          Point<2> vector_2 = p3-p2;
          double angle = atan2(vector_1[0]*vector_2[1] - vector_1[1]*vector_2[0] , vector_1[0]*vector_2[0] + vector_1[1]*vector_2[1]);
          std::cout << "=====> cp " << cp_i << " angle diff = " << angle  << ":" <<  angle*180./M_PI<< std::endl;
          if (std::abs(angle) > Consts::PI*0.5)
            {
              // check if the control point is on different sides of the line
              const int side_of_line_1 =  (p1[0] - p2[0]) * (control_points[cp_i][1] - p1[1])
                                          - (p1[1] - p2[1]) * (control_points[cp_i][0] - p1[0])
                                          < 0 ? -1.0 : 1.0;
              const int side_of_line_2 =  (p2[0] - p3[0]) * (control_points[cp_i][1] - p3[1])
                                          - (p2[1] - p3[1]) * (control_points[cp_i][0] - p3[0])
                                          < 0 ? -1.0 : 1.0;
              std::cout << "=====> side of line " << side_of_line_1 << ": " << side_of_line_2  << std::endl;

              if (side_of_line_1 == side_of_line_2)
                {
                  // mirror point
                  const double A = p2[1] - p1[1];
                  const double B = -(p2[0] - p1[0]);
                  const double C = -A * p1[0] - B * p1[1];
                  const double M = sqrt(A * A + B * B);
                  const double A_norm = A / M;
                  const double B_norm = B / M;
                  const double C_norm = C / M;
                  const double D = A_norm * control_points[cp_i][0] + B_norm * control_points[cp_i][1] + C_norm;

                  Point<2> mirror_point = Point<2>(control_points[cp_i][0] - 2 * A_norm * D, control_points[cp_i][1]  - 2 * B_norm * D, cartesian);
                  std::cout << "control_points[cp_i]=" << control_points[cp_i] << ", mirror_point=" << mirror_point << std::endl;
                  control_points[cp_i] = mirror_point;

                  // update angle
                  //const double &value_x = x_spline.m[cp_i][3];
                  //const double &value_y = y_spline.m[cp_i][3];
                  //const double value_x_h = mirror_point[0];
                  //const double value_y_h = mirror_point[1];
                  //angles[cp_i] = std::atan2(value_y_h-value_y,value_x_h-value_x);
                }
            }
        }

      // compute length of segment
      lengths[n_points-2] = arc_length(p[n_points-2],control_points[n_points-2],p[n_points-1]);

      for (size_t i = 0; i < n_points-1; ++i)
        std::cout << "============> " << i << ": p1=" << p[i] << ", control=" << control_points[i] << ", p2=" << p[i+1] << std::endl;

      WBAssert(!std::isnan(lengths[n_points-2]),
               "lengths[n_points-2] = " << lengths[n_points-2] << ", n_points:" << n_points << ", p[n_points-2]: " << p[n_points-2]
               << ", control_points[n_points-2]: "<< control_points[n_points-2]<< ", p[n_points-1]: " << p[n_points-1]);*/
    }


    Point<2>
    BezierCurve::operator()(const size_t i, const double t) const
    {
      return (1-t)*(1-t)*(1-t)*points[i] + 3*(1-t)*(1-t)*t*control_points[i][0] + 3.*(1-t)*t*t*control_points[i][1]+t*t*t*points[i+1]; //(1-t)*(1-t)*points[i] + 2*t*(1-t)*control_points[i] + t*t*points[i+1];
    }

    double
    BezierCurve::arc_length(const Point<2> &a, const Point<2>  &b, const Point<2> &c) const
    {
      // compute the arc length of the Bezier curve
      // see https://gamedev.stackexchange.com/questions/6009/bezier-curve-arc-length
      Point<2> v = a;
      Point<2> w = b;
      v[0] = 2.*(b[0] - a[0]);
      v[1] = 2.*(b[1] - a[1]);
      w[0] = c[0] - 2.*b[0] + a[0];
      w[1] = c[1] - 2.*b[1] + a[1];

      const double uu = 4.*(w[0]*w[0] + w[1]*w[1]);

      // check if all points are in a line
      if (uu < 0.00001 || (a-b).norm() < 1e-6 || (b-c).norm() < 1e-6 || fabs(a[0] * (b[1] - c[1]) + b[0] * (c[1] - a[1]) + c[0] * (a[1] - b[1])) <= 1e-9)
        {
          return std::sqrt((c[0] - a[0])*(c[0] - a[0]) + (c[1] - a[1])*(c[1] - a[1]));
        }

      const double vv = 4.*(v[0]*w[0] + v[1]*w[1]);
      const double ww = v[0]*v[0] + v[1]*v[1];

      const double t1 = (2.*std::sqrt(uu*(uu + vv + ww)));
      const double t2 = 2.*uu+vv;
      const double t3 = vv*vv - 4.*uu*ww;
      const double t4 = (2.*std::sqrt(uu*ww));

      WBAssert(!std::isnan(((t1*t2 - t3*std::log(t2+t1) -(vv*t4 - t3*std::log(vv+t4))) / (8*std::pow(uu, 1.5)))), "result is nan: " << ((t1*t2 - t3*std::log(t2+t1) -(vv*t4 - t3*std::log(vv+t4))) / (8*std::pow(uu, 1.5)))
               << ", t1=" << t1 << ", t2=" << t2 << ", t3=" << t3 << ", t4=" << t4 << ", v = " << v << ", w=" << w << ", vv=" << vv << ", uu=" << uu
               << ", a=" << a << ", b=" << b << ", c=" << c
               << ", std::log(t2+t1)=" << std::log(t2+t1) << ", std::log(vv+t4)=" << std::log(vv+t4) << ", std::pow(uu, 1.5)=" << std::pow(uu, 1.5) << ", std::log(vv+t4)=" << std::log(vv+t4));

      return ((t1*t2 - t3*std::log(t2+t1) -(vv*t4 - t3*std::log(vv+t4))) / (8*std::pow(uu, 1.5)));
    }


    double
    BezierCurve::arc_length(const size_t index, const double fraction) const
    {
      /*// This method uses a similair approach as https://malczak.info/blog/quadratic-bezier-curve-length
      // but instead of the full length, we integrate the full equaion (https://www.wolframalpha.com/input?i=integrate+sqrt%28A*t%5E2%2BB*t%2Bc%29+dt)
      // leaving t in there. Then we compute the integration constant by stating that the length at t=0 should
      // be zero.
      auto dt = points[index]-2.*control_points[index]+points[index+1];
      auto ddt = 2.*control_points[index]-2.*points[index];
      const double a = 4*(dt[0]*dt[0]+dt[1]*dt[1]);
      const double c = ddt[0]*ddt[0]+ddt[1]*ddt[1];

      if (a < 5e-4 * c || points[index] == control_points[index] || control_points[index] == points[index+1])
        {
          // all points are on a line
          return sqrt(((points[index+1]-points[index])*fraction)*((points[index+1]-points[index])*fraction));//std::sqrt((points[index][0] + dx1*t)*(points[index][0] + dx1*t) + (points[index][1] + dy1*t)*(points[index][1] + dy1*t));
        }

      const double b = 4*(dt[0]*ddt[0]+dt[1]*ddt[1]);

      const double inv_4a = 1./(4*a);
      const double u = (b*b)*inv_4a;
      const double k = (4*a*c-b*b)*inv_4a;
      const double sqrt_a = sqrt(a);
      const double sqrt_c = sqrt(c);
      const double inv_8_sqrt_a_a_a = 1./(8.*sqrt(a*a*a));
      const double sqrt_c_t_b_at = sqrt(c+fraction*(b+a*fraction));
      double x = fraction*sqrt_a+sqrt(u);

      // todo: optimize
      const double integral = ((b+2.*a*fraction)*sqrt_c_t_b_at)*inv_4a - ((b*b-4.*a*c)*log(b+2.*a*fraction+2.*sqrt_a*sqrt_c_t_b_at))*inv_8_sqrt_a_a_a;
      const double constant = (b*sqrt_c)*inv_4a - ((b*b-4.*a*c)*log(b+2.*sqrt_a*sqrt_c))*inv_8_sqrt_a_a_a;

      return integral-constant;*/
      return 0;
    }





    //ClosestPointOnCurve
    //BezierCurve::closest_point_on_curve_segment(const Point<2> &check_point) const
    //{

    /*//if (std::fabs(check_point[0]+31500) < 1e-1 && std::fabs(check_point[1]-243750)<1e-1) //std::fabs(check_point[0]+198000) < 1e-1 && std::fabs(check_point[1]-425000)<1e-1 || std::fabs(check_point[0]+315000) < 1e-1 && std::fabs(check_point[1]-431250)<1e-1)
    if (std::fabs(check_point[0]-31500) < 1e-1 && std::fabs(check_point[1]-256250)<1e-1)
      std::cout << "start closest_point_on_curve_segment" << std::endl;

    // go through each section and find all roots in domain 0 to 1 and choose the smallest
    ClosestPointOnCurve closest_point_on_curve;


    double min_estimate_solution = -1.;
    double minimum_distance_to_reference_point = std::numeric_limits<double>::infinity();

    Point<2> P1 = points[0];
    Point<2> P2 = P1;

    const size_t number_of_points = points.size();

    for ( size_t i_estimate = 0; i_estimate < control_points.size(); ++i_estimate)
      {
        //double real_roots[4] = {NaN::DSNAN,NaN::DSNAN,NaN::DSNAN,0.};
        //if (std::fabs(check_point[0]+31500) < 1e-1 && std::fabs(check_point[1]-243750)<1e-1) //std::fabs(check_point[0]+198000) < 1e-1 && std::fabs(check_point[1]-425000)<1e-1 || std::fabs(check_point[0]+315000) < 1e-1 && std::fabs(check_point[1]-431250)<1e-1)
        if (std::fabs(check_point[0]-31500) < 1e-1 && std::fabs(check_point[1]-256250)<1e-1)
          std::cout << "start i_estimate " << i_estimate << std::endl;

        // Go one point pair up.
        P1 = P2;
        P2 = points[i_estimate+1];

        // Compute distance to first point on the line.
        if (P1 == P2)
          {
            continue;
          }

        double minimum_distance_to_reference_point_tmp = P1.cheap_relative_distance_cartesian(check_point);
        {
          const Point<2> A = control_points[i_estimate]-points[i_estimate];
          const Point<2> B = points[i_estimate+1]-control_points[i_estimate]-A;
          const double a = B*B;
          const double b = 3*A*B;
          const double c = 2.*A*A+(points[i_estimate]-check_point)*B;
          const double d = (points[i_estimate]-check_point)*A;
          //const double d = (points[i_estimate]-check_point)*A;
          const double fraction_2 = 0;
    const double fraction_4 = 0;
            const double function_value = a*fraction_2*fraction_2*fraction_2 + b*fraction_2*fraction_2 + c*fraction_2+d;
            const double derivative_value = (3.*a*fraction_2*fraction_2 + 2.*b*fraction_2 + c);
            const double second_derivative_value = 6*a*fraction_2 + 2*b;
          minimum_distance_to_reference_point_tmp = function_value;



            if (minimum_distance_to_reference_point_tmp < minimum_distance_to_reference_point)
              {
                //if (std::abs(newton_root-real_roots[root_i]) > 1e-5)
                //  std::cout << "newton_root = " << newton_root << ", real_roots = " << real_roots[root_i] << ", diff = " << newton_root-real_roots[root_i] << std::endl;
                // the sign is rotating the derivative by 90 degrees.
                // When moving in the direction of increasing t, left is positve and right is negative.
                const double &t = fraction_4;
                const Point<2> derivative_point = ((-2+2*t)*points[i_estimate] + (2-4*t)*control_points[i_estimate] + 2*t*points[i_estimate+1]);
                const Point<2> second_derivative_point = ((2)*points[i_estimate] + (-4)*control_points[i_estimate] + 2*points[i_estimate+1]);
                Point<2> tangent_point = derivative_point - P1;
                // if angle between check point and tangent point is larget than 90 degrees, return a negative distance
                const double dot_product = (tangent_point*(check_point-P1));
                const double sign = dot_product < 0. ? -1. : 1.;
                tangent_point = Point<2>(-tangent_point[1],tangent_point[0],tangent_point.get_coordinate_system());

                closest_point_on_curve.distance = sign*P1.distance(check_point);
                closest_point_on_curve.parametric_fraction = fraction_4;
                closest_point_on_curve.interpolation_fraction = NaN::DSNAN; //arc_length(i,real_roots[root_i])/lengths[i];
                closest_point_on_curve.index = i_estimate;
                closest_point_on_curve.point = P1;

                if (std::fabs(check_point[0]-31500) < 1e-1 && std::fabs(check_point[1]-256250)<1e-1)
                  std::cout << "set new distance!" << std::endl;
              }
        }
        //if (minimum_distance_to_reference_point_tmp < minimum_distance_to_reference_point)
        //  {
        //    minimum_distance_to_reference_point = minimum_distance_to_reference_point_tmp;
        //    min_estimate_solution = i_estimate;
        //  }

        // Compute fraction of a straight line between the two points where the check point is closest.
        Point<2> P1P2 = (P2)-(P1);
        Point<2> P1Pc = check_point-(P1);

        double c_1 = P1Pc*P1P2;


        //if (std::fabs(check_point[0]+31500) < 1e-1 && std::fabs(check_point[1]-243750)<1e-1) //std::fabs(check_point[0]+198000) < 1e-1 && std::fabs(check_point[1]-425000)<1e-1 || std::fabs(check_point[0]+315000) < 1e-1 && std::fabs(check_point[1]-431250)<1e-1)
        if (std::fabs(check_point[0]-31500) < 1e-1 && std::fabs(check_point[1]-256250)<1e-1)
          std::cout << "c_1 =  " << c_1<< std::endl;


        if ( c_1 < 0 )
          {
            // closest point to segment before P1. Continue
            // preventing a rounded corner
            min_estimate_solution = i_estimate == 0 ? -1 : min_estimate_solution;
            continue;
          }
        double c_2 = P1P2*P1P2;


        //if (std::fabs(check_point[0]+31500) < 1e-1 && std::fabs(check_point[1]-243750)<1e-1) //std::fabs(check_point[0]+198000) < 1e-1 && std::fabs(check_point[1]-425000)<1e-1 || std::fabs(check_point[0]+315000) < 1e-1 && std::fabs(check_point[1]-431250)<1e-1)
        if (std::fabs(check_point[0]-31500) < 1e-1 && std::fabs(check_point[1]-256250)<1e-1)
          std::cout << "c_1 =  " << c_1 << ", c_2 = " << c_2 << ", c_1 / c_2 = " << c_1 / c_2 << std::endl;

        if ( c_2 < c_1 )
          {
            // closest point to segment after P2. Continue
            continue;
          }

        double fraction = c_1 / c_2;

        //if(fraction > 1.5)
        //{
        //  continue;
        //}


        //if (std::fabs(check_point[0]+31500) < 1e-1 && std::fabs(check_point[1]-243750)<1e-1) //std::fabs(check_point[0]+198000) < 1e-1 && std::fabs(check_point[1]-425000)<1e-1 || std::fabs(check_point[0]+315000) < 1e-1 && std::fabs(check_point[1]-431250)<1e-1)
        if (std::fabs(check_point[0]-31500) < 1e-1 && std::fabs(check_point[1]-256250)<1e-1)
          std::cout << "fraction =  " << fraction << std::endl;


        // Compute distance of check point to the spline using the computed fraction from the straight line.
        double min_estimate_solution_tmp = (i_estimate+fraction);
        WBAssert(min_estimate_solution_tmp>=0 && min_estimate_solution_tmp <=number_of_points, "message");

        const double index_1 = (size_t)min_estimate_solution_tmp;
        const double fraction_1 = min_estimate_solution_tmp - index_1;
        const Point<2> A = control_points[index_1]-points[index_1];
        const Point<2> B = points[index_1+1]-control_points[index_1]-A;
        const double a = B*B;
        const double b = 3*A*B;
        const double c = 2.*A*A+(points[index_1]-check_point)*B;
        const double d = (points[index_1]-check_point)*A;
        const double minimum_distance_to_reference_point_start = a*fraction_1*fraction_1*fraction_1 + b*fraction_1*fraction_1 + c*fraction_1+d;//((*this)(index_1,fraction_1)-check_point).norm();

        double new_distance_tmp = -1;
        size_t i_newton_iteration = 0;
        // Newton iteration (modified to use the sign of the derivative instead of sign of derivative/second dervative)
        // to compute the point on the spine which is closest to the check point. The modification allows for a larger
        // area of convergence since it will only converge to minima, and not to maxima.
        while (true)
          {

            const double index_2 = (size_t)min_estimate_solution_tmp;
            const double fraction_2 = min_estimate_solution_tmp - index_2;
            // compute a,b,c and d in the cubic equation describing the distance from point p to the local quadratic Bezier curve.
            // using https://blog.gludion.com/2009/08/distance-to-quadratic-bezier-curve.html
            // todo: I should also take a look at: https://iquilezles.org/articles/distfunctions2d/
            const Point<2> A = control_points[index_2]-points[index_2];
            const Point<2> B = points[index_2+1]-control_points[index_2]-A;
            const double a = B*B;
            const double b = 3*A*B;
            const double c = 2.*A*A+(points[index_2]-check_point)*B;
            const double d = (points[index_2]-check_point)*A;


            const double function_value = a*fraction_2*fraction_2*fraction_2 + b*fraction_2*fraction_2 + c*fraction_2+d;
            const double derivative_value = (3.*a*fraction_2*fraction_2 + 2.*b*fraction_2 + c);
            const double second_derivative_value = 6*a*fraction_2 + 2*b;
            if (std::fabs(derivative_value) < 1e-8)
              {
                minimum_distance_to_reference_point_tmp = function_value;//((*this)(index_2,fraction_2)-check_point).norm();
                //if (std::fabs(check_point[0]+31500) < 1e-1 && std::fabs(check_point[1]-243750)<1e-1) //std::fabs(check_point[0]+198000) < 1e-1 && std::fabs(check_point[1]-425000)<1e-1 || std::fabs(check_point[0]+315000) < 1e-1 && std::fabs(check_point[1]-431250)<1e-1)
                if (std::fabs(check_point[0]-31500) < 1e-1 && std::fabs(check_point[1]-256250)<1e-1)
                  std::cout << "break because derivatvie zero" << std::endl;
                break;
              }

            // We take the Newton derivative between some limits and only the sign of the
            // derivative, not the derivative. This ensures that we converge to the
            // min value and not the max value.
            const double update = std::min(0.25,std::max(-0.25,function_value/(derivative_value)));//derivative_value/(second_derivative_value);//std::min(0.25,std::max(-0.25,derivative_value/std::abs(second_derivative_value)));


            if (std::fabs(update) < 1e-4)
              {
                minimum_distance_to_reference_point_tmp = function_value;//((*this)(index_1,fraction_1)-check_point).norm();
                WBAssertThrow(std::fabs(minimum_distance_to_reference_point_tmp) < std::fabs(minimum_distance_to_reference_point_start),
                              "Failed to converge on spline. Initial guess " << std::setprecision(16) << minimum_distance_to_reference_point_start
                              << ") smaller than final result(" << minimum_distance_to_reference_point_tmp << ", diff = " << minimum_distance_to_reference_point_start-minimum_distance_to_reference_point_tmp
                              << ") for point " << check_point << ", cp = " << check_point << ", min_estimate_solution_tmp = " << min_estimate_solution_tmp << ", update = " << update << ", i_newton_iteration = " << i_newton_iteration
                              //<< ", new dist expensive = " << Point<2>(x_spline(min_estimate_solution_tmp),y_spline(min_estimate_solution_tmp),cartesian).cheap_relative_distance_cartesian(check_point_surface_2d)
                              <<  ".");
                break;
              }

            //minimum_distance_to_reference_point_tmp = ((*this)(index_2,fraction_2)-check_point).norm();
            double update_scaling = 1;
            // only the first few iterations need some guidence from line search
            bool line_search_finished = false;
            if (std::fabs(update) > 1e-4)
              {
                const Point<2> A = control_points[index_2]-points[index_2];
                const Point<2> B = points[index_2+1]-control_points[index_2]-A;
                const double a = B*B;
                const double b = 3*A*B;
                const double c = 2.*A*A+(points[index_2]-check_point)*B;
                const double d = (points[index_2]-check_point)*A;
                minimum_distance_to_reference_point_tmp = a*fraction_2*fraction_2*fraction_2 + b*fraction_2*fraction_2 + c*fraction_2+d;//((*this)(index_2,fraction_2)-check_point).norm();//a*fraction_2*fraction_2*fraction_2 + b*fraction_2*fraction_2 + c*fraction_2+d;

                //if (std::fabs(check_point[0]+31500) < 1e-1 && std::fabs(check_point[1]-243750)<1e-1) //std::fabs(check_point[0]+198000) < 1e-1 && std::fabs(check_point[1]-425000)<1e-1 || std::fabs(check_point[0]+315000) < 1e-1 && std::fabs(check_point[1]-431250)<1e-1)
                if (std::fabs(check_point[0]-31500) < 1e-1 && std::fabs(check_point[1]-256250)<1e-1)
                  {
                    const Point<2> A = control_points[0]-points[0];
                    const Point<2> B = points[0+1]-control_points[0]-A;
                    const double a = B*B;
                    const double b = 3*A*B;
                    const double c = 2.*A*A+(points[0]-check_point)*B;
                    const double d = (points[0]-check_point)*A;
                    std::cout << "exact solution 1: " <<  a*1.*1.*1. + b*1.*1. + c*1.+d  << ", derivative = " << (3.*a*1.*1. + 2.*b*1. + c) << ", sd = " <<  6.*a*1. + 2.*b << std::endl;//((*this)(index_2,fraction_2)-check_point).norm();//a*fraction_2*fraction_2*fraction_2 + b*fraction_2*fraction_2 + c*fraction_2+d;
                    std::cout << a << "*x^3+" << b << "*z^2+" << c << "*x+" << d << std::endl;
                    std::cout << "points = " << points[0] << ", " << points[1] << ", " << control_points[0] << std::endl;
                  }


                //if (std::fabs(check_point[0]+31500) < 1e-1 && std::fabs(check_point[1]-243750)<1e-1) //std::fabs(check_point[0]+198000) < 1e-1 && std::fabs(check_point[1]-425000)<1e-1 || std::fabs(check_point[0]+315000) < 1e-1 && std::fabs(check_point[1]-431250)<1e-1)
                if (std::fabs(check_point[0]-31500) < 1e-1 && std::fabs(check_point[1]-256250)<1e-1)
                  {
                    const Point<2> A = control_points[1]-points[1];
                    const Point<2> B = points[1+1]-control_points[1]-A;
                    const double a = B*B;
                    const double b = 3*A*B;
                    const double c = 2.*A*A+(points[1]-check_point)*B;
                    const double d = (points[1]-check_point)*A;
                    std::cout << "exact solution 2: " <<  d  << ", derivative = " << c << ", sd = " <<  b << std::endl;//((*this)(index_2,fraction_2)-check_point).norm();//a*fraction_2*fraction_2*fraction_2 + b*fraction_2*fraction_2 + c*fraction_2+d;
                    std::cout << a << "*x^3+" << b << "*z^2+" << c << "*x+" << d << std::endl;
                    std::cout << "points = " << points[1] << ", " << points[2] << ", " << control_points[1] << std::endl;
                  }

                //if (std::fabs(check_point[0]+31500) < 1e-1 && std::fabs(check_point[1]-243750)<1e-1) //std::fabs(check_point[0]+198000) < 1e-1 && std::fabs(check_point[1]-425000)<1e-1 || std::fabs(check_point[0]+315000) < 1e-1 && std::fabs(check_point[1]-431250)<1e-1)
                if (std::fabs(check_point[0]-31500) < 1e-1 && std::fabs(check_point[1]-256250)<1e-1)
                  std::cout << "abcd: " << a *fraction_2 *fraction_2 *fraction_2 + b *fraction_2 *fraction_2 + c *fraction_2+d << ", this: " << ((*this)(index_2,fraction_2)-check_point).norm() << std::endl;
                // Do a line search
                for (unsigned int i_line_search = 0; i_line_search < 25; ++i_line_search)
                  {
                    const double test_x = min_estimate_solution_tmp - update_scaling*update;
                    //if (std::fabs(check_point[0]+31500) < 1e-1 && std::fabs(check_point[1]-243750)<1e-1) //std::fabs(check_point[0]+198000) < 1e-1 && std::fabs(check_point[1]-425000)<1e-1 || std::fabs(check_point[0]+315000) < 1e-1 && std::fabs(check_point[1]-431250)<1e-1)
                    if (std::fabs(check_point[0]-31500) < 1e-1 && std::fabs(check_point[1]-256250)<1e-1)
                      std::cout << "testx = " << test_x << ", min_estimate_solution_tmp = " << min_estimate_solution_tmp << ", update_scaling = " << update_scaling << ", update = " << update << ", fraction_2 = " << fraction_2 << std::endl;
                    {
                      const double index_3 = (size_t)test_x;
                      const double fraction_3 = test_x - index_3;

                      //if (std::fabs(check_point[0]+31500) < 1e-1 && std::fabs(check_point[1]-243750)<1e-1) //std::fabs(check_point[0]+198000) < 1e-1 && std::fabs(check_point[1]-425000)<1e-1 || std::fabs(check_point[0]+315000) < 1e-1 && std::fabs(check_point[1]-431250)<1e-1)
                      if (std::fabs(check_point[0]-31500) < 1e-1 && std::fabs(check_point[1]-256250)<1e-1)
                        std::cout << "index_3: " << index_3 << ", fraction_3:" << fraction_3 << std::endl;
                      // compute a,b,c and d in the cubic equation describing the distance from point p to the local quadratic Bezier curve.
                      // using https://blog.gludion.com/2009/08/distance-to-quadratic-bezier-curve.html
                      // todo: I should also take a look at: https://iquilezles.org/articles/distfunctions2d/
                      const Point<2> A = control_points[index_3]-points[index_3];
                      const Point<2> B = points[index_3+1]-control_points[index_3]-A;
                      const double a = B*B;
                      const double b = 3*A*B;
                      const double c = 2.*A*A+(points[index_3]-check_point)*B;
                      const double d = (points[index_3]-check_point)*A;
                      //if (std::fabs(check_point[0]+31500) < 1e-1 && std::fabs(check_point[1]-243750)<1e-1) //std::fabs(check_point[0]+198000) < 1e-1 && std::fabs(check_point[1]-425000)<1e-1 || std::fabs(check_point[0]+315000) < 1e-1 && std::fabs(check_point[1]-431250)<1e-1)
                      if (std::fabs(check_point[0]-31500) < 1e-1 && std::fabs(check_point[1]-256250)<1e-1)
                        std::cout << a << "*x^3+" << b << "*z^2+" << c << "*x+" << d << std::endl;

                      //const size_t index_2 = 0;
                      //const Point<2> A = control_points[index_2]-points[index_2];
                      //const Point<2> B = points[index_2+1]-control_points[index_2]-A;
                      //const double a = B*B;
                      //const double b = 3*A*B;
                      //const double c = 2.*A*A+(points[index_2]-check_point)*B;
                      //const double d = (points[index_2]-check_point)*A;
                      //const double fraction_2=1.;
                      const double function_value = a*fraction_3*fraction_3*fraction_3 + b*fraction_3*fraction_3 + c*fraction_3+d;
                      const double derivative_value = (3.*a*fraction_3*fraction_3 + 2.*b*fraction_3 + c);
                      const double second_derivative_value = 6*a*fraction_3 + 2*b;

                      const double minimum_distance_to_reference_point_line_search = a*fraction_3*fraction_3*fraction_3 + b*fraction_3*fraction_3 + c*fraction_3+d;//((*this)(index_3,fraction_3)-check_point).norm();
                      //if (std::fabs(check_point[0]+31500) < 1e-1 && std::fabs(check_point[1]-243750)<1e-1) //std::fabs(check_point[0]+198000) < 1e-1 && std::fabs(check_point[1]-425000)<1e-1 || std::fabs(check_point[0]+315000) < 1e-1 && std::fabs(check_point[1]-431250)<1e-1)
                      if (std::fabs(check_point[0]-31500) < 1e-1 && std::fabs(check_point[1]-256250)<1e-1)
                        std::cout << "  i_line_search = " << i_line_search << ", mdrp_point_line_search = " << minimum_distance_to_reference_point_line_search
                                  << ", mdrp_start = " << minimum_distance_to_reference_point_start << ", mdrp_tmp = " << minimum_distance_to_reference_point_tmp << std::endl
                                  //<< "-- update_scaling = " << update_scaling << std::setprecision(8)<< ", update = " << update << ", min_estimate_solution_tmp = " << min_estimate_solution_tmp << ", testx = " << test_x
                                  //<< ", derivative_value = " << derivative_value << ", function_value = " << function_value << ", second_derivative_value = " << second_derivative_value
                                  //<< ", (this-check_point).norm()=" << ((*this)(1,0)-check_point).norm() << std::endl
                                  //<< "-- index_2=" << index_2 << ", fraction_2=" << fraction_2 << ", index_3 =" << index_3 << ", fraction_3=" << fraction_3
                                  //<< ", a=" << a << ", b=" << b << ", c=" << c << ", d=" << d << ", (*this)(index_3,fraction_3) = " << (*this)(index_3,fraction_3)
                                  << std::endl;

                      if (std::fabs(minimum_distance_to_reference_point_line_search)<std::fabs(minimum_distance_to_reference_point_start)
                          && std::fabs(minimum_distance_to_reference_point_line_search) < std::fabs(minimum_distance_to_reference_point_tmp))
                        {
                          line_search_finished = true;
                          break;
                        }

                    }
                    update_scaling*=2./3.;

                  }

                //if (std::fabs(check_point[0]+31500) < 1e-1 && std::fabs(check_point[1]-243750)<1e-1) //std::fabs(check_point[0]+198000) < 1e-1 && std::fabs(check_point[1]-425000)<1e-1 || std::fabs(check_point[0]+315000) < 1e-1 && std::fabs(check_point[1]-431250)<1e-1)
                if (std::fabs(check_point[0]-31500) < 1e-1 && std::fabs(check_point[1]-256250)<1e-1)
                  std::cout << "line search finished!" << std::endl;

                if (!line_search_finished)
                  {
                    // could not find a resonable solution. We are aparently not close, so just move on
                    //if (std::fabs(check_point[0]+31500) < 1e-1 && std::fabs(check_point[1]-243750)<1e-1) //std::fabs(check_point[0]+198000) < 1e-1 && std::fabs(check_point[1]-425000)<1e-1 || std::fabs(check_point[0]+315000) < 1e-1 && std::fabs(check_point[1]-431250)<1e-1)
                    if (std::fabs(check_point[0]-31500) < 1e-1 && std::fabs(check_point[1]-256250)<1e-1)
                      std::cout << "but could not find a solution..." << std::endl;
                    break;
                  }

              }

            min_estimate_solution_tmp = min_estimate_solution_tmp - update_scaling*update;

            if (min_estimate_solution_tmp < 0 || min_estimate_solution_tmp > number_of_points-1)
              {
                break;
              }

            //if (std::fabs(check_point[0]+31500) < 1e-1 && std::fabs(check_point[1]-243750)<1e-1) //std::fabs(check_point[0]+198000) < 1e-1 && std::fabs(check_point[1]-425000)<1e-1 || std::fabs(check_point[0]+315000) < 1e-1 && std::fabs(check_point[1]-431250)<1e-1)
            if (std::fabs(check_point[0]-31500) < 1e-1 && std::fabs(check_point[1]-256250)<1e-1)
              std::cout << i_newton_iteration <<  " (" << i_estimate << "): min_estimate_solution_tmp = " << min_estimate_solution_tmp
                        << ", update = " << update << ", minimum_distance_to_reference_point_tmp = " << minimum_distance_to_reference_point_tmp
                        << ", minimum_distance_to_reference_point = " << minimum_distance_to_reference_point
                        << ", fv = " << function_value << ", dv = " << derivative_value << std::endl;
            ++i_newton_iteration;

            WBAssertThrow(i_newton_iteration < 35+4*number_of_points,
                          "The spline solver doesn't seem to have finished on a reasonable ammount of Newton "
                          << "iterations ("<< 35+4*number_of_points << "). Please check whether your coordinates are resonable, "
                          << "or contact the maintainers. Newton interations = " << i_newton_iteration
                          << ", min_estimate_solution_tmp = " << min_estimate_solution_tmp
                          << ", last update_scaling = " << update_scaling << ", last update = " << update
                          << ", P1 = " << P1 << ", P2 = " << P2 << ", i_estimate = " << i_estimate
                          << ", minimum_distance_to_reference_point_tmp = " << minimum_distance_to_reference_point_tmp
                          << ", cp = " << check_point << ", norm = " << check_point.norm() << ", cpol = " << (*this)(index_2,fraction_2)
                          << ".");
          }

        if (minimum_distance_to_reference_point_tmp < minimum_distance_to_reference_point)
          {
            //if (std::fabs(check_point[0]+31500) < 1e-1 && std::fabs(check_point[1]-243750)<1e-1) //std::fabs(check_point[0]+198000) < 1e-1 && std::fabs(check_point[1]-425000)<1e-1 || std::fabs(check_point[0]+315000) < 1e-1 && std::fabs(check_point[1]-431250)<1e-1)
            if (std::fabs(check_point[0]-31500) < 1e-1 && std::fabs(check_point[1]-256250)<1e-1)
              std::cout << "update minimum_distance_to_reference_point!" << std::endl;
            minimum_distance_to_reference_point = minimum_distance_to_reference_point_tmp;
            min_estimate_solution = min_estimate_solution_tmp;


            const double index_4 = (size_t)min_estimate_solution;
            const double fraction_4 = min_estimate_solution - index_4;
            const Point<2> point_on_curve = (*this)(index_4,fraction_4);

            if (point_on_curve.distance(check_point) < std::abs(closest_point_on_curve.distance) && min_estimate_solution > 0 && min_estimate_solution < points.size()-1)
              {
                //if (std::abs(newton_root-real_roots[root_i]) > 1e-5)
                //  std::cout << "newton_root = " << newton_root << ", real_roots = " << real_roots[root_i] << ", diff = " << newton_root-real_roots[root_i] << std::endl;
                // the sign is rotating the derivative by 90 degrees.
                // When moving in the direction of increasing t, left is positve and right is negative.
                const double &t = fraction_4;
                const Point<2> derivative_point = ((-2+2*t)*points[index_4] + (2-4*t)*control_points[index_4] + 2*t*points[index_4+1]);
                const Point<2> second_derivative_point = ((2)*points[index_4] + (-4)*control_points[index_4] + 2*points[index_4+1]);
                Point<2> tangent_point = derivative_point - point_on_curve;
                // if angle between check point and tangent point is larget than 90 degrees, return a negative distance
                const double dot_product = (tangent_point*(check_point-point_on_curve));
                const double sign = dot_product < 0. ? -1. : 1.;
                tangent_point = Point<2>(-tangent_point[1],tangent_point[0],tangent_point.get_coordinate_system());

                closest_point_on_curve.distance = sign*point_on_curve.distance(check_point);
                closest_point_on_curve.parametric_fraction = fraction_4;
                closest_point_on_curve.interpolation_fraction = NaN::DSNAN; //arc_length(i,real_roots[root_i])/lengths[i];
                closest_point_on_curve.index = index_4;
                closest_point_on_curve.point = point_on_curve;

                if (std::fabs(check_point[0]-31500) < 1e-1 && std::fabs(check_point[1]-256250)<1e-1)
                  std::cout << "set new distance!" << std::endl;
              }
          }*/

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



    // first try to solve it with a Newton iteration
    // at^3+bt^2+ct+d = 0 and t >= 0. && t <= 1.
    // x_{n+1} = x_n - f(x)/f'(x)
    // f(x) = ax^3+bx^2+cx+d
    // f'(x) = 3ax^2+2bx+c

    //real_roots[0] = 0.5;
    //std::cout << "update 1 = " << (a*real_roots[0]*real_roots[0]*real_roots[0] + b*real_roots[0]*real_roots[0] + c*real_roots[0]+d)
    ///(3.*a*real_roots[0]*real_roots[0] + 2.*b*real_roots[0] + c)<< std::endl;
    /*for (size_t iteration_i = 0; iteration_i < 200; iteration_i++)
    {
    const double update = (a*real_roots[0]*real_roots[0]*real_roots[0] + b*real_roots[0]*real_roots[0] + c*real_roots[0]+d)
    /(3.*a*real_roots[0]*real_roots[0] + 2.*b*real_roots[0] + c);
    real_roots[0] = real_roots[0] - update;

    if (std::abs(update) < 1e-16)
    {
    //if(real_roots[0] > 0. && real_roots[0] < 1.)
    //  std::cout << "update size = " << update << ", iteration_i = " << iteration_i << "newton root = " << real_roots[0] << std::endl;
        break;
      }
      }*/

    /*double root_start_zero = 0.25;
    bool newton_failed = false;
    for (size_t iteration_i = 0; iteration_i < 200; iteration_i++)
    {
    double update = 0;
    //const double update = (a*root_start_zero*root_start_zero*root_start_zero + b*root_start_zero*root_start_zero + c*root_start_zero+d)
    //                      /(3.*a*root_start_zero*root_start_zero + 2.*b*root_start_zero + c);

    update = (root_start_zero*(root_start_zero*(a*root_start_zero + b) + c)+d)/(root_start_zero*(3.*a*root_start_zero + 2.*b) + c);
    root_start_zero = root_start_zero - update;

    //const double update_2 = (a*root_start_zero*root_start_zero*root_start_zero + b*root_start_zero*root_start_zero + c*root_start_zero+d)
    //                      /(3.*a*root_start_zero*root_start_zero + 2.*b*root_start_zero + c);
    //root_start_zero = root_start_zero - update_1;

    //if(iteration_i > 10 && update > 1e-3){
    //  newton_failed = true;
    //  break;
    //  }
    //std::cout << "0:
    update size = " << update_1 << ", iteration_i = " << iteration_i << ", root = " << root_start_zero << std::endl;
        if (std::abs(update) < 1e-12)// || root_start_zero < -0.25 || root_start_zero > 1.25)
        {
        //if(root_start_zero > 0. && root_start_zero < 1.)
        //if(iteration_i > 100)
      //  std::cout << "0:
        update size = " << update_1 << ", iteration_i = " << iteration_i << ", root = " << root_start_zero << std::endl;
            break;
          }
            if(iteration_i > 199){
            newton_failed = true;
            break;
          }
          }
            double root_start_half = 0.5;
            if(!newton_failed)
            for (size_t iteration_i = 0; iteration_i < 200; iteration_i++)
            {

            double update = 0;
            for (size_t iteration_j = 0; iteration_j < 3; iteration_j++){
            update = (root_start_half*(root_start_half*(a*root_start_half + b) + c)+d)/(root_start_half*(3.*a*root_start_half + 2.*b) + c);
            root_start_half = root_start_half - update;
          }

            //if(iteration_i > 2 && update > 1e-3){
            //  newton_failed = true;
            //  break;
            //  }
            if (std::abs(update) < 1e-12)// || root_start_half < -1.5 || root_start_half > 2.5)
            {
            //if(root_start_half > 0. && root_start_half < 1.)
            //if(iteration_i > 100)
          //  std::cout << "0.5:
            update size = " << update << ", iteration_i = " << iteration_i << ", root = " << root_start_half << std::endl;
                break;
              }
                if(iteration_i > 199){
                newton_failed = true;
                break;
              }
              }
                double root_start_one = 0.75;

                if(!newton_failed){
                for (size_t iteration_i = 0; iteration_i < 200; iteration_i++)
                {
                double update = 0;
                //const double update = (a*root_start_one*root_start_one*root_start_one + b*root_start_one*root_start_one + c*root_start_one+d)
                //                      /(3.*a*root_start_one*root_start_one + 2.*b*root_start_one + c);

                update = (root_start_zero*(root_start_zero*(a*root_start_zero + b) + c)+d)/(root_start_zero*(3.*a*root_start_zero + 2.*b) + c);
                root_start_one = root_start_one - update;

                //if(iteration_i > 10 && update > 1e-3){
                //  newton_failed = true;
                //  break;
                //  }
                if (std::abs(update) < 1e-12)// || root_start_one < -0.25 || root_start_one > 1.25)
                {
                //if(root_start_one > 0. && root_start_one < 1.)
                //if(iteration_i > 100)
              //  std::cout << "1:
                update size = " << update << ", iteration_i = " << iteration_i  << ", root = " << root_start_one << std::endl;
                    break;
                  }
                    if(iteration_i > 199){
                    newton_failed = true;
                    break;
                  }
                  }

                  //std::cout << root_start_zero << ":" << root_start_half << "
                  :" << root_start_one << "
                    , real_roots[0] = " << real_roots[0] << std::endl;

                        //const double newton_root = real_roots[0];
                        real_roots[0] = root_start_zero;
                        real_roots[1] = root_start_half;
                        real_roots[2] = root_start_one;
                        real_roots[3] = 3;
                      }
                        else {

                        // if the Newton iteration didn't work, fall back on a more robust method.
                        this->solve_cubic_equation_real(a,b,c,d,real_roots);
                      }*/

    //this->solve_cubic_equation_real(a,b,c,d,real_roots);
    /*bool newton_sucess = false;
    for (double starting_guess = 0; starting_guess <= 1.0; starting_guess+=0.2)
    {
    double root = starting_guess;
    for (size_t iteration_i = 0; iteration_i < 5; iteration_i++)
    {

    double update = 0;
    for (size_t iteration_j = 0; iteration_j < 3; iteration_j++)
    {
    update = (root*(root*(a*root + b) + c)+d)/(root*(3.*a*root + 2.*b) + c);
    root = root - update;
    }

    if (std::abs(update) < 1e-6 && root >= 0. && root <= 1)
    {

    const Point<2> point_on_curve = (*this)(i,root);

    if (point_on_curve.distance(check_point) < std::abs(closest_point_on_curve.distance))
    {
    //if (std::abs(newton_root-real_roots[root_i]) > 1e-5)
    //  std::cout << "newton_root = " << newton_root << ", real_roots = " << real_roots[root_i] << ", diff = " << newton_root-real_roots[root_i] << std::endl;
        // the sign is rotating the derivative by 90 degrees.
        // When moving in the direction of increasing t, left is positve and right is negative.
        const double &t = root;
        const Point<2> derivative_point = ((-2+2*t)*points[i] + (2-4*t)*control_points[i] + 2*t*points[i+1]);
        const Point<2> second_derivative_point = ((2)*points[i] + (-4)*control_points[i] + 2*points[i+1]);
        Point<2> tangent_point = derivative_point - point_on_curve;
        // if angle between check point and tangent point is larget than 90 degrees, return a negative distance
        const double dot_product = (tangent_point*(check_point-point_on_curve));
        const double sign = dot_product < 0. ? -1. : 1.;
        tangent_point = Point<2>(-tangent_point[1],tangent_point[0],tangent_point.get_coordinate_system());

        closest_point_on_curve.distance = sign*point_on_curve.distance(check_point);
        closest_point_on_curve.parametric_fraction = root;
        closest_point_on_curve.interpolation_fraction = NaN::DSNAN; //arc_length(i,real_roots[root_i])/lengths[i];
        closest_point_on_curve.index = i;
        closest_point_on_curve.point = point_on_curve;
      }

        newton_sucess = true;
        break;
      }
      }
      }*/
    //if(!newton_sucess)
    /*{
    this->solve_cubic_equation_real(a,b,c,d,real_roots);
    for (size_t root_i = 0; root_i < real_roots[3]; ++root_i)
    {
    if (real_roots[root_i] >= 0. && real_roots[root_i] <= 1.)
    {

    const Point<2> point_on_curve = (*this)(i_estimate,real_roots[root_i]);

    if (point_on_curve.distance(check_point) < std::abs(closest_point_on_curve.distance))
    {
    //if (std::abs(newton_root-real_roots[root_i]) > 1e-5)
    //  std::cout << "newton_root = " << newton_root << ", real_roots = " << real_roots[root_i] << ", diff = " << newton_root-real_roots[root_i] << std::endl;
        // the sign is rotating the derivative by 90 degrees.
        // When moving in the direction of increasing t, left is positve and right is negative.
        const double &t = real_roots[root_i];
        const Point<2> derivative_point = ((-2+2*t)*points[i_estimate] + (2-4*t)*control_points[i_estimate] + 2*t*points[i_estimate+1]);
        const Point<2> second_derivative_point = ((2)*points[i_estimate] + (-4)*control_points[i_estimate] + 2*points[i_estimate+1]);
        Point<2> tangent_point = derivative_point - point_on_curve;
        // if angle between check point and tangent point is larget than 90 degrees, return a negative distance
        const double dot_product = (tangent_point*(check_point-point_on_curve));
        const double sign = dot_product < 0. ? -1. : 1.;
        tangent_point = Point<2>(-tangent_point[1],tangent_point[0],tangent_point.get_coordinate_system());

        closest_point_on_curve.distance = sign*point_on_curve.distance(check_point);
        closest_point_on_curve.parametric_fraction = real_roots[root_i];
        closest_point_on_curve.interpolation_fraction = NaN::DSNAN; //arc_length(i,real_roots[root_i])/lengths[i];
        closest_point_on_curve.index = i_estimate;
        closest_point_on_curve.point = point_on_curve;
      }
      }
      }
      }*/
    /*Point<2> other_point(-40500,168750,cartesian);
    if(check_point == other_point){
    std::cout << "=========>>> hhelloo!!! i_estimate = " << i_estimate << ", cp=" << check_point << std::endl;
      }
        // Go one point pair up.
        P1 = P2;
        P2 = points[i_estimate+1];

        // Compute distance to first point on the line.
        if (P1 == P2)
        {
        continue;
      }
        if(check_point == other_point){
      std::cout << "===> flag 1 " << P1 << ":" << P2 << std::endl;
      }
        double minimum_distance_to_reference_point_tmp = P1.cheap_relative_distance_cartesian(check_point);
        if (minimum_distance_to_reference_point_tmp < minimum_distance_to_reference_point)
        {
        minimum_distance_to_reference_point = minimum_distance_to_reference_point_tmp;
        min_estimate_solution = i_estimate;
      }

        // Compute fraction of a straight line between the two points where the check point is closest.
        Point<2> P1P2 = (P2)-(P1);
        Point<2> P1Pc = check_point-(P1);

        double c_1 = P1Pc*P1P2;
        //if ( c_1 < 0 )
        //  {
        //    // closest point to segment beore P1. Continue
        //    // preventing a rounded corner
        //    min_estimate_solution = i_estimate == 0 ? -1 : min_estimate_solution;
        //    continue;
        //  }

        double c_2 = P1P2*P1P2;
        //if ( c_2 < c_1 )
        //  {
        //    // closest point to segment after P2. Continue
        //    continue;
        //  }

        double fraction = c_1 / c_2;

        if(check_point == other_point){
        std::cout << "
      ===> flag 2 fraction: " << fraction << std::endl;
      }
        double root = fraction;
        bool newton_sucess = false;
        for (size_t iteration_i = 0; iteration_i < 50; iteration_i++)
        {

        double update = 0;
        for (size_t iteration_j = 0; iteration_j < 3; iteration_j++)
        {
        update = (root*(root*(a*root + b) + c)+d)/(root*(3.*a*root + 2.*b) + c);
        root = root - update;
      }

        if (std::abs(update) < 1e-6)// && root >= 0. && root <= 1)
        {

        const Point<2> point_on_curve = (*this)(i_estimate,root);

        if (point_on_curve.distance(check_point) < std::abs(closest_point_on_curve.distance))
        {
        //if (std::abs(newton_root-real_roots[root_i]) > 1e-5)
        //  std::cout << "
        newton_root = " << newton_root << ", real_roots = " << real_roots[root_i] << ", diff = " << newton_root-real_roots[root_i] << std::endl;
            // the sign is rotating the derivative by 90 degrees.
            // When moving in the direction of increasing t, left is positve and right is negative.
            const double &t = root;
            const Point<2> derivative_point = ((-2+2*t)*points[i_estimate] + (2-4*t)*control_points[i_estimate] + 2*t*points[i_estimate+1]);
            const Point<2> second_derivative_point = ((2)*points[i_estimate] + (-4)*control_points[i_estimate] + 2*points[i_estimate+1]);
            Point<2> tangent_point = derivative_point - point_on_curve;
            // if angle between check point and tangent point is larget than 90 degrees, return a negative distance
            const double dot_product = (tangent_point*(check_point-point_on_curve));
            const double sign = dot_product < 0. ? -1. : 1.;
            tangent_point = Point<2>(-tangent_point[1],tangent_point[0],tangent_point.get_coordinate_system());

            closest_point_on_curve.distance = sign*point_on_curve.distance(check_point);
            closest_point_on_curve.parametric_fraction = root;
            closest_point_on_curve.interpolation_fraction = NaN::DSNAN; //arc_length(i,real_roots[root_i])/lengths[i_estimate];
            closest_point_on_curve.index = i_estimate;
            closest_point_on_curve.point = point_on_curve;
          }

            newton_sucess = true;
            break;
          }
          }*/

    //}


    //sdouble minimum_distance_to_reference_point_tmp = P2.cheap_relative_distance_cartesian(check_point);
    //sif (minimum_distance_to_reference_point_tmp < minimum_distance_to_reference_point)
    //s  {
    //s        std::cout << "update minimum_distance_to_reference_point!" << std::endl;
    //s    minimum_distance_to_reference_point = minimum_distance_to_reference_point_tmp;
    //s    min_estimate_solution = number_of_points-1;
    //s  }


    //std::cout << "return closest point index = " << closest_point_on_curve.index << ", fraction "<<  closest_point_on_curve.parametric_fraction << ", points.size()=" << points.size() << std::endl;

    //  return closest_point_on_curve;
    //}


    ClosestPointOnCurve
    BezierCurve::closest_point_on_curve_segment(const Point<2> &check_point) const
    {
      ClosestPointOnCurve closest_point_on_curve;
      const Point<2> &p = check_point;
      double min_squared_distance = std::numeric_limits<double>::infinity();
      for ( size_t cp_i = 0; cp_i < control_points.size(); ++cp_i)
        {
          const Point<2> &p1 = points[cp_i];
          const Point<2> &p2 = points[cp_i+1];
          // first check if one of the points or control points is closer than the current closest point
          {
            //if (!((check_point-p1).norm_square() <= min_squared_distance
            //      && (check_point-p2).norm_square() <= min_squared_distance
            //      && (check_point-control_points[cp_i][0]).norm_square() <= min_squared_distance
            //      && (check_point-control_points[cp_i][1]).norm_square() <= min_squared_distance))
            //  {
            //    min_squared_distance = std::min(std::min(min_squared_distance,(check_point-p1).norm_square()),(check_point-p1).norm_square());
            //    continue;
            //  }

            min_squared_distance = std::min(std::min(min_squared_distance,(check_point-p1).norm_square()),(check_point-p1).norm_square());
          }
          // Todo: check if any of the points or control points is closer than the current estimate.
          // Getting an estimate for where the closest point is with a linear approximation
          Point<2> P1P2 = p2-p1;
          Point<2> P1Pc = check_point-p1;

          double est =  std::min(1.,std::max(0.,(P1Pc*P1P2) / (P1P2*P1P2))); // est=estimate of solution
          bool found = false;
          std::stringstream output;
          //output << "cp_i=" << cp_i << ", init est = " << est << ", min_squared_distance = " << min_squared_distance << std::endl;
          Point<2> a = 3.*control_points[cp_i][0]-3.*control_points[cp_i][1]+points[cp_i+1]-points[cp_i];
          Point<2> b = 3.*points[cp_i] - 6.*control_points[cp_i][0]+3.*control_points[cp_i][1];
          Point<2> c = -3.*points[cp_i] + 3.*control_points[cp_i][0];
          Point<2> d = points[cp_i];

          //output << "wlfrm= (" << a[0] << "*x^3+" << b[0] << "*x^2+"<< c[0] << "*x+" << d[0] << "-" << p[0] << ")^2+(" << a[1] << "*x^3+" << b[1] << "*x^2+"<< c[1] << "*x+" << d[1] << "-" << p[1] << ")^2" << std::endl;
          double min_squared_distance_cartesian_temp = (a[0]*est*est*est+b[0]*est*est+c[0]*est+d[0]-p[0])*(a[0]*est*est*est+b[0]*est*est+c[0]*est+d[0]-p[0])
                                                       +(a[1]*est*est*est+b[1]*est*est+c[1]*est+d[1]-p[1])*(a[1]*est*est*est+b[1]*est*est+c[1]*est+d[1]-p[1]);
          for (size_t newton_i = 0; newton_i < 50; newton_i++)
            {
              // based on https://stackoverflow.com/questions/2742610/closest-point-on-a-cubic-bezier-curve
              a = 3.*control_points[cp_i][0]-3.*control_points[cp_i][1]+points[cp_i+1]-points[cp_i];
              b = 3.*points[cp_i] - 6.*control_points[cp_i][0]+3.*control_points[cp_i][1];
              c = -3.*points[cp_i] + 3.*control_points[cp_i][0];
              d = points[cp_i];
              //output << "  wlfrm= (" << a[0] << "*x^3+" << b[0] << "*x^2+"<< c[0] << "*x+" << d[0] << "-" << p[0] << ")^2+(" << a[1] << "*x^3+" << b[1] << "*x^2+"<< c[1] << "*x+" << d[1] << "-" << p[1] << ")^2 with x=" << est << std::endl;
              const double squared_distance_cartesian = (a[0]*est*est*est+b[0]*est*est+c[0]*est+d[0]-p[0])*(a[0]*est*est*est+b[0]*est*est+c[0]*est+d[0]-p[0])
                                                        +(a[1]*est*est*est+b[1]*est*est+c[1]*est+d[1]-p[1])*(a[1]*est*est*est+b[1]*est*est+c[1]*est+d[1]-p[1]);

              const double squared_distance_cartesian_derivative = 2.0*(3.0*a[0]*est*est+2.0*b[0]*est+c[0])*(a[0]*est*est*est+b[0]*est*est+c[0]*est+d[0]-p[0])
                                                                   + 2.0*(3.0*a[1]*est*est+2.0*b[1]*est+c[1])*(a[1]*est*est*est+b[1]*est*est+c[1]*est+d[1]-p[1]);
              const double squared_distance_cartesian_second_derivative  = 2.0*(6.0*a[0]*est + 2.0*b[0])*(a[0]*est*est*est+b[0]*est*est+c[0]*est+d[0]-p[0])
                                                                           + 2.0*(3.0*a[0]*est*est + 2.0*b[0]*est + c[0])*(3.0*a[0]*est*est + 2.0*b[0]*est + c[0])
                                                                           + 2.0*(6.0*a[1]*est + 2.0*b[1])*(a[1]*est*est*est+b[1]*est*est+c[1]*est+d[1]-p[1])
                                                                           + 2.0*(3.0*a[1]*est*est + 2.0*b[1]*est + c[1])*(3.0*a[1]*est*est + 2.0*b[1]*est + c[1]) ;
              //output << " ---->> squared_distance_cartesian = " << squared_distance_cartesian << ", squared_distance_cartesian_derivative=" << squared_distance_cartesian_derivative << ", squared_distance_cartesian_second_derivative= " << squared_distance_cartesian_second_derivative << ", est=" << est<< std::endl;
              // the local minimum is where  squared_distance_cartesian_derivative=0 and squared_distance_cartesian_derivative>=0

              const double update = squared_distance_cartesian_derivative/std::fabs(squared_distance_cartesian_second_derivative);//std::min(0.25,std::max(0.25,squared_distance_cartesian_derivative/std::fabs(squared_distance_cartesian_second_derivative)));
              double line_search = 1.;
              double est_test = est-update*line_search;
              double squared_distance_cartesian_test = (a[0]*est_test*est_test*est_test+b[0]*est_test*est_test+c[0]*est_test+d[0]-p[0])*(a[0]*est_test*est_test*est_test+b[0]*est_test*est_test+c[0]*est_test+d[0]-p[0])
                                                       +(a[1]*est_test*est_test*est_test+b[1]*est_test*est_test+c[1]*est_test+d[1]-p[1])*(a[1]*est_test*est_test*est_test+b[1]*est_test*est_test+c[1]*est_test+d[1]-p[1]);
              //output << " ->> init squared_distance_cartesian = " << squared_distance_cartesian_test << ", min_squared_distance_cartesian_temp = " << min_squared_distance_cartesian_temp << ", diff = " << squared_distance_cartesian_test-min_squared_distance_cartesian_temp << ", est = " << est << ", est_test = " << est_test << ", update: " << update << ", ls = " << line_search << ", up*ls: " << update*line_search << std::endl;

              for (unsigned int i = 0; i < 10; i++)
                {
                  est_test = est-update*line_search;

                  squared_distance_cartesian_test = (a[0]*est_test*est_test*est_test+b[0]*est_test*est_test+c[0]*est_test+d[0]-p[0])*(a[0]*est_test*est_test*est_test+b[0]*est_test*est_test+c[0]*est_test+d[0]-p[0])
                                                    +(a[1]*est_test*est_test*est_test+b[1]*est_test*est_test+c[1]*est_test+d[1]-p[1])*(a[1]*est_test*est_test*est_test+b[1]*est_test*est_test+c[1]*est_test+d[1]-p[1]);
                  //output << "    ->> test squared_distance_cartesian = " << squared_distance_cartesian_test << ", ls =" << line_search << ", est_test=" << est_test<< std::endl;
                  //double squared_distance_cartesian_derivative_test = 2.0*(3.0* a[0]*est_test*est_test + 2.0*b[0]*est_test + c[0])*(a[0]*est_test*est_test*est_test+b[0]*est_test*est_test+c[0]*est_test+d[0]-p[0])
                  //                                               + 2.0*(a[1]*est_test*est_test+b[1]*est_test+c[1])*(a[1]*est_test*est_test*est_test+b[1]*est_test*est_test+c[1]*est_test+d[1]-p[1]);

                  //double squared_distance_cartesian_second_derivative_test  = 2.0*(6.0*a[0]*est_test + 2.0*b[0])*(a[0]*est_test*est_test*est_test+b[0]*est_test*est_test+c[0]*est_test+d[0]-p[0])
                  //                                                       + 2.0*(3.0*a[0]*est_test*est_test + 2.0*b[0]*est_test + c[0])*(3.0*a[0]*est_test*est_test + 2.0*b[0]*est_test + c[0])
                  //                                                       + 2.0*(6.0*a[1]*est_test + 2.0*b[1])*(a[1]*est_test*est_test*est_test+b[1]*est_test*est_test+c[1]*est_test+d[1]-p[1])
                  //                                                       + 2.0*(3.0*a[1]*est_test*est_test + 2.0*b[1]*est_test + c[1])*(3.0*a[1]*est_test*est_test + 2.0*b[1]*est_test + c[1]) ;
                  //output << "    ls: " << newton_i<< ": i=" << cp_i << ", est_test=" << est_test << ", update=" << update*0.5 << ", deriv=" << squared_distance_cartesian_derivative_test << ", sec_deriv=" << squared_distance_cartesian_second_derivative_test << ", p1=" << p1 << ", p2= " << p2 << ", poc= " << a *est_test *est_test *est_test+b *est_test *est_test+c *est_test+d << ", cp= " <<  check_point << ", ds:" << ((a*est_test*est_test*est_test+b*est_test*est_test+c*est_test+d)-check_point).norm_square() << ":" << min_squared_distance_cartesian_temp << ", diff = " << squared_distance_cartesian_test-min_squared_distance_cartesian_temp<< std::endl;
                  if (squared_distance_cartesian_test < min_squared_distance_cartesian_temp)
                    break;

                  line_search *= 2./3.;
                }

              est -= update*line_search;

              min_squared_distance_cartesian_temp =  (a[0]*est*est*est+b[0]*est*est+c[0]*est+d[0]-p[0])*(a[0]*est*est*est+b[0]*est*est+c[0]*est+d[0]-p[0])
                                                     +(a[1]*est*est*est+b[1]*est*est+c[1]*est+d[1]-p[1])*(a[1]*est*est*est+b[1]*est*est+c[1]*est+d[1]-p[1]);

              //output << ", " << newton_i<< ": i=" << cp_i << ", est=" << est << ", update=" << update << ", deriv=" << squared_distance_cartesian_derivative << ", sec_deriv=" << squared_distance_cartesian_second_derivative << ", p1=" << p1 << ", p2= " << p2 << ", poc= " << a *est *est *est+b *est *est+c *est+d << ", cp= " <<  check_point << ", ds:" << ((a*est*est*est+b*est*est+c*est+d)-check_point).norm_square() << ":" << min_squared_distance_cartesian_temp<< std::endl;
              //output += "" +std::to_string( newton_i)+ ": i=" +std::to_string( cp_i )+", est=" +std::to_string( est )+", update=" +std::to_string( update )+", deriv="+std::to_string( squared_distance_cartesian_derivative )+", sec_deriv=" +std::to_string( squared_distance_cartesian_second_derivative )+ ", p1=" +std::to_string( P1 )+", p2= "+std::to_string( P2 )+", poc= " +std::to_string( a*est*est*est+b*est*est+c*est+d )+", cp= " +std::to_string(  check_point ) + std::endl;


              if (std::fabs(update) < 1e-4)
                {
                  found = true;
                  if (min_squared_distance_cartesian_temp < min_squared_distance)
                    {
                      if (est >= -1e-8 && est-1. <= 1e-8)
                        {
                          min_squared_distance = squared_distance_cartesian;
                          const Point<2> point_on_curve = a*est*est*est+b*est*est+c*est+d;

                          // the sign is rotating the derivative by 90 degrees.
                          // When moving in the direction of increasing t, left is positve and right is negative.
                          //const double &t = est;
                          //const Point<2> point_on_curve = a*est*est*est+b*est*est+c*est+d;
                          //(1-t)*(1-t)*(1-t)*points[cp_i] + 3*(1-t)*(1-t)*t *control_points[cp_i][0] + 3.*(1-t)*t *t *control_points[cp_i][1]+t *t *t *points[cp_i+1];
                          //const Point<2> derivative_point = ((-2+2*est)*points[cp_i] + (2-4*est)*control_points[cp_i] + 2*est*points[cp_i+1]);
                          //const Point<2> second_derivative_point = ((2)*points[cp_i] + (-4)*control_points[cp_i] + 2*points[cp_i+1]);
                          // https://www.wolframalpha.com/input?i=d%2Fdt+%281-t%29*%281-t%29*%281-t%29*a+%2B+3*%281-t%29*%281-t%29*t*b+%2B+3.*%281-t%29*t*t*c%2Bt*t*t*d
                          const Point<2> derivative_point = points[cp_i]*((6.-3.*est)*est-3.) + control_points[cp_i][0]*(est*(9*est-12)+3)
                                                            + control_points[cp_i][1]*(6.-9.*est)*est + points[cp_i+1]*3.*est*est;
                          //https://www.wolframalpha.com/input?i=d%2Fdt+d%2Fdt+%281-t%29*%281-t%29*%281-t%29*a+%2B+3*%281-t%29*%281-t%29*t*b+%2B+3.*%281-t%29*t*t*c%2Bt*t*t*d
                          const Point<2> second_derivative_point = points[cp_i]*(6-6*est) + control_points[cp_i][0]*(18*est-12)
                                                                   + control_points[cp_i][1]*-18*est+ points[cp_i+1]*6*est;

                          Point<2> tangent_point = derivative_point - point_on_curve;
                          // if angle between check point and tangent point is larget than 90 degrees, return a negative distance
                          const double dot_product = (tangent_point*(check_point-point_on_curve));
                          const double sign = dot_product < 0. ? -1. : 1.;
                          tangent_point = Point<2>(-tangent_point[1],tangent_point[0],tangent_point.get_coordinate_system());

                          closest_point_on_curve.distance = sign*std::sqrt(min_squared_distance);
                          closest_point_on_curve.parametric_fraction = est;
                          closest_point_on_curve.interpolation_fraction = NaN::DSNAN; //arc_length(i,real_roots[root_i])/lengths[i];
                          closest_point_on_curve.index = cp_i;
                          closest_point_on_curve.point = point_on_curve;
                          Point<2> normal = point_on_curve;
                          {
                            Point<2> derivative = a*est*est+b*est+c;
                            normal=derivative;
                            double normal_size = derivative.norm();
                            normal[0] = derivative[1]/normal_size;
                            normal[1] = -derivative[0]/normal_size;
                          }
                          closest_point_on_curve.normal = normal;
                        }
                    }
                  break;
                }
            }
          //if (std::fabs(check_point[0]-22500) < 1e-1 && std::fabs(check_point[1]-90625) < 1e-1)
          //if (std::fabs(check_point[0]-85500) < 1e-1 && std::fabs(check_point[1]-81250) < 1e-1)
          //if (std::fabs(check_point[0]-139500) < 1e-1 && std::fabs(check_point[1]-56250) < 1e-1)
          //if (std::fabs(check_point[0]-101250) < 1e-1 && std::fabs(check_point[1]-35156.3) < 1e-1)
          if (std::fabs(check_point[0]-147375) < 1e-1 && std::fabs(check_point[1]-71875) < 1e-1)
            {
              std::cout << "cp= " << check_point << ", point= " << closest_point_on_curve.point << ", fraction=" << closest_point_on_curve.parametric_fraction << ", index= " << closest_point_on_curve.index << ", normal = " << closest_point_on_curve.normal << std::endl << output.str();
            }
          WBAssertThrow(found, "Could not find a good solution. " << output.str());
        }
      if (std::fabs(check_point[0]-147375) < 1e-1 && std::fabs(check_point[1]-71875) < 1e-1)
        {
          std::cout << "--> cp= " << check_point << ", point= " << closest_point_on_curve.point << ", fraction=" << closest_point_on_curve.parametric_fraction << ", index= " << closest_point_on_curve.index << ", normal = " << closest_point_on_curve.normal << std::endl;
        }
      return closest_point_on_curve;
    }

    /*ClosestPointOnCurve
    BezierCurve::closest_point_on_curve_segment(const Point<2> &check_point) const
    {
      // go through each section and find all roots in domain 0 to 1 and choose the smallest
      ClosestPointOnCurve closest_point_on_curve;
      if (std::fabs(check_point[0]-31500) < 1e-1 && std::fabs(check_point[1]-256250)<1e-1)
        std::cout << "in point" << std::endl;
      for ( size_t i = 0; i < control_points.size(); ++i)
        {
          double real_roots[4] = {NaN::DSNAN,NaN::DSNAN,NaN::DSNAN,0.};
          // compute a,b,c and d in the cubic equation describing the distance from point p to the local quadratic Bezier curve.
          // using https://blog.gludion.com/2009/08/distance-to-quadratic-bezier-curve.html
          // todo: I should also take a look at: https://iquilezles.org/articles/distfunctions2d/
          const Point<2> A = control_points[i]-points[i];
          const Point<2> B = points[i+1]-control_points[i]-A;
          const double a = B*B;
          const double b = 3*A*B;
          const double c = 2.*A*A+(points[i]-check_point)*B;
          const double d = (points[i]-check_point)*A;
          this->solve_cubic_equation_real(a,b,c,d,real_roots);
          if (std::fabs(check_point[0]-31500) < 1e-1 && std::fabs(check_point[1]-256250)<1e-1)
            std::cout << "roots (" << real_roots[3] << "): ";

          for (size_t root_i = 0; root_i < real_roots[3]; ++root_i)
            {
              if (std::fabs(check_point[0]-31500) < 1e-1 && std::fabs(check_point[1]-256250)<1e-1)
                std::cout << root_i  << "=" << real_roots[root_i] << ", ";
              if (real_roots[root_i] >= 0. && real_roots[root_i] <= 1.)
                {

                  const Point<2> point_on_curve = (*this)(i,real_roots[root_i]);

                  if (point_on_curve.distance(check_point) < std::abs(closest_point_on_curve.distance))
                    {
                      // the sign is rotating the derivative by 90 degrees.
                      // When moving in the direction of increasing t, left is positve and right is negative.
                      const double &t = real_roots[root_i];
                      const Point<2> derivative_point = ((-2+2*t)*points[i] + (2-4*t)*control_points[i] + 2*t*points[i+1]);
                      const Point<2> second_derivative_point = ((2)*points[i] + (-4)*control_points[i] + 2*points[i+1]);
                      Point<2> tangent_point = derivative_point - point_on_curve;
                      // if angle between check point and tangent point is larget than 90 degrees, return a negative distance
                      const double dot_product = (tangent_point*(check_point-point_on_curve));
                      const double sign = dot_product < 0. ? -1. : 1.;
                      tangent_point = Point<2>(-tangent_point[1],tangent_point[0],tangent_point.get_coordinate_system());

                      closest_point_on_curve.distance = sign*point_on_curve.distance(check_point);
                      closest_point_on_curve.parametric_fraction = real_roots[root_i];
                      closest_point_on_curve.interpolation_fraction = NaN::DSNAN; //arc_length(i,real_roots[root_i])/lengths[i];
                      closest_point_on_curve.index = i;
                      closest_point_on_curve.point = point_on_curve;
                    }
                }
            }
          if (std::fabs(check_point[0]-31500) < 1e-1 && std::fabs(check_point[1]-256250)<1e-1)
            std::cout << std::endl;
        }
      return closest_point_on_curve;
    }*/

    ClosestPointOnCurve
    BezierCurve::closest_point_on_curve(const Point<2> &check_point) const
    {
      ClosestPointOnCurve closest_point_on_curve = this->closest_point_on_curve_segment(check_point);
      if (std::isnan(closest_point_on_curve.point[0]))
        {
          // The closest point is one of the edges, find which one
          const double distance_to_first_point = (points[0]-check_point).norm();
          const double distance_to_last_point = (points[points.size()-1]-check_point).norm();
          if (distance_to_first_point < distance_to_last_point)
            {
              closest_point_on_curve.distance = NaN::DSNAN; // the distance computed above is not following the curve. Could be computed if needed.
              closest_point_on_curve.index = 0;
              closest_point_on_curve.parametric_fraction = 0;
              closest_point_on_curve.interpolation_fraction = 0;
              closest_point_on_curve.point = points[0];
            }
          else
            {
              closest_point_on_curve.distance = NaN::DSNAN; // the distance computed above is not following the curve. Could be computed if needed.
              closest_point_on_curve.index = points.size()-1;
              closest_point_on_curve.parametric_fraction = 0;
              closest_point_on_curve.interpolation_fraction = 0;
              closest_point_on_curve.point = points[points.size()-1];
            }
        }
      return closest_point_on_curve;
    }


    double *
    BezierCurve::solve_cubic_equation_real(const double a_original,const double b_original,const double c_original,const double d_original, double real_roots[4])
    {
      // todo: look at https://stackoverflow.com/questions/2003465/fastest-numerical-solution-of-a-real-cubic-polynomial
      //real_roots = {NaN::DSNAN,NaN::DSNAN,NaN::DSNAN,0.};
      //real_roots.reserve(3);
      //constexpr double one_third = 1./3.;
      size_t index = real_roots[3];
      constexpr double tolerance = 1e-10;
      if (std::abs(a_original) <= tolerance)
        {
          if (std::abs(b_original) <= tolerance)
            {
              // linear equation
              const double &a = c_original;
              const double &b = d_original;
              real_roots[index] = -b/a;
              index++;
              real_roots[3] = index;
              return real_roots;
            }
          // quadratic equation
          const double &a = b_original;
          const double &b = c_original;
          const double &c = d_original;
          const double discriminant = b*b -4.*a*c;
          if (std::abs(discriminant) <= tolerance)
            {
              real_roots[index] = (-b+sqrt(discriminant))/(2*a);
              index++;
              real_roots[3] = index;
              return real_roots;
            }
          else if (discriminant > 0)
            {
              real_roots[index] = (-b + sqrt(discriminant))/(2.*a);
              real_roots[index+1] = (-b - sqrt(discriminant))/(2.*a);
              index += 2;
              real_roots[3] = index;
              return real_roots;
            }
          return real_roots;
        }
      else
        {
          // based on https://quarticequations.com/Cubic.pdf
          // divide by a
          //constexpr double one_third = 1./3.;
          const double b = b_original/a_original; // a_2
          const double c = c_original/a_original; // a_1
          const double d = d_original/a_original; //a_0

          const double q = (c/3.)-(b*b/9.);
          const double r = (c*b-3*d)/6. - b*b*b/27.;
          const double discriminant = r*r+q*q*q;

          if (discriminant > 0)
            {
              // only one real solution
              const double A = std::pow(std::abs(r) + sqrt(discriminant),1./3.);
              const double t = r >= 0 ? A-q/A : q/A-A;
              real_roots[0] = t-b/3.;
              index++;
            }
          else
            {
              // three real solutions
              // std::pow(-q,3./2.) == sqrt(-q*q*q)
              const double phi_1 = std::abs(q) <= tolerance ? 0 : (acos(r/sqrt(-q*q*q))/3.);
              const double phi_2 = phi_1 - 2.*Consts::PI/3.;
              const double phi_3 = phi_1 + 2.*Consts::PI/3.;
              const double sqrt_q_3 = 2*sqrt(-q);
              const double b_t_one_third = b/3.;
              const double value_1 = sqrt_q_3 * FT::cos(phi_1)-b_t_one_third;
              const double value_2 = sqrt_q_3 * FT::cos(phi_2)-b_t_one_third;
              const double value_3 = sqrt_q_3 * FT::cos(phi_3)-b_t_one_third;

              real_roots[0] = value_1;
              index++;
              if (std::abs(value_1 - value_2) > tolerance)
                {
                  real_roots[index] = value_2;
                  index++;
                }
              // we don't have to check value 1 and 3 because z3 <= z2 <= z1
              // so if z2 and z3 are not equal, z1 and z3 are not either.
              if (std::abs(value_2 - value_3) > tolerance)
                {
                  real_roots[index] = value_3;
                  index++;
                }
            }
        }
      real_roots[3] = index;
      return real_roots;
    }

  }
}