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

#ifndef WORLD_BUILDER_TYPES_CONTOURS_H
#define WORLD_BUILDER_TYPES_CONTOURS_H


#include "world_builder/types/plugin_system.h"


namespace WorldBuilder
{
  class Parameters;

  namespace Types
  {

    /**
     * This class represents a contours value with documentation
     */
    class Contours : public Interface
    {
      public:
        /**
         * A constructor
         */
        Contours(const double default_length,
                 const WorldBuilder::Point<2> &default_thickness,
                 const WorldBuilder::Point<2> &default_top_truncation,
                 const WorldBuilder::Point<2> &default_angle,
                 const Types::Interface &temperature_plugin_system,
                 const Types::Interface &composition_plugin_system,
                 const Types::Interface &grains_plugin_system_);

        /**
         * A constructor for the load_entry function
         */
        Contours(double default_length,
                 WorldBuilder::Point<2> default_thickness,
                 WorldBuilder::Point<2> default_angle,
                 std::string description);

        /**
         * Copy constructor
         */
        Contours(Contours const &other);


        /**
         * Destructor
         */
        ~Contours() override;

        /**
         * Todo
         */
        void write_schema(Parameters &prm,
                          const std::string &name,
                          const std::string &documentation) const override final;


        double value_length;
        double default_length;
        WorldBuilder::Point<2> value_thickness;
        WorldBuilder::Point<2> default_thickness;
        WorldBuilder::Point<2> default_top_truncation;
        WorldBuilder::Point<2> value_angle;
        WorldBuilder::Point<2> default_angle;
        std::unique_ptr<Types::Interface> temperature_plugin_system;
        std::unique_ptr<Types::Interface> composition_plugin_system;
        std::unique_ptr<Types::Interface> grains_plugin_system;

      protected:
        Contours *clone_impl() const override final
        {
          return new Contours(*this);
        };
      private:

    };
  } // namespace Types
  /*
    namespace Objects
    {

      / **
        * This class represents an actual contours
        * /
      template <class A, class B, class C>
      class Contours final: public Types::Interface
      {
        public:

          / **
           * A constructor for the clone and set_entry function
           * /
          Contours(const double default_length,
                  const WorldBuilder::Point<2> &default_thickness,
                  const WorldBuilder::Point<2> &default_top_truncation,
                  const WorldBuilder::Point<2> &default_angle,
                  std::vector<std::shared_ptr<A> > temperature_systems,
                  std::vector<std::shared_ptr<B> > composition_systems,
                  std::vector<std::shared_ptr<C> > grains_systems);

          / **
           * Copy constructor
           * /
          Contours(Contours const &other);

          / **
           * Destructor
           * /
          ~Contours() override final;

          / **
           * Todo
           * /
          void write_schema(Parameters &prm,
                            const std::string &name,
                            const std::string &documentation) const override final;


          double value_length;
          double default_length;
          WorldBuilder::Point<2> value_thickness;
          WorldBuilder::Point<2> value_top_truncation;
          WorldBuilder::Point<2> value_angle;
          std::vector<std::shared_ptr<A> > temperature_systems;
          std::vector<std::shared_ptr<B> > composition_systems;
          std::vector<std::shared_ptr<C> > grains_systems;

        protected:
          Contours *clone_impl() const override final
          {
            return new Contours(*this);
          };
        private:

      };
    }*/ // namespace Objects
} // namespace WorldBuilder

#endif
