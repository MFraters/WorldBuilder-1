/*
  Copyright (C) 2018-2024 by the authors of the World Builder code.

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

#ifndef WORLD_BUILDER_FEATURES_CONTINENTAL_PLATE_MODELS_COMPOSITION_LITHO1_0_H
#define WORLD_BUILDER_FEATURES_CONTINENTAL_PLATE_MODELS_COMPOSITION_LITHO1_0_H

#include "world_builder/features/continental_plate_models/composition/interface.h"
#include "world_builder/features/feature_utilities.h"
#include "world_builder/objects/natural_coordinate.h"
#include "world_builder/objects/surface.h"

namespace WorldBuilder {
namespace Features {
using namespace FeatureUtilities;
namespace ContinentalPlateModels {
namespace Composition {
/**
 * This class represents a continental plate and can implement submodules
 * for temperature and composition. These submodules determine what
 * the returned temperature or composition of the temperature and composition
 * functions of this class will be.
 */
class Litho1_0 final : public Interface {
public:
  /**
   * constructor
   */
  Litho1_0(WorldBuilder::World *world);

  /**
   * Destructor
   */
  ~Litho1_0() override final;

  /**
   * declare and read in the world builder file into the parameters class
   */
  static void declare_entries(Parameters &prm,
                              const std::string &parent_name = "");

  /**
   * declare and read in the world builder file into the parameters class
   */
  void parse_entries(Parameters &prm,
                     const std::vector<Point<2>> &coordinates) override final;

  /**
   * Returns compositional laybel from litho 1.0
   */
  unsigned int get_litho1_0_composition(
      const Objects::NaturalCoordinate &position_in_natural_coordinates) const;

  /**
   * Returns a composition based on the given position, depth in the model,
   * gravity and current composition.
   */
  double get_composition(
      const Point<3> &position,
      const Objects::NaturalCoordinate &position_in_natural_coordinates,
      const double depth, const unsigned int composition_number,
      double composition, const double feature_min_depth,
      const double feature_max_depth) const override final;

private:
  // Litho1.0 composition submodule parameters
  enum LayerType {
    ASTHENOSPHERE,
    LITHOSPHERE,
    CRUST3,
    CRUST2,
    CRUST1,
    SEDIMENT3,
    SEDIMENT2,
    SEDIMENT1,
    ICE,
    WATER,
    MAX_LAYERTYPE
  };
  class earthLayers {
  public:
    float depth;
    float pvel;
    float svel;
    float density;
    float qkappa;
    float qshear;
    float pvel2;
    float svel2;
    float eta;
    char layertype[20];
  };

  class earthModel {
  public:
    int numlayers;
    int num_ic_layers;
    int num_oc_layers;
    earthLayers layers[250];
  };

  double min_depth;
  Objects::Surface min_depth_surface;
  double max_depth;
  Objects::Surface max_depth_surface;
  std::vector<double> fractions;
  Operations operation;
  unsigned int compositions[LayerType::MAX_LAYERTYPE];
  unsigned int composition_asthenopshere;
  unsigned int composition_lithopshere;
  unsigned int composition_crust_3;
  unsigned int composition_crust_2;
  unsigned int composition_crust_1;
  unsigned int composition_sediment_3;
  unsigned int composition_sediment_2;
  unsigned int composition_sediment_1;
  unsigned int composition_ice;
  unsigned int composition_water;
  std::string model_path_data;
};
} // namespace Composition
} // namespace ContinentalPlateModels
} // namespace Features
} // namespace WorldBuilder

#endif
