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

#include "world_builder/features/continental_plate_models/composition/litho1.0.h"

#include "world_builder/assert.h"
#include "world_builder/consts.h"
#include "world_builder/nan.h"
#include "world_builder/objects/natural_coordinate.h"
#include "world_builder/types/array.h"
#include "world_builder/types/double.h"
#include "world_builder/types/object.h"
#include "world_builder/types/one_of.h"
#include "world_builder/types/unsigned_int.h"
#include "world_builder/types/value_at_points.h"
#include <limits>

namespace WorldBuilder {

using namespace Utilities;

namespace Features {
using namespace FeatureUtilities;
namespace ContinentalPlateModels {
namespace Composition {

unsigned int Litho1_0::get_litho1_0_composition(
    const Objects::NaturalCoordinate &position_in_natural_coordinates) const {
  int level = 7;
  int mode = 1; // point mode
  float depth0 = position_in_natural_coordinates.get_depth_coordinate();
  float lon0 =
      position_in_natural_coordinates
          .get_surface_coordinates()[0]; // TODO: needs to be in radians, check
  float lat0 =
      position_in_natural_coordinates
          .get_surface_coordinates()[1]; // TODO: needs to be in radians, check

  unsigned int n1 = 4 * n1 - 6;
  float minlat1 = 0, minlon1 = 0, mindist1 = 0;
  int minnode1 = 0;

  float minlat2 = 0, minlon2 = 0, mindist2 = 0;
  int minnode2 = 0;

  float minlat3 = 0, minlon3 = 0, mindist3 = 0;
  int minnode3 = 0;

  mindist1 = 1.0e5;
  mindist2 = 2.0e5;
  mindist3 = 3.0e5;

  char tessfile[200];
  sprintf(tessfile, "%s/Icosahedron_Level7_LatLon_mod.txt",
          model_path_data.c_str());
  FILE *fp = nullptr;
  if ((fp = fopen(tessfile, "r")) == nullptr) {
    WBAssertThrow(false, "ERROR: Could not open file " << tessfile);
  }

  int node = 0;
  float latitude = 0;
  float glatitude = 0;
  float longitude = 0;
  float lat1 = 0;
  float lon1 = 0;
  float dlat = 0;
  float dlon = 0;
  float radius = 6371.; // TODO: check if I want to connect this to the radius
                        // of the planet in the input file.
  // TODO: preprocess the file
  while (fscanf(fp, "%f %f %f", &latitude, &glatitude, &longitude) != EOF) {

    /* read in a point and convert it to radians to compare */

    ++node;
    /* fprintf(stdout,"%d %f %f %f\n", node, latitude, glatitude, longitude); */

    lat1 = latitude * Consts::PI / 180.;
    lon1 = longitude * Consts::PI / 180.;

    /* calculate the delta lat and lon from the target point */
    dlat = lat0 - lat1;
    dlon = lon0 - lon1;
    float a = sin(dlat / 2) * sin(dlat / 2) +
              cos(lat1) * cos(lat0) * sin(dlon / 2) * sin(dlon / 2);
    float dist = radius * 2 * atan2(sqrt(a), sqrt(1 - a));

    if (dist < mindist1 && node <= n1) {
      mindist3 = mindist2;
      minlat3 = minlat2;
      minlon3 = minlon2;
      minnode3 = minnode2;

      mindist2 = mindist1;
      minlat2 = minlat1;
      minlon2 = minlon1;
      minnode2 = minnode1;

      mindist1 = dist;
      minlat1 = lat1;
      minlon1 = lon1;
      minnode1 = node;

    } else if (dist < mindist2 && node <= n1) {
      mindist3 = mindist2;
      minlat3 = minlat2;
      minlon3 = minlon2;
      minnode3 = minnode2;

      mindist2 = dist;
      minlat2 = lat1;
      minlon2 = lon1;
      minnode2 = node;
    } else if (dist < mindist3 && node <= n1) {
      mindist3 = dist;
      minlat3 = lat1;
      minlon3 = lon1;
      minnode3 = node;
    }
  }

  // TODO: why are we converting back? can't we stay in radians?
  // minlat1 /= Consts::PI / 180.0;
  // minlon1 /= Consts::PI / 180.0;

  // minlat2 /= Consts::PI / 180.0;
  // minlon2 /= Consts::PI / 180.0;

  // minlat3 /= Consts::PI / 180.0;
  // minlon3 /= Consts::PI / 180.0;

  /* get weights of Barycentric coordinate system */

  double lambda1, lambda2, lambda3;

  lambda1 = ((minlon2 - minlon3) * (lat0 - minlat3) +
             (minlat3 - minlat2) * (lon0 - minlon3)) /
            ((minlon2 - minlon3) * (minlat1 - minlat3) +
             (minlat3 - minlat2) * (minlon1 - minlon3));
  lambda2 = ((minlon3 - minlon1) * (lat0 - minlat3) +
             (minlat1 - minlat3) * (lon0 - minlon3)) /
            ((minlon2 - minlon3) * (minlat1 - minlat3) +
             (minlat3 - minlat2) * (minlon1 - minlon3));
  lambda3 = 1 - lambda1 - lambda2;

  // now that we know the correct files to look in, open them:
  FILE *fp1, *fp2, *fp3;
  int nlayers1, nlayers2, nlayers3;
  char modelfile1[200];
  char modelfile2[200];
  char modelfile3[200];
  earthModel model1, model2, model3;

  sprintf(modelfile1, "%s/node%d.model", model_path_data.c_str(), minnode1);
  if ((fp1 = fopen(modelfile1, "r")) == nullptr) {
    fprintf(stdout, "ERROR 1: Could not open file %s\n", modelfile1);
    exit(1);
  }
  fscanf(fp1, "%*s %*s %d", &nlayers1);
  model1.numlayers = nlayers1;

  sprintf(modelfile2, "%s/node%d.model", model_path_data.c_str(), minnode2);
  if ((fp2 = fopen(modelfile2, "r")) == nullptr) {
    fprintf(stdout, "ERROR 2: Could not open file %s\n", modelfile2);
    exit(1);
  }
  fscanf(fp2, "%*s %*s %d", &nlayers2);
  model2.numlayers = nlayers2;

  sprintf(modelfile3, "%s/node%d.model", model_path_data.c_str(), minnode3);
  if ((fp3 = fopen(modelfile3, "r")) == nullptr) {
    fprintf(stdout, "ERROR 3: Could not open file %s\n", modelfile3);
    exit(1);
  }
  fscanf(fp3, "%*s %*s %d", &nlayers3);
  model3.numlayers = nlayers3;

  int i = 0;
  while (fscanf(fp1, "%f %f %f %f %f %f %f %f %f %s", &model1.layers[i].depth,
                &model1.layers[i].density, &model1.layers[i].pvel,
                &model1.layers[i].svel, &model1.layers[i].qkappa,
                &model1.layers[i].qshear, &model1.layers[i].pvel2,
                &model1.layers[i].svel2, &model1.layers[i].eta,
                model1.layers[i].layertype) != EOF) {
    ++i;
  }
  i = 0;
  while (fscanf(fp2, "%f %f %f %f %f %f %f %f %f %s", &model2.layers[i].depth,
                &model2.layers[i].density, &model2.layers[i].pvel,
                &model2.layers[i].svel, &model2.layers[i].qkappa,
                &model2.layers[i].qshear, &model2.layers[i].pvel2,
                &model2.layers[i].svel2, &model2.layers[i].eta,
                model2.layers[i].layertype) != EOF) {
    ++i;
  }
  i = 0;
  while (fscanf(fp3, "%f %f %f %f %f %f %f %f %f %s", &model3.layers[i].depth,
                &model3.layers[i].density, &model3.layers[i].pvel,
                &model3.layers[i].svel, &model3.layers[i].qkappa,
                &model3.layers[i].qshear, &model3.layers[i].pvel2,
                &model3.layers[i].svel2, &model3.layers[i].eta,
                model3.layers[i].layertype) != EOF) {
    ++i;
  }

  fclose(fp1);
  fclose(fp2);
  fclose(fp3);

  /* use weights in interpolating all values */

  /* lets try and find a specific layer */

  char **layertype;
  layertype = (char **)calloc(250, sizeof(char *));
  for (i = 0; i < 250; ++i) {
    layertype[i] = (char *)calloc(20, sizeof(char));
  }

  int k = 0;

  char string[20];
  for (i = 0; i <= 24; ++i) {
    sprintf(string, "IC%d", i);
    strcpy(layertype[++k], string);
  }
  for (i = 0; i <= 45; ++i) {
    sprintf(string, "OC%d", i);
    strcpy(layertype[++k], string);
  }
  for (i = 0; i <= 71; ++i) {
    sprintf(string, "M%d", i);
    strcpy(layertype[++k], string);
  }

  strcpy(layertype[++k], "A-BOTTOM");
  strcpy(layertype[++k], "A-TOP");

  strcpy(layertype[++k], "ASTHENO-BOTTOM");
  strcpy(layertype[++k], "ASTHENO-TOP");

  strcpy(layertype[++k], "LID-BOTTOM");
  strcpy(layertype[++k], "LID-TOP");

  strcpy(layertype[++k], "CRUST3-BOTTOM");
  strcpy(layertype[++k], "CRUST3-TOP");
  strcpy(layertype[++k], "CRUST2-BOTTOM");
  strcpy(layertype[++k], "CRUST2-TOP");
  strcpy(layertype[++k], "CRUST1-BOTTOM");
  strcpy(layertype[++k], "CRUST1-TOP");

  strcpy(layertype[++k], "SEDS3-BOTTOM");
  strcpy(layertype[++k], "SEDS3-TOP");
  strcpy(layertype[++k], "SEDS2-BOTTOM");
  strcpy(layertype[++k], "SEDS2-TOP");
  strcpy(layertype[++k], "SEDS1-BOTTOM");
  strcpy(layertype[++k], "SEDS1-TOP");

  strcpy(layertype[++k], "ICE-BOTTOM");
  strcpy(layertype[++k], "ICE-TOP");

  strcpy(layertype[++k], "WATER-BOTTOM");
  strcpy(layertype[++k], "WATER-TOP");

  float sum, depth, den, pvel, svel, qkappa, qshear, pvel2, svel2, eta;
  float tmp_depth = 0, tmp_den = 0, tmp_pvel = 0, tmp_svel = 0, tmp_qkappa = 0,
        tmp_qshear = 0, tmp_pvel2 = 0, tmp_svel2 = 0, tmp_eta = 0;
  float tmp1_depth = 0, tmp1_den = 0, tmp1_pvel = 0, tmp1_svel = 0,
        tmp1_qkappa = 0, tmp1_qshear = 0, tmp1_pvel2 = 0, tmp1_svel2 = 0,
        tmp1_eta = 0;
  float tmp2_depth = 0, tmp2_den = 0, tmp2_pvel = 0, tmp2_svel = 0,
        tmp2_qkappa = 0, tmp2_qshear = 0, tmp2_pvel2 = 0, tmp2_svel2 = 0,
        tmp2_eta = 0;
  float tmp3_depth = 0, tmp3_den = 0, tmp3_pvel = 0, tmp3_svel = 0,
        tmp3_qkappa = 0, tmp3_qshear = 0, tmp3_pvel2 = 0, tmp3_svel2 = 0,
        tmp3_eta = 0;
  int tmp1_layer = 0, tmp2_layer = 0, tmp3_layer = 0;
  int tmp1_flag = 0, tmp2_flag = 0, tmp3_flag = 0;

  int j;

  for (j = 0; j <= k; ++j) {
    tmp1_flag = 0;
    tmp2_flag = 0;
    tmp3_flag = 0;
    /* if layer does not exist, use depth from previous layer */
    /* use tmp_flag to make sure you don't use the parameter values */

    for (i = 0; i < model1.numlayers; ++i) {
      if (strcmp(model1.layers[i].layertype, layertype[j]) == 0) {
        tmp1_depth = model1.layers[i].depth;
        tmp1_den = model1.layers[i].density;
        tmp1_pvel = model1.layers[i].pvel;
        tmp1_svel = model1.layers[i].svel;
        tmp1_qkappa = model1.layers[i].qkappa;
        tmp1_qshear = model1.layers[i].qshear;
        tmp1_pvel2 = model1.layers[i].pvel2;
        tmp1_svel2 = model1.layers[i].svel2;
        tmp1_eta = model1.layers[i].eta;

        tmp1_flag = 1;
        tmp1_layer = i;
      }
    }
    for (i = 0; i < model2.numlayers; ++i) {
      if (strcmp(model2.layers[i].layertype, layertype[j]) == 0) {
        tmp2_depth = model2.layers[i].depth;
        tmp2_den = model2.layers[i].density;
        tmp2_pvel = model2.layers[i].pvel;
        tmp2_svel = model2.layers[i].svel;
        tmp2_qkappa = model2.layers[i].qkappa;
        tmp2_qshear = model2.layers[i].qshear;
        tmp2_pvel2 = model2.layers[i].pvel2;
        tmp2_svel2 = model2.layers[i].svel2;
        tmp2_eta = model2.layers[i].eta;

        tmp2_flag = 1;
        tmp2_layer = i;
      }
    }
    for (i = 0; i < model3.numlayers; ++i) {
      if (strcmp(model3.layers[i].layertype, layertype[j]) == 0) {
        tmp3_depth = model3.layers[i].depth;
        tmp3_den = model3.layers[i].density;
        tmp3_pvel = model3.layers[i].pvel;
        tmp3_svel = model3.layers[i].svel;
        tmp3_qkappa = model3.layers[i].qkappa;
        tmp3_qshear = model3.layers[i].qshear;
        tmp3_pvel2 = model3.layers[i].pvel2;
        tmp3_svel2 = model3.layers[i].svel2;
        tmp3_eta = model3.layers[i].eta;

        tmp3_flag = 1;
        tmp3_layer = i;
      }
    }

    sum = (lambda1 * tmp1_flag + lambda2 * tmp2_flag + lambda3 * tmp3_flag);

    depth =
        (lambda1 * tmp1_depth + lambda2 * tmp2_depth + lambda3 * tmp3_depth);
    den = (lambda1 * tmp1_flag * tmp1_den + lambda2 * tmp2_flag * tmp2_den +
           lambda3 * tmp3_flag * tmp3_den) /
          sum;
    pvel = (lambda1 * tmp1_flag * tmp1_pvel + lambda2 * tmp2_flag * tmp2_pvel +
            lambda3 * tmp3_flag * tmp3_pvel) /
           sum;
    svel = (lambda1 * tmp1_flag * tmp1_svel + lambda2 * tmp2_flag * tmp2_svel +
            lambda3 * tmp3_flag * tmp3_svel) /
           sum;
    qkappa =
        (lambda1 * tmp1_flag * tmp1_qkappa + lambda2 * tmp2_flag * tmp2_qkappa +
         lambda3 * tmp3_flag * tmp3_qkappa) /
        sum;
    qshear =
        (lambda1 * tmp1_flag * tmp1_qshear + lambda2 * tmp2_flag * tmp2_qshear +
         lambda3 * tmp3_flag * tmp3_qshear) /
        sum;
    pvel2 =
        (lambda1 * tmp1_flag * tmp1_pvel2 + lambda2 * tmp2_flag * tmp2_pvel2 +
         lambda3 * tmp3_flag * tmp3_pvel2) /
        sum;
    svel2 =
        (lambda1 * tmp1_flag * tmp1_svel2 + lambda2 * tmp2_flag * tmp2_svel2 +
         lambda3 * tmp3_flag * tmp3_svel2) /
        sum;
    eta = (lambda1 * tmp1_flag * tmp1_eta + lambda2 * tmp2_flag * tmp2_eta +
           lambda3 * tmp3_flag * tmp3_eta) /
          sum;

    if ((strcmp(layertype[j], "IC0") == 0) &&
        ((tmp1_flag == 0) || (tmp2_flag == 0) || (tmp3_flag == 0))) {
      /* throw an error if there is no IC0 layer, it means that one of the nodes
       * is missing */
      if (tmp1_flag == 0)
        fprintf(stderr, "ERROR: Missing node = %d\n", minnode1);
      if (tmp2_flag == 0)
        fprintf(stderr, "ERROR: Missing node = %d\n", minnode2);
      if (tmp3_flag == 0)
        fprintf(stderr, "ERROR: Missing node = %d\n", minnode3);
      return -1;
    }
    if (sum > 0.0 && (depth / 1000. <= depth0) &&
        (tmp_depth / 1000. > depth0)) {
      /* interpolate to get the results */
      tmp_den = tmp_den + (den - tmp_den) * (depth0 * 1000. - tmp_depth) /
                              (depth - tmp_depth);
      tmp_pvel = tmp_pvel + (pvel - tmp_pvel) * (depth0 * 1000. - tmp_depth) /
                                (depth - tmp_depth);
      tmp_svel = tmp_svel + (svel - tmp_svel) * (depth0 * 1000. - tmp_depth) /
                                (depth - tmp_depth);
      tmp_qkappa = tmp_qkappa + (qkappa - tmp_qkappa) *
                                    (depth0 * 1000. - tmp_depth) /
                                    (depth - tmp_depth);
      tmp_qshear = tmp_qshear + (qshear - tmp_qshear) *
                                    (depth0 * 1000. - tmp_depth) /
                                    (depth - tmp_depth);
      tmp_pvel2 = tmp_pvel2 + (pvel2 - tmp_pvel2) *
                                  (depth0 * 1000. - tmp_depth) /
                                  (depth - tmp_depth);
      tmp_svel2 = tmp_svel2 + (svel2 - tmp_svel2) *
                                  (depth0 * 1000. - tmp_depth) /
                                  (depth - tmp_depth);
      tmp_eta = tmp_eta + (eta - tmp_eta) * (depth0 * 1000. - tmp_depth) /
                              (depth - tmp_depth);

      fprintf(stdout,
              "%7.0f. %8.2f %8.2f %8.2f %7.2f %7.2f %8.2f %8.2f %7.5f %s %s\n",
              depth0 * 1000., tmp_den, tmp_pvel, tmp_svel, tmp_qkappa,
              tmp_qshear, tmp_pvel2, tmp_svel2, tmp_eta, layertype[j - 1],
              layertype[j]);
    }
    if (layertype[j - 1] == "ASTHENO-BOTTOM") {
      return compositions[LayerType::ASTHENOSPHERE];
    } else if (layertype[j - 1] == "LID-BOTTOM") {
      return compositions[LayerType::LITHOSPHERE];

    } else if (layertype[j - 1] == "CRUST3-BOTTOM") {
      return compositions[LayerType::CRUST3];
    } else if (layertype[j - 1] == "CRUST2-BOTTOM") {
      return compositions[LayerType::CRUST2];
    } else if (layertype[j - 1] == "CRUST1-BOTTOM") {
      return compositions[LayerType::CRUST1];
    } else if (layertype[j - 1] == "SEDS3-BOTTOM") {
      return compositions[LayerType::SEDIMENT3];
    } else if (layertype[j - 1] == "SEDS2-BOTTOM") {
      return compositions[LayerType::SEDIMENT2];
    } else if (layertype[j - 1] == "SEDS1-BOTTOM") {
      return compositions[LayerType::SEDIMENT1];
    } else if (layertype[j - 1] == "ICE-BOTTOM") {
      return compositions[LayerType::ICE];
    } else if (layertype[j - 1] == "WATER-BOTTOM") {
      return compositions[LayerType::WATER];
    }
    // if not found, just continue
    /* for point mode, save previous layer before moving on */
    if (sum > 0.0) {
      tmp_depth = depth;
      tmp_den = den;
      tmp_pvel = pvel;
      tmp_svel = svel;
      tmp_qkappa = qkappa;
      tmp_qshear = qshear;
      tmp_pvel2 = pvel2;
      tmp_svel2 = svel2;
      tmp_eta = eta;
    }
  }
  return std::numeric_limits<unsigned int>::quiet_NaN();
}

Litho1_0::Litho1_0(WorldBuilder::World *world_)
    : min_depth(NaN::DSNAN), max_depth(NaN::DSNAN) {
  this->world = world_;
  this->name = "Litho1.0";
}

Litho1_0::~Litho1_0() = default;

void Litho1_0::declare_entries(Parameters &prm,
                               const std::string & /*unused*/) {
  // Document plugin and require entries if needed.
  // Add compositions to the required parameters.
  // prm.declare_entry(
  //    "", Types::Object({"compositions"}),
  //    "Uniform compositional model. Sets constant compositional field.");

  // Declare entries of this plugin
  prm.declare_entry("min depth",
                    Types::OneOf(Types::Double(0),
                                 Types::Array(Types::ValueAtPoints(0., 2.))),
                    "The depth in meters from which the composition of this "
                    "feature is present.");

  prm.declare_entry(
      "max depth",
      Types::OneOf(Types::Double(std::numeric_limits<double>::max()),
                   Types::Array(Types::ValueAtPoints(
                       std::numeric_limits<double>::max(), 2.))),
      "The depth in meters to which the composition of this feature is "
      "present.");

  prm.declare_entry("composition asthenosphere", Types::UnsignedInt(0),
                    "The composition label of the astenophere.");
  prm.declare_entry("composition lithosphere", Types::UnsignedInt(0),
                    "The composition label of the lithosphere (no crust).");
  prm.declare_entry("composition crust 3", Types::UnsignedInt(0),
                    "The composition label of the lowest crust.");
  prm.declare_entry("composition crust 2", Types::UnsignedInt(0),
                    "The composition label of the middle crust.");
  prm.declare_entry("composition crust 1", Types::UnsignedInt(0),
                    "The composition label of the upper crust.");
  prm.declare_entry("composition sediment 3", Types::UnsignedInt(0),
                    "The composition label of the lower sediment.");
  prm.declare_entry("composition sediment 2", Types::UnsignedInt(0),
                    "The composition label of the middle sediment.");
  prm.declare_entry("composition sediment 1", Types::UnsignedInt(0),
                    "The composition label of the upper sediment.");
  prm.declare_entry("composition ice", Types::UnsignedInt(0),
                    "The composition label of the ice.");
  prm.declare_entry("composition water", Types::UnsignedInt(0),
                    "The composition label of the water.");
  prm.declare_entry(
      "model data path", Types::String(""),
      "The path to the Litho1.0 model data, which can be downloaded from...");

  prm.declare_entry("fractions", Types::Array(Types::Double(1.0), 1),
                    "TA list of compositional fractions corresponding to the "
                    "compositions list.");

  prm.declare_entry(
      "operation",
      Types::String("replace",
                    std::vector<std::string>{"replace", "replace defined only",
                                             "add", "subtract"}),
      "Whether the value should replace any value previously defined at this "
      "location (replace) or "
      "add the value to the previously define value. Replacing implies that "
      "all compositions not "
      "explicitly defined are set to zero. To only replace the defined "
      "compositions use the replace only defined option.");
}

void Litho1_0::parse_entries(Parameters &prm,
                             const std::vector<Point<2>> &coordinates) {
  min_depth_surface = Objects::Surface(prm.get("min depth", coordinates));
  min_depth = min_depth_surface.minimum;
  max_depth_surface = Objects::Surface(prm.get("max depth", coordinates));
  max_depth = max_depth_surface.maximum;

  composition_asthenopshere =
      prm.get<unsigned int>("composition asthenosphere");
  composition_lithopshere = prm.get<unsigned int>("composition lithosphere");
  composition_crust_3 = prm.get<unsigned int>("composition crust 3");
  composition_crust_2 = prm.get<unsigned int>("composition crust 2");
  composition_crust_1 = prm.get<unsigned int>("composition crust 1");
  composition_sediment_1 = prm.get<unsigned int>("composition sediment 1");
  composition_sediment_2 = prm.get<unsigned int>("composition sediment 2");
  composition_sediment_3 = prm.get<unsigned int>("composition sediment 3");
  composition_ice = prm.get<unsigned int>("composition ice");
  composition_water = prm.get<unsigned int>("composition water");
  model_path_data = prm.get<std::string>("model data path");
  fractions = prm.get_vector<double>("fractions");
  operation = string_operations_to_enum(prm.get<std::string>("operation"));
}

double Litho1_0::get_composition(
    const Point<3> & /*position_in_cartesian_coordinates*/,
    const Objects::NaturalCoordinate &position_in_natural_coordinates,
    const double depth, const unsigned int composition_number,
    double composition, const double /*feature_min_depth*/,
    const double /*feature_max_depth*/) const {
  if (depth <= max_depth && depth >= min_depth) {
    const double min_depth_local =
        min_depth_surface.constant_value
            ? min_depth
            : min_depth_surface
                  .local_value(
                      position_in_natural_coordinates.get_surface_point())
                  .interpolated_value;
    const double max_depth_local =
        max_depth_surface.constant_value
            ? max_depth
            : max_depth_surface
                  .local_value(
                      position_in_natural_coordinates.get_surface_point())
                  .interpolated_value;
    if (depth <= max_depth_local && depth >= min_depth_local) {
      unsigned int composition_of_layer =
          this->get_litho1_0_composition(position_in_natural_coordinates);
      // for (unsigned int i = 0; i < compositions.size(); ++i) {
      if (composition_of_layer == composition_number) {
        return apply_operation(operation, composition, 1.0);
      }
      //}

      if (operation == Operations::REPLACE)
        return 0.0;
    }
  }
  return composition;
}
WB_REGISTER_FEATURE_CONTINENTAL_PLATE_COMPOSITION_MODEL(Litho1_0, litho1.0)
} // namespace Composition
} // namespace ContinentalPlateModels
} // namespace Features
} // namespace WorldBuilder
