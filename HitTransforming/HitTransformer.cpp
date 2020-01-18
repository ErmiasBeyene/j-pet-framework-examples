/**
 *  @copyright Copyright 2020 The J-PET Framework Authors. All rights reserved.
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may find a copy of the License in the LICENCE file.
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 *  @file HitTransformer.cpp
 */

#include "../LargeBarrelAnalysis/EventCategorizerTools.h"
#include <JPetHit/JPetHit.h>
#include "HitTransformer.h"
#include "JPetPhysHit.h"
#include <iostream>
#include <string>

using namespace jpet_options_tools;

using namespace std;

HitTransformer::HitTransformer(const char *name) : JPetUserTask(name) {}

bool HitTransformer::init() {
  INFO("Transforming hits started");

  std::string output_class_name = "JPetPhysHit";
  fOutputEvents = new JPetTimeWindow(output_class_name.c_str());

  return true;
}

bool HitTransformer::exec() {
  if (auto timeWindow = dynamic_cast<const JPetTimeWindow *const>(fEvent)) {
    for (uint i = 0; i < timeWindow->getNumberOfEvents(); i++) {
      const auto& hit = dynamic_cast<const JPetHit&>(timeWindow->operator[](i));

      auto tot = EventCategorizerTools::calculateTOT(hit);
      JPetPhysHit physHit(
        hit.getTime(), tot, hit.getPos(), hit.getBarrelSlot().getTheta()
      );
      fOutputEvents->add<JPetPhysHit>(physHit);
    }

  } else {
    return false;
  }
  return true;
}

bool HitTransformer::terminate() {
  INFO("Transformig Hits finished.");
  return true;
}
