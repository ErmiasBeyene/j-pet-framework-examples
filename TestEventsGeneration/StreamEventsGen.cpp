/**
 *  @copyright Copyright 2018 The J-PET Framework Authors. All rights reserved.
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
 *  @file StreamEventsGen.cpp
 */

#include <JPetOptionsTools/JPetOptionsTools.h>
#include <JPetWriter/JPetWriter.h>
#include "StreamEventsGen.h"

using namespace jpet_options_tools;
using namespace std;

StreamEventsGen::StreamEventsGen(const char* name): JPetUserTask(name) {}

bool StreamEventsGen::init()
{
  INFO("Test event generation started.");
  fOutputEvents = new JPetTimeWindow("JPetEvent");
  return true;
}

bool StreamEventsGen::exec()
{
  INFO("Execute.");

  auto map1 = getParamBank().getBarrelSlots();
  auto map2 = getParamBank().getBarrelSlots();

  for (auto mapEntry1 : map1) {
    for (auto mapEntry2 : map2) {
      if(mapEntry1.first == mapEntry2.first) { continue; }
        auto hit1 = new JPetHit();
        auto hit2 = new JPetHit();
        hit1->setBarrelSlot(*(mapEntry1.second));
        hit2->setBarrelSlot(*(mapEntry2.second));
        auto event = new JPetEvent();
        event->addHit(*hit1);
        event->addHit(*hit2);
        saveEvent(*event);
        delete hit1;
        delete hit2;
        delete event;
    }
  }
  return true;
}

bool StreamEventsGen::terminate()
{
  INFO("Test event generation ended.");
  return true;
}

void StreamEventsGen::saveEvent(const JPetEvent& event)
{
  fOutputEvents->add<JPetEvent>(event);
}
