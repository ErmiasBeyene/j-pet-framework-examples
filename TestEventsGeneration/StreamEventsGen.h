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
 *  @file StreamEventsGen.h
 */

#ifndef STREAMEEVENTSGEN_H
#define STREAMEEVENTSGEN_H

#include <JPetBarrelSlot/JPetBarrelSlot.h>
#include <JPetUserTask/JPetUserTask.h>
#include <JPetEvent/JPetEvent.h>
#include <vector>

class JPetWriter;

#ifdef __CINT__
#define override
#endif

/**
 * @brief StreamEventsGen
 */
class StreamEventsGen: public JPetUserTask
{
public:
  StreamEventsGen(const char * name);
  virtual ~StreamEventsGen(){}
  virtual bool init() override;
  virtual bool exec() override;
  virtual bool terminate() override;

protected:
  void saveEvent(const JPetEvent& event);
};

#endif /* !STREAMEEVENTSGEN_H */
