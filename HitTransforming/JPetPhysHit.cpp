/**
 *  @copyright Copyright 2020 The J-PET Framework Authors. All rights reserved.
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may find a copy of the License in the LICENCE file.
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *  @file JPetPhysHit.cpp
 */

#include "JPetPhysHit.h"

ClassImp(JPetPhysHit);

JPetPhysHit::JPetPhysHit() : TObject() { /**/ }


JPetPhysHit::JPetPhysHit(
  float time, float tot, TVector3 position, float slotTheta
): fTime(time), fTOT(tot), fPosition(position),fSlotTheta(slotTheta) {}

float JPetPhysHit::getTime() const
{
  return fTime;
}

float JPetPhysHit::getTOT() const
{
  return fTOT;
}

const TVector3& JPetPhysHit::getPos() const
{
  return fPosition;
}

float JPetPhysHit::getSlotTheta() const
{
  return fSlotTheta;
}

void JPetPhysHit::setTime(float time)
{
  fTime = time;
}

void JPetPhysHit::setTOT(float tot)
{
  fTOT = tot;
}

void JPetPhysHit::setPos(float x, float y, float z)
{
  fPosition.SetXYZ(x, y, z);
}
