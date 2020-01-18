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
 *  @file JPetPhysHit.h
 *  @brief Example of a user-defined data class to store 2-gamma annihilations
 */

#ifndef JPETPHYSHIT_H
#define JPETPHYSHIT_H

#include <TVector3.h>

/**
 * @brief JPetPhysHit class
 *
 */

class JPetPhysHit : public TObject {

public:
  JPetPhysHit();
  JPetPhysHit(
    float time, float tot, TVector3 position, float slotTheta
  );

  float getTime() const;
  float getTOT() const;
  const TVector3& getPos() const;
  float getSlotTheta() const;
  void setTime(float time);
  void setTOT(float tot);
  void setPos(float x, float y, float z);

protected:
  float fTime = 0.0f;
  float fTOT = 0.0f;
  TVector3 fPosition;
  float fSlotTheta = -1.0f;

  ClassDef(JPetPhysHit, 1);
};
#endif /*  !JPETPHYSHIT_H */
