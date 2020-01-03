/**
 *  @copyright Copyright 2016 The J-PET Framework Authors. All rights reserved.
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
 *  @file JPetFilterInterface.h
 */

#ifndef _JPetFilterInterface_H_
#define _JPetFilterInterface_H_

/*! \brief Interface that all filters should implement.
*/
class JPetFilterInterface
{
public:
  JPetFilterInterface() {};
  virtual ~JPetFilterInterface(){};
  /*! @brief Returns rescale factor, that should be variable at pos scaled
      @par pos Position of variable rescaled to [0, 1] 
   */
  virtual double operator()(double pos) = 0;
};

#endif /*  !_JPetFilterInterface_H_ */
