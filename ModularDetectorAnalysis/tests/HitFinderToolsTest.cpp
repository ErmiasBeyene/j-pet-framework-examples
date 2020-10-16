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
 *  @file HitFinderToolsTest.cpp
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE HitFinderToolsTest

#include <JPetMatrixSignal/JPetMatrixSignal.h>
#include <JPetWLS/JPetWLS.h>

#include <boost/test/unit_test.hpp>
#include "../HitFinderTools.h"

BOOST_AUTO_TEST_SUITE(SignalFinderTestSuite)

BOOST_AUTO_TEST_CASE(matchSignals_test)
{
  JPetLayer layer1(1, "black module", 11.1);
  JPetLayer layer2(2, "red module", 22.2);
  JPetLayer layer3(3, "wls layer", 33.3);

  JPetSlot slot1(1, 0.0, JPetSlot::Module);
  slot1.setLayer(layer1);
  JPetSlot slot2(2, 90.0, JPetSlot::Module);
  slot2.setLayer(layer2);
  JPetSlot slot3(3, 90.0, JPetSlot::WLS);
  slot3.setLayer(layer3);

  JPetScin scin1(1, 50.0, 10.0, 0.5, 1.1, 2.2, 3.3);
  scin1.setSlot(slot1);
  JPetScin scin2(2, 50.0, 10.0, 0.5, 1.1, 2.2, 3.3);
  scin2.setSlot(slot2);
  JPetWLS wls1(1, 1.1, 2.2, 3.3, 1.1, 2.2, 3.3);
  wls1.setSlot(slot3);
  JPetWLS wls2(2, 1.1, 2.2, 3.3, 1.1, 2.2, 3.3);
  wls2.setSlot(slot3);

  // Matrix objects - scintillator side A and B, two WLS matrices
  std::vector<int> a1ids = {1,2,3,4};
  std::vector<int> b1isd = {5,6,7,8};
  std::vector<int> a2ids = {9,10,11,12};
  std::vector<int> b2ids = {13,14,15,16};
  std::vector<int> wls1ids = {17,18};
  std::vector<int> wls2ids = {19,20,21};
  JPetMatrix mtx1(1, "SideA", a1ids);
  mtx1.setScin(scin1);
  JPetMatrix mtx2(2, "SideB", b1isd);
  mtx2.setScin(scin1);
  JPetMatrix mtx3(3, "SideA", a2ids);
  mtx3.setScin(scin2);
  JPetMatrix mtx4(4, "SideB", b2ids);
  mtx4.setScin(scin2);
  JPetMatrix mtx5(5, "WLS", wls1ids);
  mtx5.setWLS(wls1);
  JPetMatrix mtx6(6, "WLS", wls2ids);
  mtx6.setWLS(wls2);

  // Matrix signals for hit creation from scintillator side A and B
  // hit 1 from layer 2
  JPetMatrixSignal sig0(10.0);
  sig0.setMatrix(mtx3);
  JPetMatrixSignal sig1(12.0);
  sig1.setMatrix(mtx4);

  // hit 2 from layer 1
  JPetMatrixSignal sig2(20.0);
  sig2.setMatrix(mtx2);
  JPetMatrixSignal sig3(21.2);
  sig3.setMatrix(mtx2);
  JPetMatrixSignal sig4(21.4);
  sig4.setMatrix(mtx2);
  JPetMatrixSignal sig5(22.0);
  sig5.setMatrix(mtx1);

  // no hit - not in time coincidece with others
  JPetMatrixSignal sig6a(30.0);
  JPetMatrixSignal sig6b(31.0);
  JPetMatrixSignal sig6c(32.0);
  JPetMatrixSignal sig6d(33.0);
  sig6a.setMatrix(mtx4);
  sig6b.setMatrix(mtx4);
  sig6c.setMatrix(mtx4);
  sig6d.setMatrix(mtx4);

  // hit 3 from layer 2
  JPetMatrixSignal sig7(40.0);
  sig7.setMatrix(mtx3);
  JPetMatrixSignal sig8(42.0);
  sig8.setMatrix(mtx4);
  JPetMatrixSignal sig9(42.5);
  sig9.setMatrix(mtx4);
  JPetMatrixSignal sig10(43.0);
  sig10.setMatrix(mtx4);

  // Matrix signals from WLS
  // for hit 1
  JPetMatrixSignal sig11(11.0);
  sig11.setMatrix(mtx5);
  // no fit
  JPetMatrixSignal sig12(19.0);
  sig12.setMatrix(mtx5);
  JPetMatrixSignal sig13(23.0);
  sig13.setMatrix(mtx6);
  JPetMatrixSignal sig14(24.0);
  sig14.setMatrix(mtx5);
  // for hit 3
  JPetMatrixSignal sig15(41.0);
  sig15.setMatrix(mtx6);
  // no fit
  JPetMatrixSignal sig16(67.0);
  sig16.setMatrix(mtx5);

  // adding signals to std::vectors, undered in time
  std::vector<JPetMatrixSignal> scinSignals;
  scinSignals.push_back(sig2);
  scinSignals.push_back(sig6d);
  scinSignals.push_back(sig5);
  scinSignals.push_back(sig6b);
  scinSignals.push_back(sig4);
  scinSignals.push_back(sig0);
  scinSignals.push_back(sig10);
  scinSignals.push_back(sig8);
  scinSignals.push_back(sig3);
  scinSignals.push_back(sig6a);
  scinSignals.push_back(sig6c);
  scinSignals.push_back(sig1);
  scinSignals.push_back(sig9);
  scinSignals.push_back(sig7);

  std::vector<JPetMatrixSignal> wlsSignals;
  wlsSignals.push_back(sig13);
  wlsSignals.push_back(sig11);
  wlsSignals.push_back(sig16);
  wlsSignals.push_back(sig15);
  wlsSignals.push_back(sig14);
  wlsSignals.push_back(sig12);

  JPetStatistics stats;
  auto result = HitFinderTools::matchSignals(scinSignals, wlsSignals, 0.0, 4.0, stats, false);

  double epsilon = 0.0001;
  BOOST_REQUIRE_EQUAL(result.size(), 3);
  BOOST_REQUIRE_CLOSE(result.at(0).getTime(), 11.0, epsilon);
  BOOST_REQUIRE_CLOSE(result.at(1).getTime(), 21.0, epsilon);
  BOOST_REQUIRE_CLOSE(result.at(2).getTime(), 41.0, epsilon);

  BOOST_REQUIRE(result.at(0).isSignalASet());
  BOOST_REQUIRE(result.at(0).isSignalBSet());
  BOOST_REQUIRE(result.at(0).isSignalWLSSet());

  BOOST_REQUIRE_CLOSE(result.at(0).getSignalA().getTime(), 10.0, epsilon);
  BOOST_REQUIRE_CLOSE(result.at(0).getSignalB().getTime(), 12.0, epsilon);
  BOOST_REQUIRE_CLOSE(result.at(0).getSignalWLS().getTime(), 11.0, epsilon);

  BOOST_REQUIRE(result.at(1).isSignalASet());
  BOOST_REQUIRE(result.at(1).isSignalBSet());
  BOOST_REQUIRE(!result.at(1).isSignalWLSSet());

  BOOST_REQUIRE_CLOSE(result.at(1).getSignalB().getTime(), 20.0, epsilon);
  BOOST_REQUIRE_CLOSE(result.at(1).getSignalA().getTime(), 22.0, epsilon);

  BOOST_REQUIRE(result.at(2).isSignalASet());
  BOOST_REQUIRE(result.at(2).isSignalBSet());
  BOOST_REQUIRE(result.at(2).isSignalWLSSet());

  BOOST_REQUIRE_CLOSE(result.at(2).getSignalA().getTime(), 40.0, epsilon);
  BOOST_REQUIRE_CLOSE(result.at(2).getSignalB().getTime(), 42.0, epsilon);
  BOOST_REQUIRE_CLOSE(result.at(2).getSignalWLS().getTime(), 41.0, epsilon);
}

BOOST_AUTO_TEST_SUITE_END()
