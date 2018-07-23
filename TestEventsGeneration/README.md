# Test events generation example

## Aim
Purose of this example program is to generate sample of J-PET Framework data with JPetEvent objects, that can be later used in testing the reconstruction and streaming prcedures.

## Execution
You need a dummy input file (since Framework requires any input file to start) - to obtain it, contact with the author or make your own. It should be a ROOT file with JPetTimeWindow object in the tree - preferably one object, then the test events generation will be done once (you do not need more). 

Run:
```
./TestEventsGeneration.x -t root -f test..root -l detectorSetupRun2_3_4_5.json -i 2
```
Note the double dot in the file name. It generates the events with pairs of hits - one for every pair of scintillators in the Large Barrel. 


## What next?
The reulst of the execution is a `*.unk.evt.root` file, that you can later use in testing the `EventCategorizer` module.


## Author
[Krzysztof Kacprzak](https://github.com/kkacprzak)

Please report any bugs and suggestions of corrections though [email](k.kacprzak@yahoo.com)