ddsim --inputFiles gun_3to5.hepmc \
 --numberOfEvents 1000 --outputFile ALFA_o1_v00-onlyCALO_gun_3to5_1000evt.root \
 --random.enableEventSeed --random.seed 1121 \
 --steeringFile ../example/SteeringFile_ALFA_o1_v00_hitpf.py \
 --compactFile ../FCCee/ALFA/compact/ALFA_o1_v00/ALFA_o1_v00-onlyCALO.xml \
 --part.minimalKineticEnergy "1*GeV" --printLevel "ERROR"


k4run run_digi.py
