INPUT=data/sim_itk/muon100
OUTPUT=output/seed_itk_muon100
RESPONSE=tgeo-atlas-itk-hgtd.response
DIGICONFIG=itk-pixel-digitization.json
SPCONFIG=geoSelection-ITk.json


./../build/bin/ActsExampleSeedingTGeo \
    --input-dir=$INPUT \
    --output-dir=$OUTPUT \
    --bf-constant-tesla=0:0:2 \
    --response-file=$RESPONSE \
    --geo-selection-config-file=$SPCONFIG \
    --digi-config-file $DIGICONFIG \
    2>&1 | tee log.txt
