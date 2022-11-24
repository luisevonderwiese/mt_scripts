#!/bin/bash
#NAMES=(zhang_et_al2019 lee_and_hasegawa2013 robinson_and_holton2012 michael_et_al2015 walker_and_ribeiro2011 sicoli_and_holton2014 sagart_et_al2019 lee_and_hasegawa2011 lee2015 kolipakam_et_al2018 kitchen_et_al2009 hruschka_et_al2015 honkola_et_al2013 grollemund_et_al2015 greenhill2015 gray_et_al2009 dunn_et_al2011 defilippo_et_al2012 chang_et_al2015 bowern_and_atkinson2012 bouckaert_et_al2018 bouckaert_et_al2012 atkinson2006 birchall_et_al2016 chacon_and_list2015)
NAMES=(kolipakam_et_al2018 kitchen_et_al2009)
mkdir data

for NAME in ${NAMES[@]}; do
    git clone https://github.com/phlorest/$NAME.git
    cp $NAME/cldf/data.nex data/$NAME.nex
    rm -rf $NAME
done
