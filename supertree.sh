#!/bin/bash

while getopts "bp" opt
do
    case "${opt}" in
        b) ./fasturec -G data/tree_for_urec/tree_urec_bootstrap.nwk -Y;;
        p) ./fasturec -G data/tree_for_urec/tree_urec_paralogs.nwk -Y;;
    esac
done

./fasturec -G data/tree_for_urec/tree_urec.nwk -Y