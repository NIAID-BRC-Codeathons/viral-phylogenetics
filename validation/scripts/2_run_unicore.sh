#!/bin/bash

TMP=$(mktemp -d)

# run unicore using default ProstT5 model
unicore easy-core fasta unicore/core-prostt5/ unicore/model/ $TMP --max-len 2000 --no-inference
