#!/bin/bash

sbatch comet.sh caicos 0.01
sbatch comet.sh caicos 0.03
sbatch comet.sh caicos 0.1
sbatch comet.sh caicos 0.3
sbatch comet.sh caicos 1.0
sbatch comet.sh caicos 3.0
sbatch comet.sh caicos 10.0
sbatch comet.sh disc 0.01
sbatch comet.sh disc 0.03
sbatch comet.sh disc 0.1
sbatch comet.sh disc 0.3
sbatch comet.sh disc 1.0
sbatch comet.sh disc 3.0
sbatch comet.sh disc 10.0