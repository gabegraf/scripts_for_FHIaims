#!/usr/bin/env python 3
import numpy as np

def parse_geometry_file(initial_geometry_file):
  with open(f"{initial_geometry_file}", "r") as file:
    lines = file.readlines()

    for line in lines:
      parts = line.split()
      
