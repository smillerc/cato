# -*- coding: utf-8 -*-
from configparser import ConfigParser

parser = ConfigParser()

input_file = "input.ini"
parser.read(input_file)
parser.set("scheme", "limiter", "superbee")

# Writing our configuration file to 'example.ini'
with open(input_file, "w") as configfile:
    parser.write(configfile)
