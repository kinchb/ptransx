import json

config_file_name = 'config.json'
with open(config_file_name, 'r') as f:
    config = json.load(f)

import os
import sys

cwd = os.getcwd()
sys.path.insert(0, cwd + '/ptxlib')
sys.path.insert(0, cwd + '/ptxlib/fsolver')
sys.path.insert(0, cwd + '/ptxlib/pandurata')
sys.path.insert(0, cwd + '/ptxlib/response')
sys.path.insert(0, cwd + '/ptxlib/xstarcomm')

failed = False

try:
    import _fsolver
except:
    print("!!! FAILED TO IMPORT FSOLVER !!!")
    failed = True

try:
    import _pandurata
except:
    print("!!! FAILED TO IMPORT PANDURTA !!!")
    failed = True

try:
    import _response
except:
    print("!!! FAILED TO IMPORT RESPONSE !!!")
    failed = True

try:
    import _xstarcomm
except:
    print("!!! FAILED TO IMPORT XSTARCOMM !!!")
    failed = True

if failed:
    print("!!! THE BUILD PROCEDURE FAILED !!!")
else:
    print("Congratulations! You can now run PTRANSX+Pandurata!")
