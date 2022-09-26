import os 
from pathlib import Path
import sys


def findfile(name, path):
    _li_file_path = []
    for (dirpath, dir, files) in os.walk(path):
        for filename in files:
            if name in filename:
                _li_file_path.append(dirpath)

    return _li_file_path

if sys.platform.startswith('win32'):
  li_file_path = findfile("Rscript", "C:/") + findfile("Rscript", "D:/") + findfile("Rscript", "/")

  for i in range(len(li_file_path)):
    try:
      file_path = li_file_path[i]
      file_path = file_path.replace('/', '\\')
      path_dir = os.path.dirname(os.path.abspath(__file__))
      os.system('"'+ os.path.join(file_path, "Rscript") + '"' + ' ' + os.path.join(path_dir, "input\haploR.R"))
      
    except:
      print("There's a problem")


else:
  li_file_path = findfile("Rscript", "/")

  for i in range(len(li_file_path)):
    try:
      file_path = li_file_path[i]
      path_dir = os.path.dirname(os.path.abspath(__file__))
      os.system('"'+ os.path.join(file_path, "Rscript") + '"' + ' ' + os.path.join(path_dir, "input/haploR.R"))
      
    except:
      print("There's a problem")