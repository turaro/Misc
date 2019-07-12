import sys, re

#usage: python blast-fix.py [input_file] [output_file]
#use on blast-generated xml files to correct for -parse_deflines
#takes an xml file and moves the query ID to the beginning of the defline for each query ID
def main():
  counter=0
  with open(sys.argv[1]) as f:
    with open(sys.argv[2], 'w') as f2:
      for line in f.readlines():
        if "query-ID>" in line:
          counter+=1
          id = re.search(r'>(.+?)<', line)
          if id:
            id=id.group(1)
          f2.write(line)
        elif "query-def>" in line:
          defline = re.search(r'>(.+?)<', line)
          if defline:
            defline=defline.group(1)
            if id:
              new_defline=id+ " " + defline
              f2.write(line.replace(defline, new_defline))
            else:
              print("error: "), line
        else:
          f2.write(line)

if __name__ == '__main__':
  main()