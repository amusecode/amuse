def read_log(filename):
  infile = open(filename, "r")		# open file
  data = [] 				# declare data container
  while True:				# while data to read
    line = infile.readline()		# read data
    if not line:
      break
    line = line.split()			# distinguish columns
    data.append(line)  			# put data in container
  infile.close()			# close file
  data = zip(*data)			# inverse rows and columns
  return data[2]			# return data

def read_out(filename):
  infile = open(filename, "r")		# open file
  data = [] 				# declare data container
  while True:				# while data to read
    line = infile.readline()		# read data
    if not line:
      break
    line = line.split()			# distinguish columns
    data.append(line)  			# put data in container
  infile.close()			# close file
  data_line = data[4001]
  x0 = float(data_line[1])
  vx0 = float(data_line[4])
  return [x0, vx0]			# return data

def read_xy(filename):
  infile = open(filename, "r")		# open file

  data1 = [] 				# declare data container
  data2 = []
  data3 = []

  numLine = 0
  p1line = 1
  p2line = 2
  p3line = 3
  dp = 4

  while True:				# while data to read
    line = infile.readline()		# read data
    if not line:
      break
    line = line.split()			# distinguish columns

    if numLine == p1line:
      data1.append(line)  		# put data in container
      p1line += dp
    elif numLine == p2line:
      data2.append(line)
      p2line += dp
    elif numLine == p3line:
      data3.append(line)
      p3line += dp

    numLine += 1

  infile.close()			# close file

  data1 = zip(*data1)
  data2 = zip(*data2)
  data3 = zip(*data3)

  x1 = data1[1]
  y1 = data1[2]
  x2 = data2[1]
  y2 = data2[2]
  x3 = data3[1]
  y3 = data3[2]

  x1 = [float(k) for k in x1]
  y1 = [float(k) for k in y1]
  x2 = [float(k) for k in x2]
  y2 = [float(k) for k in y2]
  x3 = [float(k) for k in x3]
  y3 = [float(k) for k in y3]

  return x1, y1, x2, y2, x3, y3 	# return data


