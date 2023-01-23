import math
import numpy as np
from scipy import linalg as spla

# TODO:  Implement AC sweep. Plotting can be done using MatPlotLib

"""
Main circuit class with support for both DC and AC circuits
"""
class Circuit:
	"""
	Class member variables.

	N = number of nodes in the circuit
	Ns = number of sources in the circuit
	w = voltage source frequency in rad/s. Currently all sources must be of the same frequency. For DC, w = 0
	G = conductance matrix, which is built as resistors, inductors, and capacitors are added
	Vin = Voltage source values
	inNode = Nodes to which the sources are connected. Currently one node must be connected to ground.
	vNode = Node voltages which are calculated in the calcNodeVoltages function
	"""
	__slots__ = 'N', 'Ns', 'w', 'G', 'Vin', 'inNode', 'vNode'

	"""
	Create an instance of the Circuit class.
	This creates a new circuit to which sources and
	resistors can be added. As a new source or resistor
	is added the input voltage vector and conductance
	matrix, respectively, are updated accordingly.
	"""
	def __init__ (self, freq):
		self.N = int (0)
		self.Ns = int (0)
		self.w = freq
		self.Vin = np.zeros (0, dtype = np.float64)
		self.inNode = np.zeros (0, dtype = np.int64)
		self.vNode = np.zeros (0, dtype = np.complex128)
		self.G = np.zeros ((0, 0), dtype = np.complex128)

	"""
	Add a DC voltage source to the circuit.
	One of the nodes must be connected to ground (node 0)
	"""
	def addSource (self, n1, n2, value):
		# Check that one of the nodes is connected to ground
		if (n1 != 0) and (n2 != 0):
			raise RuntimeError ("Voltage source must have one node connected to ground.")
#			print ("Voltage source must have one node connected to ground.")
#			return -1

		elif n1 == n2:
			raise RuntimeError ("Positive and negative terminals of a voltage source must be connected to two different terminals")

		# Make n1 the larger node and n2 the ground
		if n1 == 0:
			n1 = n2
			n2 = 0
			value *= -1

		# Increment the value of Ns, which is the number of sources
		self.Ns += 1

		# Add the source to the Vin vector and update the inNode vector
		self.Vin = np.resize (self.Vin, new_shape = self.Ns)
		self.inNode = np.resize (self.inNode, new_shape = self.Ns)

		self.Vin[self.Ns - 1] = value
		self.inNode[self.Ns - 1] = n1

	"""
	Add a resistor to the circuit and update the conductance matrix
	"""
	def addResistor (self, n1, n2, value):
		# If the resistance is zero or negative we cannot proceed
		# so we raise an exception
		if value <= 0:
			raise RuntimeError ("Zero or negative resistance between nodes %d and %d" % (n1, n2))

		elif n1 == n2:
			raise RuntimeError ("Positive and negative terminals of a resistor must be connected to two different nodes")

		# Update the value of N, which is the number of nodes in the circuit
		# If either n1 or n2 is greater than N then update its value to
		# the maximum of n1 and n2. Then resize the conductance matrix and
		# set the new elements to zero.
		if (n1 > self.N) or (n2 > self.N):
			size = np.shape (self.G)[0]
			self.N = max (n1, n2)

			self.G = np.pad (self.G, ((0, self.N - size), (0, self.N - size)), 'constant', constant_values = 0)

		# Check if the new value of N is greater than the number of rows in
		# the conductance matrix. If so then resize the conductance matrix.
#		if self.N > np.shape (self.G)[0]:
#			self.G = np.resize (self.G, new_shape =(self.N, self.N))

		# Add the new resistor to the conductance matrix
		# If one of the nodes is 0, meaning ground, then we
		# only need to add it to its corresponding diagonal element.
		if n1 == 0:
			self.G[n2-1][n2-1] += (1 / value)

		elif n2 == 0:
			self.G[n1-1][n1-1] += (1 / value)

		# Otherwise we need to add it to two diagonals and
		# two off-diagonal elements
		else:
			self.G[n1-1][n1-1] += (1 / value)
			self.G[n2-1][n2-1] += (1 / value)

			self.G[n1-1][n2-1] -= (1 / value)
			self.G[n2-1][n1-1] -= (1 / value)

	"""
	Add an inductor between nodes n1 and n2
	"""
	def addInductor (self, n1, n2, value):
		# First check if the circuit has DC or AC sources.
		# If they are DC then we can't add the capacitor.
		if self.w == 0:
			raise RuntimeError ("Cannot add an inductor to a DC circuit.")

		# If the resistance is zero or negative we cannot proceed
		# so we raise an exception
		if value <= 0:
			raise RuntimeError ("Zero or negative inductance between nodes %d and %d" % (n1, n2))

		elif n1 == n2:
			raise RuntimeError ("Positive and negative terminals of a resistor must be connected to two different nodes")

		# Update the value of N, which is the number of nodes in the circuit
		# If either n1 or n2 is greater than N then update its value to
		# the maximum of n1 and n2. Then resize the conductance matrix and
		# set the new elements to zero.
		if (n1 > self.N) or (n2 > self.N):
			size = np.shape (self.G)[0]
			self.N = max (n1, n2)

			self.G = np.pad (self.G, ((0, self.N - size), (0, self.N - size)), 'constant', constant_values = 0)


		Z = complex (0, self.w * value)

		if n1 == 0:
			self.G[n2-1][n2-1] += (1 / Z)

		elif n2 == 0:
			self.G[n1-1][n1-1] += (1 / Z)

		else:
			self.G[n1-1][n1-1] += (1 / Z)
			self.G[n2-1][n2-1] += (1 / Z)

			self.G[n1-1][n2-1] -= (1 / Z)
			self.G[n2-1][n1-1] -= (1 / Z)

	"""
	Add a capacitor between nodes 1 and 2
	"""
	def addCapacitor (self, n1, n2, value):
		# First check if the circuit has DC or AC sources.
		# If they are DC then we can't add the capacitor.
		if self.w == 0:
			raise RuntimeError ("Cannot add a capacitor to a DC circuit.")

		# If the resistance is zero or negative we cannot proceed
		# so we raise an exception
		if value <= 0:
			raise RuntimeError ("Zero or negative resistance between nodes %d and %d" % (n1, n2))

		elif n1 == n2:
			raise RuntimeError ("Positive and negative terminals of a capacitor must be connected to two different nodes")

		# Update the value of N, which is the number of nodes in the circuit
		# If either n1 or n2 is greater than N then update its value to
		# the maximum of n1 and n2. Then resize the conductance matrix and
		# set the new elements to zero.
		if (n1 > self.N) or (n2 > self.N):
			size = np.shape (self.G)[0]
			self.N = max (n1, n2)

			self.G = np.pad (self.G, ((0, self.N - size), (0, self.N - size)), 'constant', constant_values = 0)


		Z = complex (0, -1 / (self.w * value))

		if n1 == 0:
			self.G[n2-1][n2-1] += (1 / Z)

		elif n2 == 0:
			self.G[n1-1][n1-1] += (1 / Z)

		else:
			self.G[n1-1][n1-1] += (1 / Z)
			self.G[n2-1][n2-1] += (1 / Z)

			self.G[n1-1][n2-1] -= (1 / Z)
			self.G[n2-1][n1-1] -= (1 / Z)

	"""
	Calculate the node voltages by solving the linear system
	of equations created by the conductance matrix and the input voltage vector
	"""
	def calcNodeVoltages (self):
		nodes = self.N
		self.vNode = np.resize (self.vNode, new_shape = (nodes, 1))

		self.vNode[:] = 0

		# Copy the values in Vin to vVector
		for n in range (0, self.Ns):
			i = self.inNode[n] - 1

			self.vNode[i] = self.Vin[n]

			self.G[i][:] = 0
			self.G[i][i] = 1

		# Solve the system of equations and
		# check for exceptions
		try:
			self.vNode = spla.solve (self.G, self.vNode)

		# Check if the matrix is singular
		except LinAlgError:
			raise ("Failed to calculate node voltages. Conductance matrix is singular.")

	"""
	Perform a frequency sweep of the circuit and plot the frequency
	response. Plotting is done using MatPlotLib and output has
	units of dB.

	def freqSweep (wLow, wHigh, n):
		# Check that wLow is a power of 10
		# If it is not then make it the largest
		# power of 10 less than wLow
		l = math.floor (math.log10 (wLow))
		if 10**l != wLow:
			wLow = 10**l

		# Now do the same thing for wHigh
		l = math.floor (math.log10 (wHigh))
		if 10**l != wHigh:
			wHigh = 10**l

		# Determine the number of decades covered by the
		# interval [wLow, wHigh]
		nDecades = math.log10 (wHigh / wLow)

		# Starting frequency for the sweep
		w = wLow

		for i in range (0, nDecades):
			# Determine the sample points using the size of
			# the current range and the specified number
			# of samples per decade
			dp = (10 / n) * wLow

			# Starting frequency for the current decade
			p0 = wLow

			for j in range (0, n):
				w = p0 + dp * j
	"""
	Print the result vector
	"""
	def printNodeVoltages (self):
		i = 1
		for v in self.vNode:
			# If the real part of the voltage is zero then the voltage is purely imaginary
			if v.real == 0:
				# If the imaginary part is also zero then the voltage at the node is zero
				if v.imag == 0:
					print ("V%d = 0<0" % i)

				# If the imaginary part is negative then there is a -90 (270) degree phase shift
				elif v.imag < 0:
					print ("V%d = %f<-90" % (i, -1 * v.imag))

				# Otherwise the imaginary part is positive
				else:
					print ("V%d = %f<90" % (i, v.imag))


			# If the real part of the voltage is negative we can print it normally
			elif v.real < 0:
				# If the imaginary part is zero then we don't need to calculate the phase shift
				if v.imag == 0:
					print ("V%d = %f<0" % (i, v.real))

				# If the imaginary part is negative we can print both normally
				elif v.imag < 0:
					# Calculate the magnitude
					mag = math.sqrt (v.real**2 + v.imag**2)

					# Calculate the phase shift in radians
					angle = math.pi + math.atan (v.imag / v.real)

					# Convert the phase shift from radians to degrees
					angle *= (180 / math.pi)

					print ("V%d = %f<%f" % (i, mag, angle))

				# Otherwise the imaginary part is positive
				else:
					# Calculate the magnitude
					mag = math.sqrt (v.real**2 + v.imag**2)

					# Calculate the phase shift in radians
					angle = math.pi - math.atan (v.imag / v.real)

					# Convert the phase shift from radians to degrees
					angle *= (180 / math.pi)

					print ("V%d = %f<%f" % (i, mag, angle))

			# Otherwise the real part of the voltage is positive
			else:
				# If the imaginary part is zero then there is no phase shift
				if v.imag == 0:
					print ("V%d =  %f<0" % (i, v.real))

				# If the imaginary part is negative
				elif v.imag < 0:
					# Calculate the magnitude
					mag = math.sqrt (v.real**2 + v.imag**2)

					# Calculate the phase shift in radians
					angle = math.atan (v.imag / v.real)

					# Convert the phase shift from radians to degrees
					angle *= (180 / math.pi)

					print ("V%d = %f<%f" % (i, mag, angle))

				# Otherwise both the real and imaginary parts of the voltage are positive
				else:
					# Calculate the magnitude
					mag = math.sqrt (v.real**2 + v.imag**2)

					# Calculate the phase shift in radians
					angle = math.atan (v.imag / v.real)

					# Convert the phase shift from radians to degrees
					angle *= (180 / math.pi)

					print ("V%d = %f<%f" % (i, mag, angle))

			i += 1

