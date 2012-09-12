"""
A tool for reading files that sue the Fortran "unformatted" format.
The class FortranFile.FortranFile is a subclass of file with additional methods to read such files into numpy arrays.
See its docstring for more info.

The "unformat" is a series of records that look like:
[number of bytes in following data] [data] [number of bytes in preceding data]
where the numbers of bytes are usually recorded as unsigned 32bit integers.  This is not standardized, though, so you can select this if your fortran compiler does something different.  In particular, earlier versions of gfortran

Like all Fortran-later language interoperability code documentation should, this docstring ends with an exhortation to please stop using fortran.
"""

import numpy as np
import struct


class SentinelError(IOError):
	"""
	An IOError subclass for cases where the sentinel bytes in Fortran unformatted files are incorrect.
	"""
	pass

class FortranFile(file):
	"""
	A subclass of the buitin file type, with additional methods to read and write numpy 
	arrays from/to files in the Fortran "unformatted" format.
	
	The format does not seem to be standardized, but always seems to consist of a sentinel indicating the
	size in bytes of the following field, then the data, then the data, and then the same sentinel.
	
	For newer versions of gfortran (4.2 and later) and intel fortran, the sentinel is a 32-bit integer.

	If a read of an array fails the file is returned to where it started before the read, and an error raised.
	This is to help with exploring unknow formats.
	"""
	def __init__(self,name, mode='r', buffering=0,endian="=",sentinel=np.uint32):
		"""Open the file for reading or writing.
		Parameters
		__________
		name : the filename
		mode : the read/write mode ("r","w","rw", etc.)
		buffering : 0 for unbuffered, 1 for line buffers, >1 for buffer size
		
		endian : '=' for native (default)
				 '!' for non-native
				 '>' for big-endian
				 '<' for little-endian
				 This should depend on the architecture of the system
				 where the file was written, or is to be read.  If you are not
				 transferring files between machines then you can always leave this as the default.
		sentinel : This specifies the data type of the sentinel used in the unformatted file.
					  The default is correct at least for gfortran versions >= 4.2 and ifort, and probably lots of others.
			
		Returns
		-------
		None
				
		Raises
		------
		ValueError : If an endianness other than '=','!','>','<' is specified
		IOError : If raised by the superclass constructor.
		TypeError : If raised by the superclass constructor.
		"""
		file.__init__(self, name, mode, buffering)
		self.endian=endian
		self.sentinel=sentinel
		try:
			self.swap = {
			"=":      {True:False, False:False},  #native byte order; never swap
			"!":      {True:True,  False:True},   #non-native byte order; always swap
			">":      {True:True,  False:False},  #big endian; swap if system little-endian
			"<":      {True:False, False:True},   #little endian; swap if system big-endian
			}[endian] [np.little_endian]
		except KeyError:
			raise ValueError("Endianness must be in [=,<,>,!], not %s" % endian)
		
	def _readSentinel(self):
		"""
		Internal method to red the fortran sentinel bytes that surround a record.
		Parameters
		----------
		None
		
		Returns
		_______
		l : integer, of type self.sentinel
		
		Raises
		------
		SentinelError
			if sentinel could not be read (end of file?)
		"""
		l=np.fromfile(self,self.sentinel,1)
		if self.swap:
			l=l.byteswap()
		if not len(l)==1:
			raise SentinelError("Fortran sentinel not found.")
		return int(l[0])  #have to convert to int owing to what seems to be a bug in integer division of uint64s in numpy

	def _checkSentinel(self,x):
		"""
		Internal method to check that the fortran sentinel number is correct.
		Parameters
		----------
		x : int, the sentinel to compare
		
		Returns
		-------
		None
		
		Raises
		------
		IOError
			if unable to read from self.
		SentinelError
			if sentinel is incorrect, indicating file is not correctly formatted
		"""
		x2=self._readSentinel()
		if x2!=x:
			raise SentinelError("Incorrect sentinel in Fortran File")

	def readString(self):
		"""Read in a string from a Fortran file.
		Parameters
		----------
		None
		
		Returns
		-------
		s : string, read from file.
		
		Raises
		------
		IOError : If unable to read array from file or file is not of correct fortran format.
				  If this is raised the file position returns to its position before calling this function.
		
		"""
		pos=self.tell()
		try:
			l=self._readSentinel()
			s = self.read(l)
			self._checkSentinel(l)
		except IOError, e:
			self.seek(pos)
			raise e("String not read correctly from fortran file.  Incorrect format.")
		return s

	def writeString(self,s):
		"""Write a string to a fortran file.
		Parameters
		----------
		s : String to be written to file
		
		Returns
		-------
		None
		
		Raises
		------
		IOError : If unable to write to file.
		
		"""
		sentinel = np.array([len(s)],dtype=self.sentinel)
		if self.swap:
			sentinel=sentinel.byteswap()
		sentinel.tofile(self)
		self.write(s)
		sentinel.tofile(self)


	def readArrays(self,dtypes,shapes=None):
		"""
		Read a number of numpy arrays from a fortran file.
		If any of the reads fail the file pointer remains where it started.
		
		Parameters
		----------
		dtypes : The numpy data types of the arrays to be read.
		shapes : (optional) shapes to reshape the read arrays into
		
		Returns
		-------
		[arrays] : list of arrays of data types [dtypes] read from the file, possibly reshaped 
		
		Raises
		------
		IOError : If unable to read from file.
		ValueError : If cannot reshape array into desired shape.
		"""
		pos=self.tell()
		try:
			if shapes is None:
				shapes=[None]*len(dtypes)
			return [self.readArray(dtype,shape) for dtype,shape in zip(dtypes,shapes)]
		except IOError, e:
			self.tell(pos)
			raise e
		except ValueError, e:
			self.seek(pos)
			raise e
	def writeArrays(self,arrays):
		"""
		Write a number of numpy arrays to file
		Parameters
		----------
		arrays : iterable of ndarrays, all to be written to file.

		Returns
		-------
		None

		Raises
		------
		IOError : If unable to write to file
		"""
		for data in arrays:
			self.writeArray(data)

	def advance(self,n=1):
		"""
		Advance a number of fields in the fortran file; skip them without reading them.
		Parameters
		----------
		n : optional, number of fields to advance (default = 1)
		
		Returns
		-------
		None
		
		Raises
		------
		IOError : If unable to advance in file or fortran format not correct.
		"""
		pos=self.tell()
		try:
			for i in xrange(n):
				nb=self._readSentinel()
				self.seek(nb,1) #seek from the current position
				self._checkSentinel(nb)
		except IOError, e:
			self.seek(pos)
			raise e
			

	def readArray(self,dtype,shape=None):
		"""
		Read a numpy array from a fortran file.
		If the read fails the file pointer remains where it started
		
		Parameters
		----------
		dtype : type, The numpy data type of the array to be read.
		shape : (optional) shape to reshape the read array into
		
		Returns
		-------
		data : ndarray, of data type dtype, read from the file, possibly reshaped
		
		Raises
		------
		IOError : If unable to read from file.
		ValueError : If cannot reshape array into desired shape.
		"""
		pos=self.tell()
		try:
			nb=self._readSentinel()
			print nb
			n = nb//np.nbytes[dtype]
			
			if n*np.nbytes[dtype]!=nb:
				raise IOError("Fortran array format not correct")
			data=np.fromfile(self,dtype,n)
			if not len(data)==n:
				raise IOError("Fortran array format not correct")
			self._checkSentinel(nb)
			if self.swap:
				data=data.byteswap()
			if shape is not None:
				data=data.reshape(shape)
			return data
		except IOError, err:
			self.seek(pos)
			raise err

	def writeArray(self,data):
		"""
		Write a numpy array to file
		Parameters
		----------
		data : ndarray, data to be written to file

		Returns
		-------
		None

		Raises
		------
		IOError : If unable to write to file
		"""
		data=np.array(data)
		nb=np.array(data.nbytes,dtype=self.sentinel)
		if self.swap:
			nb = nb.byteswap()
			data = data.byteswap()
		nb.tofile(self)
		data.tofile(self)
		nb.tofile(self)

