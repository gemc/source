#!/usr/bin/env python3

import sys
import sqlite3


def compare_with_database(line, cursor, system, runno, variation, nlines, index):
	# Split the line into fields
	fields = line.split('|')

	# Extract values from the line
	name, mother, dummy1, pos, rot, color, shape, shape_params, material, magfield, dummy2, dummy3, exist, visible, style, sensitivity, hitType, identity = [
		f.strip() for f in fields]

	# Query the database using the extracted values
	query = f"SELECT * FROM geometry WHERE name=? AND mother=? AND pos=? AND rot=? AND col=? AND type=? AND dimensions=? AND material=? and magfield=? and exist=? and visible=? and style=? and sensitivity=? and hitType=? and identity=? and system=? and run=? and variation=?"
	cursor.execute(query, (
	name, mother, pos, rot, color, shape, shape_params, material, magfield, exist, visible, style,
	sensitivity, hitType, identity, system, runno, variation))

	# Fetch the result
	result = cursor.fetchone()

	# if the fetch is None, then the line is not in the database. Print the line index with nlines.
	if result is None:
		print(f"Line {index + 1} out of {nlines} is NOT in the database.", end='\r')

	# Check if the result is not None, indicating a match
	return result is not None


def main():
	file_path = sys.argv[1]
	database_path = sys.argv[2]
	system = sys.argv[3]
	runno = sys.argv[4]
	variation = sys.argv[5]

	table_name = 'geometry'

	try:
		# Connect to the database
		connection = sqlite3.connect(database_path)
		cursor = connection.cursor()
		nlines = 0

		with open(file_path, 'r') as file:
			lines = file.readlines()

		# check that the number of lines in the file is the same as the number of lines in the database
		query = f"SELECT COUNT(*) FROM {table_name} WHERE system=? and run=? and variation=?"
		cursor.execute(query, (system, runno, variation))
		result = cursor.fetchone()
		if result is not None:
			nlines = result[0]
			if nlines != len(lines):
				print(
					f"\nNumber of lines in the file ({len(lines)}) does not match the number of lines in the database ({result[0]}).")
				exit (1)
			else:
				print(
					f"\nNumber of lines in the file ({len(lines)}) matches the number of lines in the database ({result[0]}).")

		# Compare each line with the database
		match=1
		for i, line in enumerate(lines):
			if not compare_with_database(line, cursor, system, runno, variation, nlines, i):
				print(f"Line {i + 1} does not match the database.")
				match=0
				print(line)
				exit (1)
		if match==1:
			print("\nAll lines match the database.")


	except FileNotFoundError:
		print(f"File '{file_path}' not found.")
	except sqlite3.Error as e:
		print(f"SQLite error: {e}")
	finally:
		# Close the database connection
		if connection:
			connection.close()


if __name__ == "__main__":
	main()
