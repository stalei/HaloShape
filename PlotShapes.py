#  © Shahram Talei @ 2020 The University of Alabama - All rights reserved.
#you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation; either version 3 of the License, or
#(at your option) any later version.
#You should have received a copy of the GNU General Public License
#along with this program.  If not, see <http://www.gnu.org/licenses/>.


#files
fG='/media/shahram/SD/19880/g3z19880-976.0.csv'
fC='/media/shahram/SD/19880/c3z19880-932.0.csv'
hID=19880

#open and plot
import csv
with open(fG, newline='') as csvfile:
    reader = csv.reader(csvfile, delimiter=' ', quotechar='|')
    next(reader)
    next(reader)
    next(reader)
    for row in reader:
        print(row[1])
