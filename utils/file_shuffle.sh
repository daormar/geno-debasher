# Bio-PanPipe package
# Copyright (C) 2019,2020 Daniel Ortiz-Mart\'inez
#  
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public License
# as published by the Free Software Foundation; either version 3
# of the License, or (at your option) any later version.
#  
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#  
# You should have received a copy of the GNU Lesser General Public License
# along with this program; If not, see <http://www.gnu.org/licenses/>.
  
# *- bash -*

# Randomly reorders the lines of a file given a seed value

########
shuffle()
{
    # Initialize variables
    local seed=$1
    local tdir=$2
    local file=$3

    # Shuffle file
    "$AWK" -v seed=${seed} 'BEGIN{srand(seed)}{printf"%f %d %s\n",rand(),NR,$0}' "${file}" \
        | LC_ALL=C "$SORT" -k1n -k2n -T "${tdir}" | "${CUT}" -d' ' -f3-
}

########
shuffle_alt()
{
    # Alternative implementation (it has a higher spatial complexity)

    # Initialize variables
    local seed=$1
    local file=$2
    
    # Shuffle file
    "$AWK" -v seed=$seed \
        'function random(b) {return rand()*b}
  BEGIN{
         i=0
       }
       {
         # Store line of file
         line[i]=$0
         i=i+1
       }
    END{
         # Generate random numbers
         num_lines=NR
         srand(seed)
         for(i=0;i<num_lines;++i) 
         {
           vec_lines[i]=i
         }

         # Determine new ordering for the lines
         j=num_lines-1
         for(i=0;i<num_lines;++i) 
         {
           # Obtain next line number
           n=int(random(j+1))
           next_line_num=vec_lines[n]
           vec_lines[n]=vec_lines[j]
           j=j-1

           # Print line
           printf"%s\n",line[next_line_num]
         }

       }' "$file"
}

########
if [ $# -ne 2 -a $# -ne 3 ]; then
    echo "Usage: file_shuffle <seed> <tmpdir> [<file>]"
else

    # Take parameters
    seed=$1
    tdir=$2

    if [ $# -eq 3 ]; then
        file=$3
    fi

    # Invoke shuffling function
    shuffle $seed "$tdir" "$file"

fi
