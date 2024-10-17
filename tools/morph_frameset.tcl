
set path_value $env(PATHTOOLS)
puts "PATHTOOLS in $path_value"
source $path_value/morphing.tcl

# parse arguments
set options [dict create]
for {set i 0} {$i < [llength $argv]} {incr i 2} {
  set key [lindex $argv $i]
  set value [lindex $argv [expr $i + 1]]
  puts "KEY $key VALUE $value "
  dict set options $key $value
}

set pdb [dict get $options "-pdb"]
set opdb [dict get $options "-outpdb"]
set morphN [dict get $options "-morph"]

#
puts "PDB $pdb "
puts "OUTPDB $opdb"

# VMD:

mol new $pdb waitfor all

morph $morphN

animate write pdb $opdb beg 0 end -1 waitfor all

quit
