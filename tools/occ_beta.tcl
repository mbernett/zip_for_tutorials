
set path_value $env(PATHTOOLS)
puts "PATHTOOLS in $path_value"
#source $path_value/my_morph.tcl

# parse arguments
set options [dict create]
for {set i 0} {$i < [llength $argv]} {incr i 2} {
  set key [lindex $argv $i]
  set value [lindex $argv [expr $i + 1]]
  puts "KEY $key VALUE $value "
  dict set options $key $value
}

set pdb [dict get $options "-pdb"]
set traj [dict get $options "-traj"]
set protsel [dict get $options "-protsel"]
set ligsel [dict get $options "-ligsel"]
set opdb [dict get $options "-outpdb"]

# workaround for reading a selection
set  ligsel [join [split $ligsel "__"] " "]
set  protsel [join [split $protsel "__"] " "]

#
puts "PDB $pdb "
puts "TRAJ" $traj
puts "PROTSEL $protsel"
puts "LIGSEL $ligsel"
puts "OUTPDB $opdb"

# VMD:
#mol new $pdb waitfor all
#morph 5

mol new $pdb
mol addfile $traj last -1 step 1 waitfor all

set all [atomselect top "all"]
set sel_prot [atomselect top $protsel ]
set sel_lig [atomselect top $ligsel ]
$all set occupancy 0
$sel_prot set occupancy 1
$sel_lig set beta 1

animate write pdb $opdb beg 0 end -1 waitfor all

quit
