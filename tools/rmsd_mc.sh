#!/bin/bash
if [ -z "${VMDPATH+x}" ]; then
  echo "VMDPATH environment variable is not defined. It should be set to the vmd path. E.g  export VMDPATH=/Users/username/bin/vmd"
  exit
fi

MYVMD=$VMDPATH

usage()
{
  echo """

  rmsd_mc: tries to extract equally distant frames from a trajectory provided in pdb format.
  The input trajectory is in pdb format and assumed to be sequential (e.g. from a state to b state ) 
  and not stochastic.
  Steps performed: 
  	- builds an interpolation of the trajectory and is computing the total lenght.
	- given the number of output frames, it calculates the ideal length 
	- the frames are retrieved along the interpolation so that they are equally spaced along the interpolation
	- Since the equal distance along the interpolation doesn't mean that the frames are equally spaced among themselves
	  an optimization process is required
	      - the distance between adjacent frames is calculated and the frame is moved along the original interpolation
		    towards the farthest neighbor (note: we use the original interpolation otherwise the path will straighten
		    at every iteration). An empirical step (-s) is used for evolving
		  - since this is done for each frame this may give a non linear optimization so this needs to be iterated multiple
		    times 

  USAGE:
  ./rmsd_mc.sh  -h

  prints this help
  
  Example use:

  ./rmsd.sh  -f filename  -r 80 -n 30 -s 1 -m -t 50 -c 100 

  -f  filename.pdb		: input file to process
  -r 	80				: number of frames in output
  -n 	30 				: number of iterations in the optimization
  -s 	 1				: lenght of the step in the optimization
  -m 					: use monte carlo instead of steepest-descent like
  -t 	50				: read input file with a stride of 50 steps
  -d					: just dump the RMSD matrix and exit
  -c  	10				: checkpoint the reparametrized file REPARAM_i.pdb every 10 steps 

  defaults:  no mc ; -t 1 ; -s 1 

  """
}
rm -rf .l.tcl
if [ $# -eq 0 ]; then
   usage
   exit
fi
ndiv=0
stride=1
mc=0
nsteps=1
stepsize=1
justdump=0
cptidx=10
while getopts “dmhf:r:n:s:t:c:” OPTION
do
	case $OPTION in
	     h)
	            usage
	            exit 1
	     ;;
	     f)
	            filename=$(echo "$OPTARG" | sed 's/\./\\./g; s/\//\\\//g')
	     ;;
	     r)
		    ndiv=$OPTARG
	            echo "doing reparametrization in $ndiv pieces dumped on REPARAM.pdb" 
	     ;;	    
     	 t)
		    stride=$OPTARG
	            echo "reading the input file with stride $stride " 
	     ;;
     	 s)
		    stepsize=$OPTARG
	            echo "optimizing with a stepsize of $stepsize " 
	     ;;
     	 n)
		    nsteps=$OPTARG
	            echo "optimizing  in  $nsteps steps " 
		 ;;
     	 c)
		    cptidx=$OPTARG
	            echo "checkpointing the REPARAM every in  $cptidx steps " 
	     ;;
     	 m)
	            mc=1		
	            echo "doing the MC kind of optimization " 
	     ;;
	     d)
	            justdump=1		
	            echo "just dumping matrix and diagonal " 
	     ;;	
	     ?)
	             usage
	             exit
	     ;;
	esac	
done
#
# cannot do anything without n 
#
if [ "$ndiv" -eq "0"  ] ; then 
	if [ "$justdump" -eq "0"  ] ; then 
		echo "must specify the number of divisions with -r "
		usage
		exit
	fi
fi
cat >>.l.tcl<<"EOF"
#set name "small.pdb" ; puts -nonewline  ""
set nn I_NDIV
set stepsize I_STEPSIZE 
set nsteps I_NSTEPS 
set mc I_MC 
set stride I_STRIDE
set justdump I_JUSTDUMP
set cptidx I_CPTIDX
set name "I_FILENAME" ; puts -nonewline  ""
mol new $name type pdb first 0 last -1 step $stride filebonds 1 autobonds 1 waitfor all
# csv header
proc write_csv_header { fh v } {
	set mylist [list]
	for {set i 1} {$i <= [llength $v]} {incr i} {
    	lappend mylist "BASIS $i"
	}
	set head_basis [join  $mylist "," ]
	puts $fh "ENERGY,$head_basis"
}
proc write_csv_line { fh e v } {
	set head_basis [join  $v "," ]
	puts $fh "$e,$head_basis"
}



#
# reset com
#
proc reset_com { molid id } {
        set selocc [ atomselect $molid  "occupancy > 0." frame $id  ] ; puts -nonewline  ""
        set occval  [ $selocc get occupancy  ] ; puts -nonewline  "" 
        set occindex  [ $selocc get index  ] ; puts -nonewline  ""
	set gc [veczero]
    	set k 0
    	set totw 0. 
        foreach coord [ $selocc  get {x y z}] {
    		set myw [ lindex $occval $k ] 
    	        set gc [vecadd $gc [ vecscale $myw $coord ] ] 
    		set totw [ expr $totw+ $myw ]
    		incr k
        }
    	set gc  [vecscale [expr 1.0 /$totw ] $gc]
    	set newcoord []
        foreach coord [ [ atomselect $molid all  frame $id ] get {x y z}] {
    	    lappend newcoord [ vecsub $coord $gc ]	
    	}
    	[  atomselect $molid all frame $id ] set  { x y z } $newcoord 

}
#
# calculate msd and  distance vector given two frames: provides back   msd { position of id1-com} { delta (id2-com)-(id1-com)} 
#

proc calc_msd_and_dist { molid id1 id2 } {

	set plumedver 1 
	#puts "calculating msd and dist vector ... "
	#
	# save coordinates	
	#
        set store_id1 [ [ atomselect $molid  all frame $id1 ]  get { x y z } ] 
        set store_id2 [ [ atomselect $molid  all frame $id2 ]  get { x y z } ] 
	# 
	# center the com 
	# 
	#puts "reset com of $id1 and $id2 "	
	reset_com $molid $id1
	reset_com $molid $id2
	#puts "rescale the coordinate according the alignment weights "
        #puts "rescaling values ..."
	# save before the expansion but before reset of com
        set allcoor_id1 [ [ atomselect $molid  all frame $id1 ]  get { x y z } ] 
        set allcoor_id2 [ [ atomselect $molid  all frame $id2 ]  get { x y z } ] 
        #
        set selocc [ atomselect $molid  "occupancy > 0."  ] ; puts -nonewline  ""
        set occval  [ $selocc get occupancy  ] ; puts -nonewline  "" 
        set occindex  [ $selocc get index  ] ; puts -nonewline  ""
        set coord [ [  atomselect $molid "index $occindex" frame $id1 ] get { x y z } ] 
        for { set j 0 } { $j <  [llength $coord ] } { incr j } { 
                set myind [ lindex $occindex $j ]  ; puts -nonewline  ""
		if { $plumedver == 2 } { 
                	[  atomselect $molid "index $myind" frame $id1 ] set  { x y z }  [ list [ vecscale [ expr sqrt( [ lindex $occval $j ] ) ] [ lindex $coord $j ] ]    ] 
		} elseif { $plumedver == 1 } { 
                	[  atomselect $molid "index $myind" frame $id1 ] set  { x y z }  [ list [ vecscale  [ lindex $occval $j ]  [ lindex $coord $j ] ]    ] 
		}
        
        }
        set coord [ [  atomselect $molid "index $occindex" frame $id2 ] get { x y z } ] 
        for { set j 0 } { $j <  [llength $coord ] } { incr j } { 
                set myind [ lindex $occindex $j ]  ; puts -nonewline  ""
   	        if { $plumedver == 2 } {
             	   [  atomselect $molid "index $myind" frame $id2 ] set  { x y z }  [ list [ vecscale [ expr sqrt ( [ lindex $occval $j ] ) ]  [ lindex $coord $j ] ]    ] 
   	        } elseif { $plumedver == 1 } {
                   [  atomselect $molid "index $myind" frame $id2 ] set  { x y z }  [ list [ vecscale  [ lindex $occval $j ]   [ lindex $coord $j ] ]    ] 
		}
        }
       # now that it is scaled is time to calculate the alignment
       set mymat  [ measure fit [ atomselect $molid "index $occindex" frame $id2 ] [  atomselect $molid "index $occindex" frame $id1 ]  ]	
       # restore the value post resetcom	
       [ atomselect $molid all frame $id1 ] set  { x y z } [ lindex $allcoor_id1  ] ; puts -nonewline  ""
       [ atomselect $molid all frame $id2 ] set  { x y z } [ lindex $allcoor_id2  ] ; puts -nonewline  ""
       # now move id2 onto id1 
       [ atomselect $molid  all frame $id2 ] move $mymat 
       # now calculate the weighted msd
       set selbeta [ atomselect top  "beta > 0."  ] ; puts -nonewline  ""
       set betaval [ $selbeta get beta  ] ; puts -nonewline  ""
       set betaindex  [ $selbeta get index  ] ; puts -nonewline  ""
       set totw 0. 
       set totv 0. 
       for { set k 0 } { $k < [ llength $betaindex ] } { incr k } {
             set myind [ lindex $betaindex $k] ; puts -nonewline  ""	
             set ww [ lindex $betaval $k ]	 ; puts -nonewline  ""
             set pi [ lindex [ [  atomselect $molid "index $myind" frame $id2 ] get { x y z }   ] 0 ] ; puts -nonewline  ""
             set pj [ lindex [ [  atomselect $molid "index $myind" frame $id1 ] get { x y z }   ] 0 ]		 
	     if { $plumedver == 1 } {
	             set v [ expr [  veclength2 [ vecsub $pi $pj ] ] *  $ww * $ww ]  ; puts -nonewline  ""	
        	     set totw [ expr $totw + $ww*$ww ]  ; puts -nonewline  ""
             	     set totv [ expr $totv + [expr $v ] ]  ; puts -nonewline  ""
	     } elseif { $plumedver == 2 } {
	        set v [ expr [  veclength2 [ vecsub $pi $pj ] ]]  ; puts -nonewline  ""	
        	set totv [ expr $totv + ($v)*$ww  ]  ; puts -nonewline  ""
       	        set totw [ expr $totw + $ww ]  ; puts -nonewline  ""
	     }
             unset pi; unset pj ; unset v ; unset ww ; unset myind
       } 
       if { $plumedver == 1 } {
         # a-la-plumed 1
         set totv [ expr sqrt($totv / $totw) ]		
       } elseif { $plumedver == 2 } {
         set totv [ expr sqrt($totv / $totw) ]		
       }
       #puts "MSD between $id2  and $id1 : $totv "
       # 
       # give back also id1 and id2-id1 
       #
       set pi [ [  atomselect $molid all frame $id2 ] get { x y z }  ] ; puts -nonewline  "" 
       set pj [ [  atomselect $molid all frame $id1 ] get { x y z }  ]	; puts -nonewline  ""	 
       set delta [] ; puts -nonewline  ""
       for  { set k 0 } { $k < [ llength $pj ] } { incr k } {	
	 	lappend delta [ vecsub [ lindex $pi $k ] [ lindex $pj $k ] ]	; puts -nonewline  ""	
       }	
       return [ list  $totv $pj $delta  ] 
}


proc newscaledcoords { coords delta fact } {
        set newcoords []
	for { set i 0 } { $i < [llength $coords ] } { incr i } {
		lappend newcoords  [ vecadd [ lindex $coords $i ] [ vecscale $fact [ lindex $delta $i]  ] ]
	}	
	return $newcoords
}


proc test_rescaling { molid reference_coords reference_deltas  myindex myfactor  } {
	set deltas []	
	set bak1 0 
	set bak2 1 
	set save1 [ [ atomselect $molid all frame $bak1 ] get { x y z } ]
	set save2 [ [ atomselect $molid all frame $bak2 ] get { x y z } ]
	set harm 0.
	#puts "INPUT BASIS INDEX:  $myindex"
	#puts "INPUT BASIS FACTOR: $myfactor"
	
	for { set i  0 } { $i < [ expr [ llength $myindex ] -1 ] } { incr i } {
		# just pick two adjacent frames and save the coordinates
		#set bak1 [ lindex $myindex  $i  ]
		#set bak2 [expr [ lindex $myindex  $i  ] + 1]
		#set save1 [ [ atomselect $molid all frame $bak1 ] get { x y z } ]
		#set save2 [ [ atomselect $molid all frame $bak2 ] get { x y z } ]
	        #get info over the first point
		set id1 [ lindex $myindex  $i  ]
		set ff1 [lindex $myfactor $i ]
	        set cc1 [lindex $reference_coords $id1 ] 
		set dd1 [lindex $reference_deltas $id1 ]
	        #get info over the second point
		set id2 [ lindex $myindex [expr  $i +1 ]  ]
		set ff2 [lindex $myfactor [expr  $i +1 ]  ]
	        set cc2 [lindex $reference_coords $id2 ] 
		set dd2 [lindex $reference_deltas $id2 ]
		# now produce the two coordinate sets
		set c1 [ newscaledcoords $cc1 $dd1 $ff1  ]			
		set c2 [ newscaledcoords $cc2 $dd2 $ff2  ]			
		# replace them	
		[  atomselect $molid all frame $bak1 ] set  { x y z } $c1	
		[  atomselect $molid all frame $bak2 ] set  { x y z } $c2	
		# calculate the distance
	        set rr [ calc_msd_and_dist $molid $bak1 $bak2 ] ; puts -nonewline  "" 
		puts "DELTA from $i to [ expr  $i +1 ] :  [lindex $rr 0 ] "	
		lappend  deltas [lindex $rr 0 ];
		set harm [ expr $harm + [lindex $rr  0 ] **2  ]  ;
		# set back the coordinates
		#[  atomselect $molid all frame $bak1 ] set  { x y z } $save1	
		#[  atomselect $molid all frame $bak2 ] set  { x y z } $save2
	}
	#puts "HARMONIC_ENERGY $harm ";
	[  atomselect $molid all frame $bak1 ] set  { x y z } $save1	
	[  atomselect $molid all frame $bak2 ] set  { x y z } $save2
	return   [ list $deltas $harm]
}

proc move_images { step newdelta  myindex myfactor   } {
	#
	# remember: the deltas are one less than the frames
	#
	set mysum 0.
	for { set i  0 } { $i <  [ llength $newdelta ]  } { incr i } {
	  set mysum [ expr $mysum + [lindex $newdelta $i ] ] ; puts -nonewline  ""
	}
	set myavg [ expr $mysum / [llength $newdelta ] ];  puts -nonewline  ""
	#
	set prop [] ; puts -nonewline  "" 
	#
	# each delta connects frame $i with $i+1 
	# note that you move node i+1
	#
	for { set i  0 } { $i <  [expr  [ llength $newdelta ]  ] } { incr i } {
		set nf [ expr [lindex $myfactor [expr $i +1] ] *  ( $myavg / [lindex $newdelta $i ] )   ] ; 
		# alt nf calculation: use harmonic-like force
		# test 1 : nf is a force that is proportional to the distance between adjacent frames
		#set nf [ expr [lindex $myfactor [expr $i +1] ] * ( 1+ $step*( [lindex $newdelta [ expr $i+1 ]  ] -[lindex $newdelta $i ] )/$myavg )**3 ]	
		#set nf [ expr [lindex $myfactor [expr $i +1] ] * ( 1+ $step*( [lindex $newdelta [ expr $i+1 ]  ] -[lindex $newdelta $i ] )/$myavg )**1 ]	
		set nf [ expr [lindex $myfactor [expr $i +1] ] + $step*(( [lindex $newdelta [ expr $i+1 ]  ] -[lindex $newdelta $i ] )/$myavg )**1 ]	
		set perc_err [ expr sqrt(( $myavg - [lindex $newdelta $i ] )**2) / $myavg  ] ;
		puts "DISTANCE BEWTEEN  $i AND [expr $i +1]    AVG $myavg  PRESENT [lindex $newdelta $i ]  PROP [expr $myavg / [lindex $newdelta $i ] ]  PERCERR  $perc_err " 
		set newindex [lindex $myindex [expr $i +1] ]
		set newfactor $nf 
		if { $nf > 1. } {
			if { [expr [lindex $myindex [expr $i +1] ]+1] <= [ lindex $myindex end  ] } { 
				set newfactor [expr 0.01 ]
				set newindex [expr [lindex $myindex [expr $i +1] ]+1]
				puts "MOVING frame [expr $i +1]  on from snapshot [lindex $myindex [expr $i +1] ] to  adjancent snapshot : [expr [lindex $myindex [expr $i +1] ]+1]  "
			} else {
				set newfactor [expr 0.99 ]
				set newindex  [lindex $myindex [expr $i +1] ]
			}
		} elseif { $nf < 0. } {
			if { [expr [lindex $myindex [expr $i +1] ]-1] >=0 } { 
			puts "MOVING frame [expr $i +1]  on from snapshot [lindex $myindex [expr $i +1] ]  to  adjancent snapshot :  [expr [lindex $myindex [expr $i +1] ]-1]  "
			set newfactor [expr 0.99 ]
			set newindex  [expr [lindex $myindex [expr $i +1] ]-1] 
			} else {
			set newfactor [expr 0.01 ]
			set newindex  [lindex $myindex [expr $i +1] ]
			}
		}
		set nextindex [expr $i +2] 
		set previndex $i 
		set nextfactor [lindex $myfactor [expr $i +2] ] 
		set prevfactor [lindex $myfactor  $i  ] 
		
		lset myfactor [expr $i +1]  $newfactor     
		lset myindex  [expr $i +1]  $newindex 
	
		#lappend prop  [ expr ( $myavg / [lindex $newdelta $i ] ) ]  ; 	
		lappend prop  $perc_err  ; 	
	}
	return [ list $myindex $myfactor  $prop ]
}

#
# this gives a weight in term of sum of harmonics 
#
proc harmonic_potential_and_force { molid reference_coords reference_deltas  myindex myfactor } {

	set res [ test_rescaling  $molid $reference_coords $reference_deltas  $myindex $myfactor ] 
	set delta [ lindex $res 0]
	set energy [ lindex $res 1]
	set harm 0. ; 
	for { set i 0 } { $i < [ llength $delta ]} { incr i } {
		set harm [ expr $harm + [lindex $delta  $i ] **2  ]  ;
	}
	set force []  
	for { set i 0 } { $i < [ llength $delta ] } { incr i } {
		lappend force 0. ; 
	}	
	lappend force 0. ; # now thw force is equal to the number of images
	for { set i 0 } { $i < [ llength $delta ]  } { incr i } {
		set oldval [  lindex $force [expr $i + 1] ] ; 
		set newval  [ expr $oldval - [lindex $delta $i ] ] ;		
		lset force [expr $i + 1] $newval   ;
		set oldval [  lindex $force $i  ] ; 
		set newval  [ expr $oldval + [lindex $delta $i ] ] ;		
		lset force $i $newval   ;
	}
	#puts "FORCE $force ";
	# set first and last with zero force
	lset force 0 0.   ;
	lset force [ llength $delta ]   0.   ;
	return [ list $harm $force ] ;
}
proc eval_lambda {  molid reference_coords reference_deltas  myindex myfactor } {
	set res [ test_rescaling  $molid $reference_coords $reference_deltas  $myindex $myfactor ] 
	set delta [ lindex $res 0]
	set energy [ lindex $res 1]
	set avg 0. ; 
	set devstd 0. ; 
	for { set i 0 } { $i < [ llength $delta ]} { incr i } {
		set avg [ expr $avg + [lindex $delta  $i ] **2  ]  ;
		set devstd [ expr $devstd + [lindex $delta  $i ] **4  ]  ;
	}
	set avg [ expr $avg / [ llength $delta ] ] 
	set devstd [ expr $devstd / [ llength $delta ] ] 
	set devstd [  expr sqrt( $devstd -$avg*$avg ) ] 
	puts "******************************"
	puts ">>>> MSD        $avg  Ang^2 " 
	puts ">>>> MSD_STDEV  $devstd  Ang^2 " 
	set safe [ expr 2.3/( $devstd+$avg  ) ] ;puts -nonewline "" 
	set lambda [ expr  2.3/$avg  ]  ;puts -nonewline "" 
	puts ">>>> LAMBDA $lambda  \[Ang ^ -2\] "
	puts ">>>> SAFE_LIMIT $safe \[Ang ^ -2\] "
	puts ">>>> Note: if you use gromacs you have to increase these values of 100 since units are \[ nm ^-2 \]"
	puts ">>>> Note2: LAMBDA and SAFE_LIMIT are truly different your matrix could be messed up: check MSD_MATRIX FILE "
	puts "******************************"
}

proc eval_matrix {  molid reference_coords reference_deltas  myindex myfactor fname } {
	#array set arr {}
	set bak1 0
	set bak2 1 
	set save1 [ [ atomselect $molid all frame $bak1 ] get { x y z } ]
	set save2 [ [ atomselect $molid all frame $bak2 ] get { x y z } ]
	for { set i  0  } { $i < [ llength $myindex ] } { incr i } {
		puts "MATRIX_CALCULATION element $i "	
		for { set j $i } { $j < [ llength $myindex ] } { incr j } {
			# just pick two adjacent frames and save the coordinates
	        	#get info over the first point
			set id1 [ lindex $myindex  $i  ]
			set ff1 [lindex $myfactor $i ]
	        	set cc1 [lindex $reference_coords $id1 ] 
			set dd1 [lindex $reference_deltas $id1 ]
	        	#get info over the second point
			set id2 [ lindex $myindex $j   ]
			set ff2 [lindex $myfactor $j   ]
	        	set cc2 [lindex $reference_coords $id2 ] 
			set dd2 [lindex $reference_deltas $id2 ]
			# now produce the two coordinate sets
			set c1 [ newscaledcoords $cc1 $dd1 $ff1  ]			
			set c2 [ newscaledcoords $cc2 $dd2 $ff2  ]			
			# replace them	
			[  atomselect $molid all frame $bak1 ] set  { x y z } $c1	
			[  atomselect $molid all frame $bak2 ] set  { x y z } $c2	
			# calculate the distance
	        	set rr [ calc_msd_and_dist $molid $bak1 $bak2 ] ; puts -nonewline  "" 
			#puts "DELTA from $i to $j :  [lindex $rr 0 ] "	
			set arr($i,$j) [lindex $rr 0 ]
			set arr($j,$i) [lindex $rr 0 ]
			# set back the coordinates
		}
	}
	[  atomselect $molid all frame $bak1 ] set  { x y z } $save1	
	[  atomselect $molid all frame $bak2 ] set  { x y z } $save2
   	set fp [ open $fname w ] ; puts -nonewline  ""
	for { set i  0  } { $i < [ llength $myindex ] } { incr i } {
		for { set j 0 } { $j < [ llength $myindex ] } { incr j } {
		   puts $fp "[expr $i+1] [expr $j+1]  $arr($i,$j)"
		}
   		puts $fp  "  "
   	}
	close $fp
}
proc eval_diag {  molid reference_coords reference_deltas  myindex myfactor fname } {
	#array set arr {}
	set bak1 0
	set bak2 1 
	set save1 [ [ atomselect $molid all frame $bak1 ] get { x y z } ]
	set save2 [ [ atomselect $molid all frame $bak2 ] get { x y z } ]
	for { set i  0  } { $i < [ expr [ llength $myindex ] -1 ]  } { incr i } {
		#puts "DIAGONAL_CALCULATION element $i "	
			# just pick two adjacent frames and save the coordinates
	        	#get info over the first point
			set id1 [ lindex $myindex  $i  ]
			set ff1 [lindex $myfactor $i ]
	        	set cc1 [lindex $reference_coords $id1 ] 
			set dd1 [lindex $reference_deltas $id1 ]
	        	#get info over the second point
			set id2 [ lindex $myindex [expr $i+1]   ]
			set ff2 [lindex $myfactor [expr $i+1]   ]
	        	set cc2 [lindex $reference_coords $id2 ] 
			set dd2 [lindex $reference_deltas $id2 ]
			# now produce the two coordinate sets
			set c1 [ newscaledcoords $cc1 $dd1 $ff1  ]			
			set c2 [ newscaledcoords $cc2 $dd2 $ff2  ]			
			# replace them	
			[  atomselect $molid all frame $bak1 ] set  { x y z } $c1	
			[  atomselect $molid all frame $bak2 ] set  { x y z } $c2	
			# calculate the distance
	        	set rr [ calc_msd_and_dist $molid $bak1 $bak2 ] ; puts -nonewline  "" 
			#puts "DELTA from $i to $j :  [lindex $rr 0 ] "	
			set arr($i) [lindex $rr 0 ]
			# set back the coordinates
	}
	[  atomselect $molid all frame $bak1 ] set  { x y z } $save1	
	[  atomselect $molid all frame $bak2 ] set  { x y z } $save2
   	set fp [ open $fname w ] ; puts -nonewline  ""
	for { set i  0  } { $i < [ expr [ llength $myindex ] -1 ] } { incr i } {
	   puts $fp "[expr $i+1]   $arr($i)"
   	}
	close $fp
}



proc mcmove { step force myindex myfactor  } {
	#
	# each delta connects frame $i with $i+1 
	# note that you move node i+1
	#
	#
	for { set i  0 } { $i <   [ llength $force ]  } { incr i } {
		set nf [  expr   [lindex $myfactor $i ] *(1.+$step*rand()* [lindex $force $i ] )    ] ; 
		# alt nf calculation: use harmonic-like force
		if { $nf > 1. } {
			puts "MOVING frame  $i  on adjancent frame : [expr [lindex $myindex $i  ]+1]  "
			set nf [expr 0.01 ]
			lset myfactor  $i   $nf     
			lset myindex   $i   [expr [lindex $myindex $i ]+1]    
		} elseif { $nf < 0. } {
			puts "MOVING frame  $i  on adjancent frame : [expr [lindex $myindex $i ]-1]  "
			set nf [expr 0.99 ]
			lset myfactor  $i   $nf     
			lset myindex   $i   [expr [lindex $myindex $i ]-1]    
		} else {
			lset myfactor  $i $nf   
		}
		#lappend prop  [ expr ( $myavg / [lindex $newdelta $i ] ) ]  ; 	
	}
	return [ list $myindex $myfactor  ]
}

proc dump_frames { molid reference_coords reference_deltas  myindex myfactor cptidx } {
   set mydir [ exec  pwd  ];
   exec rm -rf REPARAM_${cptidx}.pdb .ff_* 
   for { set i  0 } { $i < [llength $myindex ] } { incr i } {
	set id [ lindex $myindex $i ]
	set newcoord [ newscaledcoords [ lindex $reference_coords $id ]  [ lindex $reference_deltas $id ] [ lindex $myfactor $i ] ] 
   	[  atomselect $molid all frame 0  ] set  { x y z }  $newcoord  
   	[  atomselect $molid all frame 0  ] writepdb ".ff_$i.pdb"
   	exec cat .ff_$i.pdb >>REPARAM_${cptidx}.pdb
   	exec rm -rf .ff_$i.pdb 
   }
}

#diagonal distance 
set nframes [ molinfo top get numframes ] ; puts -nonewline  ""
set cumlength []
set reference_coords []
set reference_deltas []
set totl 0.
set myfactor [ list 0. ]
set myindex [ list 0 ]
# first element
#set rr [ calc_msd_and_dist 0  0 0  ] ; puts -nonewline  "" 
#puts "TEST: ZERO DISTANCE FRAME IS [ lindex $rr 0 ] "
#lappend reference_coords [ lindex $rr 1 ] ; puts -nonewline  ""	
#lappend reference_deltas [ lindex $rr 2 ]	; puts -nonewline  ""
lappend cumlength 0. ; puts -nonewline  ""
lappend myindex  0 ; puts -nonewline  ""
lappend myfactor 0. 	; puts -nonewline  ""
for  { set i 0 } { $i <  [ expr $nframes -1 ]  } { incr i } {
   set rr [ calc_msd_and_dist 0 $i [ expr $i +1 ] ] ; puts -nonewline  "" 
   set totl [ expr $totl  + [ lindex $rr 0 ] ] 
   lappend reference_coords [ lindex $rr 1 ]	
   lappend reference_deltas [ lindex $rr 2 ]	
   lappend cumlength $totl	
   lappend myindex [expr $i  ]	
   lappend myfactor 1.	
   #puts "DISTANCE FROM FRAME $i and [ expr $i +1  ] is [ lindex $rr 0 ] "
}
if { $justdump > 0 } {
	eval_matrix  0 $reference_coords $reference_deltas  $myindex $myfactor "RMSD_MATRIX" ;  puts -nonewline  ""
	eval_diag  0 $reference_coords $reference_deltas  $myindex $myfactor "RMSD_DIAGONAL"  ;  puts -nonewline  ""
	eval_lambda  0 $reference_coords $reference_deltas  $myindex $myfactor   ;  puts -nonewline  ""
	puts "BYE..."
	quit
}
puts $totl
puts $cumlength
set step [ expr $totl / ( $nn -1 ) ]  ; puts -nonewline  ""
# now do the maragliano thing
# first image is the first index with 0 rescaling
set myfactor [ list 0. ] ;
set myindex [ list 0 ] ;
for { set i  1 } { $i < $nn } { incr i } {
	set myl [ expr $step * $i  ]
        for { set j  1 } { $j < [llength $cumlength] } { incr j } {
		if { $myl < [lindex $cumlength $j ] && $myl > [lindex $cumlength [expr  $j -1] ]  } {
			set s2 [lindex $cumlength $j ] 
			set s1 [lindex $cumlength [expr  $j -1] ]  
			set fact [ expr ($myl-$s1) / ($s2-$s1)   ] 
			lappend myfactor $fact		;	
			lappend myindex [expr  $j -1 ] ;
			puts "FRAME $i length $myl is between  [lindex $cumlength [expr  $j -1] ] ( block [expr  $j -1] )  and [lindex $cumlength $j ]   ( block [expr  $j ] ) : fact $fact"
		}
 	}	
} 
# last image is factor 1 with respect of last-but-one image
lappend myfactor  1. ;  puts -nonewline  "" 
lappend myindex [ expr $nframes -2 ];  puts -nonewline  ""
# 

set filename "output.csv"
set fileout [open $filename w]
# header
write_csv_header $fileout $myindex

puts "BASIS INDEX:  $myindex" ;  puts -nonewline  ""
puts "BASIS FACTOR: $myfactor" ;  puts -nonewline  ""

#
# do the mc evolution
#
if { $mc == 1} {
	set pot_force [ harmonic_potential_and_force 0 $reference_coords $reference_deltas  $myindex $myfactor ] 
	set oldpot [ lindex $pot_force 0 ]
	set oldforce [ lindex $pot_force 1 ]
	puts "ENTERING MCLOOP OLDPOT: $oldpot "
	for { set j 0  } { $j < $nsteps } { incr j } {
		  set new_ind_and_factor [ mcmove $stepsize  $oldforce $myindex $myfactor ]
		  set newindex 	 [lindex $new_ind_and_factor 0 ]
		  set newfactor 	 [lindex $new_ind_and_factor 1 ]
		  set pot_force [ harmonic_potential_and_force 0 $reference_coords $reference_deltas  $newindex $newfactor ] 
		  #eval_diag  0 $reference_coords $reference_deltas  $myindex $myfactor "RMSD_DIAGONAL_$j"   ;  puts -nonewline  ""
		  set pot [ lindex $pot_force 0 ]
		  set force [ lindex $pot_force 1 ]
		  if { $oldpot > $pot } {
			puts "TEST ACCEPTED !!! OLD $oldpot > NEW $pot "
			set myindex $newindex ;
			set myfactor $newfactor ;
			set oldpot $pot 
			set oldforce $force 
		  } else {
			puts "TEST REJECTED !!!  OLD $oldpot < NEW $pot "
		  }	
		  if {[expr {$j % $cptidx}] == 0 } {
		  	dump_frames 0 $reference_coords $reference_deltas  $myindex $myfactor $j
		  }
	}
	dump_frames 0 $reference_coords $reference_deltas  $myindex $myfactor "FINAL"
} else {
	for { set j 0  } { $j < $nsteps } { incr j } {
		# calculate distance
		set res [ test_rescaling 0 $reference_coords $reference_deltas  $myindex $myfactor ] ;  puts -nonewline  ""
		set newdelta [ lindex $res 0]
		set energy [ lindex $res 1]
		write_csv_line $fileout $energy $myindex
		# save matrix
		#eval_matrix  0 $reference_coords $reference_deltas  $myindex $myfactor "RMSD_MATRIX_$j"  ;  puts -nonewline  ""
		set harm [ eval_diag  0 $reference_coords $reference_deltas  $myindex $myfactor "RMSD_DIAGONAL_$j" ]  ;  puts -nonewline  ""
		# move images
		set mm [ move_images $stepsize  $newdelta $myindex $myfactor  ]
		set newindex [ lindex $mm 0 ] ;  puts -nonewline  ""
		set newfactor [ lindex $mm 1 ] ;  puts -nonewline  ""
		set prop [ lindex $mm 2 ] ;  puts -nonewline  ""
		#set go 0 
		#for { set k 0  } { $k < [ llength $prop ] } { incr k } {
		#	set val [ lindex $prop $k ] ; 
		#	if  { $val <0.1 } { set go 1 } 
		#}
		puts "REPLACING WITH NEW FRAMES" ; 
		# remember that the last one should stay fixed
		for { set k 0  } { $k < [ expr [ llength  $prop ]  ] } { incr k } {
	        	 lset myfactor $k [ lindex $newfactor $k ]; 
	        	 lset myindex  $k [ lindex $newindex  $k ]; 
		}
		#if { $go != 1 } {
		#	puts "CYCLE DONE!" ; 
		#	break
		#}
		# report
		puts "BASIS INDEX:  $myindex"
		puts "BASIS FACTOR: $myfactor"
		if {[expr {$j % $cptidx}] == 0 } {
		  	dump_frames 0 $reference_coords $reference_deltas  $myindex $myfactor $j
		}
		
	}
	dump_frames 0 $reference_coords $reference_deltas  $myindex $myfactor "FINAL"
}
#
# now calculate the lambda
#
#
#
set fname "RMSD_MATRIX_FINAL"
#eval_matrix  0 $reference_coords $reference_deltas  $myindex $myfactor $fname ;  puts -nonewline  ""
eval_diag  0 $reference_coords $reference_deltas  $myindex $myfactor "RMSD_DIAGONAL_FINAL"  ;  puts -nonewline  ""
eval_lambda  0 $reference_coords $reference_deltas  $myindex $myfactor   ;  puts -nonewline  ""
close $fileout
quit
EOF
sed  s/I_NDIV/$ndiv/ .l.tcl | sed s/I_STEPSIZE/$stepsize/ | sed s/I_NSTEPS/$nsteps/ | sed s/I_MC/$mc/ | sed s/I_STRIDE/$stride/ | sed s/I_FILENAME/$filename/ | sed s/I_JUSTDUMP/$justdump/ | sed s/I_CPTIDX/$cptidx/ > .m.tcl; mv .m.tcl .l.tcl 
$MYVMD -dispdev text -e .l.tcl
exit
# now clean the resulting REPARAM.pdb
if [ $justdump -eq 0 ] ; then
cat >.p.pl <<\EOF
#!/usr/bin/perl
my $ifile=shift;
my @index;
open (FILE, "<  $ifile " );
while (<FILE>){
        my @a=split("",$_);
        my $atom=join('',@a[0..5]);
	my $ind=join('',@a[6..10]);
	my $end=join('',@a[0..2]);
	if($end=~/END|TER/){last ;};
	if($atom=~/ATOM|HETATM/){
		$ind=~s/^s+|s+$//;
		push @index, $ind;
	}
}
close FILE;
# paste into the new file
my $ifile=shift;
open(FILE,"<$ifile");
my $i;
while(<FILE>){
	my @a=split("",$_);
	my $atom=join('',@a[0..5]);
	my $back=join('',@a[11..65]);
	if($atom=~/ATOM|HETATM/){
			my $ind=sprintf("%5d",@index[$i]); $i++	;	
			print $atom.$ind.$back."\n";
	}else{print $_;$i=0;}
}
close FILE;
EOF
chmod +x .p.pl ; ./.p.pl $filename REPARAM_FINAL.pdb |grep -v CRYS >.l ; mv .l REPARAM_FINAL.pdb ; rm -rf .p.pl
fi
