proc morph { N } {

        set nframes [molinfo top get numframes]
        set molid [molinfo top get id]

        for {set i 0} {$i < [expr $nframes -1]} {incr i} {
                puts "$i, [expr $i +1]"

                for {set j 1} {$j <= $N} {incr j} {animate dup frame $i $molid}
        }

        animate dup frame $i $molid

        animate delete beg 0 end [expr $nframes -1] $molid

        for {set i 0} {$i < [expr $nframes -1]} {incr i} {
                set sel1 [atomselect $molid "all" frame [expr $N*$i]]
                set sel2 [atomselect $molid "all" frame [expr $N*$i+$N]]
                set x1 [$sel1 get x]
                set y1 [$sel1 get y]
                set z1 [$sel1 get z]
                set x2 [$sel2 get x]
                set y2 [$sel2 get y]
                set z2 [$sel2 get z]

                puts "[lindex $x1 0], [lindex $x2 0]"

                for {set j 0} {$j < $N} {incr j} {
                        set f [expr {double($j) / double($N)}]
                        puts "Scale frame [expr $N*$i+$j] from [expr $N*$i] to [expr $N*$i+$N] by $f"
                        $sel1 frame [expr $N*$i+$j]
                        $sel1 set x [vecadd [vecscale [expr {1.0 - $f}] $x1] [vecscale $f $x2]]
                        $sel1 set y [vecadd [vecscale [expr {1.0 - $f}] $y1] [vecscale $f $y2]]
                        $sel1 set z [vecadd [vecscale [expr {1.0 - $f}] $z1] [vecscale $f $z2]]

                }
        }

}