package main


import (
    "fmt"
    "math"
    "math/rand"
    "time"
    "bufio"
    "os"
    "strings"
    "strconv"
)


type Particle struct {
    id int
    x_pos float64
    y_pos float64
    z_pos float64
    ux    float64
    uy    float64
    uz    float64
}

// UNIVERSALS
const Pi float64 = math.Pi

// SIMULATION SPECIFIC GLOBALS
var STEP_LIMIT int = 1200
var asp float64 = 10         // instead of D!
var boxL float64 = 50       // this is the TOTAL box length in one direction
var phi_target float64 = .39
var vol_step_size float64 = .000037   // controls scale of volume growth, V_new = (1+vol_step_size)*V_old
var trans_mag float64 = .96         // value from paper, supposed highest density packings
var rot_min float64 = .001*Pi/3.0
var rot_max float64 = Pi/55.0       // POSSIBLE: lowering helps, as we are so ordered that we mostly want to translate
                                    //           big rotations from random selection actually cause MORE overlaps!

func main() {
    // make sure we the right amount of arguments
    if len(os.Args) != 3 {
        fmt.Println("Not enough arguments!")
        os.Exit(1)
    }
    // need to intialize random num generator
    rand.Seed(time.Now().UnixNano())

    // check setup config
    fmt.Println("Particle apsect ratio: ", asp)
    fmt.Println("Box volume: ", math.Pow(boxL,3))
    fmt.Println("min rand ang, max rand ang: ", rot_min, rot_max)
   
    // read in particle data, grabs r_com and u_hat
    parlist := read_particles(os.Args[1])
    fmt.Println("READ FROM: ", os.Args[1])
    // this is the total particle length, not just the cylindrical part
    L_0, _ := strconv.ParseFloat(strings.Split(os.Args[1], "_")[3], 64)

    // str for unique save names!
    filenum := string(strings.Split(os.Args[1], "_")[1][4:])

    fmt.Println("Initial particle length: ", L_0)
    writefile := false
    if (os.Args[2] == "yes" || os.Args[2] == "y") {
        writefile = true
    }

    // structure: check phi, loop in steps increasing vol, run full algo
    phi_curr := float64(len(parlist))*vol_p(L_0, L_0/asp) / (math.Pow(boxL,3))
    fmt.Println("volume fraction initial: ", phi_curr)
   
    // set percentage to increase volume in one step
    scale := 1 + vol_step_size
    length_scale := math.Pow(1 + 4.5*vol_step_size, 1/3.0)   // larger initial that may (will!) be cut
    min_l_scale  := math.Pow(scale, 1/3.0)

    // set num steps... -> phi_targ/phi_curr = full scaling, divide this by size to get steps!
    num_steps := int(math.Ceil(math.Log(phi_target/phi_curr)/math.Log(scale)))  // ceil stupidly returns float!
    fmt.Println("needed steps to reach target: ", num_steps)
    L_curr := L_0
    m := make(map[int][][3]float64)

    // time the volume expansion loop
    start := time.Now()

    // initialize contact map
    contact_m := init_contact_map(len(parlist))

    // flag to write avg, activated when there ever were any contacts at a particle volume loop
    avg_flag := false

    // initalize nearest neighbor list
    near_n := get_nearest(parlist, L_curr)

    // label for the vol scaling loop
    vol_loop:
    for vol_step := 0; vol_step < num_steps; vol_step ++ {   // this steps us along in volume!
        //loop_start := time.Now()
        if int(math.Mod(float64(vol_step), 200.0)) == 0 {
            fmt.Println("On vol scaling step = ", vol_step)
            fmt.Println("vol_frac_curr = ", float64(len(parlist))*vol_p(L_curr, L_curr/asp) / (math.Pow(boxL,3)) )
            fmt.Println("took so far: ", time.Since(start))
        }
        L_curr = length_scale*L_curr
        num_unoverlap_steps := 0         // used to limit number of times adjusting overlaps
        // unoverlapping loop, if first step, or contacts before, and under limit of steps
        for ((len(m) > 0 && num_unoverlap_steps < STEP_LIMIT) || num_unoverlap_steps == 0) {
            m = make(map[int][][3]float64)  // resetting m initially
            if num_unoverlap_steps <= 100 {
                for i := 0; i < len(parlist); i ++ {
                    m, contact_m = check_overlaps(parlist[i], near_n[parlist[i].id], L_curr, m, contact_m)
                }
            } else { // otherwise, switch back to N^2 loops... these are faster than many N^2 within each vol?
                for i := 0; i < len(parlist); i ++ {
                    m, contact_m = check_overlaps(parlist[i], parlist, L_curr, m, contact_m)
                }
            }
            // if there were any contacts, flag this (to save data later)
            if len(m) > 0 {
                parlist = update_particles(parlist, m)
                // re-construct near_n_list since things moved?
                // only do it a couple times OR if many reorganizations, every so often? (presumption: not needed every rearrange)
                if (num_unoverlap_steps > 4 && num_unoverlap_steps <= 100) {
                    near_n = get_nearest(parlist, L_curr)
                }
                avg_flag = true  // could put above?
            }
            num_unoverlap_steps = num_unoverlap_steps + 1   // add frst so we dont print on 0 updates!
            if int(math.Mod(float64(num_unoverlap_steps), 100.0)) == 0 {
                fmt.Println("unoverlap tries = ", num_unoverlap_steps)
                fmt.Println("vol_frac_curr = ", float64(len(parlist))*vol_p(L_curr, L_curr/asp) / (math.Pow(boxL,3)) )
                fmt.Println("")
            }
        }
        // in case of 30 (pick num..) unoverlap, cut lengthscale until at min!
        // only need to do if not already at min
        if num_unoverlap_steps >= 40 && length_scale > min_l_scale {
            //fmt.Println("Had to cut length scaling...")
            //fmt.Println("on phi = ", float64(len(parlist))*vol_p(L_curr, L_curr/asp) / (math.Pow(boxL,3)) )
            if (length_scale - 1) / 2 >= (min_l_scale - 1) {
                length_scale = ((length_scale - 1) / 2) + 1
            } else {
                length_scale = min_l_scale
            }
            //fmt.Println("Curr length scale, min: ", length_scale, min_l_scale)
        }
        //fmt.Println("did unoverlapping for L = ", L_curr)
        // write averaged contact data to file
        if avg_flag {
            // so write if we want to at all, its later in the simulation, and there were even any contacts that were undone at all
            if writefile && float64(len(parlist))*vol_p(L_curr, L_curr/asp) / (math.Pow(boxL,3)) >= (phi_target/1.3) {
                write_contact_avg(strings.Split(os.Args[1], "_")[0]+"_avg_contacts", contact_m, L_curr, num_unoverlap_steps)
            }
            avg_flag = false
        }
        // taken too many attempts to unoveralp everything
        if num_unoverlap_steps >= STEP_LIMIT - 1 {
            fmt.Println("EXHAUSTED UNOVERLAPPING LOOPS AT L = ", L_curr)
            break vol_loop
        }
        // since vol growth variable, also check if already at request phi!
        if float64(len(parlist))*vol_p(L_curr, L_curr/asp) / (math.Pow(boxL,3)) >= phi_target {
            fmt.Println("Met vol_frac_target = ", phi_target, " early! (breaking)")
            break vol_loop
        }
        // reset map for new L
        contact_m = init_contact_map(len(parlist))
        //loop_len := time.Since(loop_start)
        //fmt.Println(loop_len)
    }

    // once done, grab time again to calc total
    elapsed := time.Since(start)
    fmt.Println("TOOK: ", elapsed)
    // now should have reached desired vol frac, check!
    phi_new := float64(len(parlist))*vol_p(L_curr, L_curr/asp) / (math.Pow(boxL,3))
    fmt.Println("FINAL VOL FRAC: ", phi_new)
    fmt.Println("FINAL L: ", L_curr)
    // for now, just write final positions (have initial)
    if writefile {
        write_position_file(strings.Split(os.Args[1], "_")[0]+"_final_positions"+filenum, parlist, L_curr)
    }
}




// --- SIMULATION SPECFIC FUNCTIONS ---

 


// based on overlap vectors, translate and rotate all particles simultaneously
func update_particles(parlist []Particle, m map[int][][3]float64) []Particle {
    // m is id -> [r_hat_tot, w_tot]
    // rprime = r_old + trans_mag*r_hat_tot
    new_parlist := make([]Particle, len(parlist))
    copy(new_parlist, parlist)  // deepcopy into new list, next update only particles that have overlaps!
    for id, vects := range m {
        // here leveraging that partcle IDs match list order, otherwise have to search
        new_parlist[id].x_pos = periodic_map(new_parlist[id].x_pos + trans_mag*vects[0][0])
        new_parlist[id].y_pos = periodic_map(new_parlist[id].y_pos + trans_mag*vects[0][1])
        new_parlist[id].z_pos = periodic_map(new_parlist[id].z_pos + trans_mag*vects[0][2])
        // to "put back in the box" do we only care about r_com ? -> don't need to fix u_hat ...
       
        new_u := vect_rotate([3]float64{new_parlist[id].ux,new_parlist[id].uy,new_parlist[id].uz},vects[1])
        new_parlist[id].ux = new_u[0]
        new_parlist[id].uy = new_u[1]
        new_parlist[id].uz = new_u[2]
    }
    return new_parlist
}


// remap particle coordinates to -boxl/2 --> boxl/2 box (boxL size centered on origin)
func periodic_map(x float64) float64 {
    if x < -boxL*0.5 {
        x = x + boxL
    }
    if x >= boxL*0.5 {
        x = x - boxL
    }
    return x
}


// check a particle for overlaps, return updated map with entry i -> {w, summed(overlap*n_hat)}
func check_overlaps(particle_i Particle, neighbors []Particle, L_curr float64, m map[int][][3]float64, contact_m map[int]int) (map[int][][3]float64, map[int]int){
    var r_hat_tot [3]float64
    var w_tot [3]float64
    var r_lam [3]float64
    var r_cij [3]float64
    var n_hat [3]float64
    //var id_save int
    contact_count := 0
    for j := 0; j < len(neighbors); j ++ {
        if (particle_i.id != neighbors[j].id) {
            r_lam, r_cij = overlap_vects(particle_i, neighbors[j], L_curr)
            if (dot(r_cij, r_cij) < (L_curr/asp)*L_curr/asp) {
                sep := norm(r_cij)
                // account for contacting!
                contact_count = contact_count + 1
                // add to totals in correct way
                if sep != 0 {
                    // so as long as non-0 vector.. should always be?
                    n_hat = mult(r_cij, 1/sep)
                } else {
                    n_hat = r_cij
                }
                r_hat_tot = add(r_hat_tot, mult(n_hat, L_curr/asp - sep))
                w_tot     = add(w_tot, mult(cross(r_lam,n_hat),L_curr/asp - sep))
   
            }
            // else: this pair isn't contacting, go ahead to next pair
        }
    }
        // ^^... time saved by cell lists/near neigh lists here! goes over ALL j otherwise
    // record number of neighbors this particle had at this unoverlap step
    contact_m[particle_i.id] = contact_m[particle_i.id] + contact_count
    if contact_count > 0 {
        m[particle_i.id] = [][3]float64{r_hat_tot, w_tot}
    }
    return m, contact_m
}


// calc r_lam and r_cij for some j (possible) neighbor of i
// match fortran implementation!
func overlap_vects(particle_i, particle_j Particle, L_curr float64) ([3]float64, [3]float64) {
    // first adjust for min image
    if (particle_i.x_pos - particle_j.x_pos) < -boxL/2 {
        particle_j.x_pos = particle_j.x_pos - boxL
    } else if (particle_i.x_pos - particle_j.x_pos) > boxL/2 {
        particle_j.x_pos = particle_j.x_pos + boxL
    }
    if (particle_i.y_pos - particle_j.y_pos) < -boxL/2 {
        particle_j.y_pos = particle_j.y_pos - boxL
    } else if (particle_i.y_pos - particle_j.y_pos) > boxL/2 {
        particle_j.y_pos = particle_j.y_pos + boxL
    }
    if (particle_i.z_pos - particle_j.z_pos) < -boxL/2 {
        particle_j.z_pos = particle_j.z_pos - boxL
    } else if (particle_i.z_pos - particle_j.z_pos) > boxL/2 {
        particle_j.z_pos = particle_j.z_pos + boxL
    } 

    rcom_A := [3]float64{particle_i.x_pos, particle_i.y_pos, particle_i.z_pos}
    rcom_B := [3]float64{particle_j.x_pos, particle_j.y_pos, particle_j.z_pos}
    // potential speed up is check if rcoms rule out contact
    // then just return vect larger than D by design, so auto ignored, skips all below
    dRcom := add(rcom_A, mult(rcom_B, -1))
    if dot(dRcom,dRcom) > (L_curr*L_curr)*1.001 {
        // return will always be large enough to say: "particles not contacting!"
        return [3]float64{1,0,0}, [3]float64{L_curr, 0, 0}
    }

    ui := [3]float64{particle_i.ux, particle_i.uy, particle_i.uz}
    uj := [3]float64{particle_j.ux, particle_j.uy, particle_j.uz}

    a := dot(ui, dRcom)
    b := dot(uj, dRcom)
    c := dot(ui, uj)
    e := 1.0 - c*c
    var xla float64
    var xniu float64
    xe    := (L_curr - L_curr/asp)/2.0
    small := math.Pow(10.0, -12.0)

    if ( e > small) {
       xla   = ( (b*c) - a) / e
       xniu  = ( b - (a*c)) / e
    } else if ( math.Abs(a) > small ) {
       xla  = math.Copysign(xe, a)
       xniu = b + (xla*c)
       if ( math.Abs(xniu) > xe ) {
            xniu = math.Copysign(xe, xniu)
            goto full_radius
       }
    } else {
       xla  = 0.0
       xniu = 0.0
    }
    if ( math.Abs(xla) <= xe && math.Abs(xniu) <= xe ) {
        goto full_radius
    }
    if ( math.Abs(xla) >= math.Abs(xniu) ) {
       xla  = math.Copysign(xe, xla)
       xniu = b + (xla*c)
    } else {
       xniu = math.Copysign(xe, xniu)
       xla  = (xniu*c) - a
    }
    if ( math.Abs(xniu) > xe ) {
        xniu = math.Copysign(xe, xniu)
    }
    if ( math.Abs(xla) > xe ) {
        xla = math.Copysign(xe, xla)
    }
full_radius:  
    //dd = rr + (xla * xla) + (xniu * xniu) + 2.0*(xla*a - xniu*b - xla*xniu*c)
    //d := math.Sqrt(dd)
    //delta := (L_curr/asp) - d
    v := add(dRcom, mult(ui, xla))
    v = add(v, mult(uj, -1*xniu))
    //vhat := mult(v, 1/d)
    //mi := mult(cross(ui, vhat), xla)

    return mult(ui, xla), v

}


func get_nearest(parlist []Particle, L_curr float64) map[int][]Particle{
    neighbor_map := make(map[int][]Particle)
    // for each particle, look in neaighborhood of L_curr, get only those?
    // so those close enough, [particle id] = append(list, that neighbor)
    // closest enough = ? -> maybe nearest approach (capture those point into circle)
    // search_rad := _*L_curr -> would modify this..
    for i := 0; i < len(parlist); i ++ {
        neighbor_map[i] = []Particle{}
        for j := 0; j < len(parlist); j ++ {
            if i != j {
                _,rclose := overlap_vects(parlist[i],parlist[j], L_curr)
                if dot(rclose,rclose) <= (1.3*1.3)*L_curr*L_curr {
                    neighbor_map[parlist[i].id] = append(neighbor_map[parlist[i].id], parlist[j])
                }
            }
        }
    }
    return neighbor_map
}


// initialize contact_map correctly
func init_contact_map(num_particles int) map[int]int {
    contact_map := make(map[int]int)
    for i := 0; i < num_particles; i++ {
        contact_map[i] = 0
    }
    return contact_map
}


// use rodrigues to rotate partcle u_hat around calculated axis
func vect_rotate(u_hat [3]float64, rot_axis [3]float64) [3]float64 {
    // normalize rot_axis first, if non-zero
    if norm(rot_axis) > 0 {
        rot_axis = mult(rot_axis, 1/norm(rot_axis))
    }
    // pick random angle, in prescribed range
    rand_ang := rot_min + rand.Float64()*(rot_max - rot_min)
    // now apply rodrig formula = v*cos(ang) + cross(axis, v)*sin(ang) + axis*(dot(axis,v))*(1-cos(ang))
    cross_part := mult(cross(rot_axis, u_hat), math.Sin(rand_ang))
    dot_part   := mult(rot_axis, (dot(rot_axis, u_hat))*(1-math.Cos(rand_ang)))
    cos_part   := mult(u_hat, math.Cos(rand_ang))
    new_u_hat  := add(dot_part, add(cos_part, cross_part))
    return new_u_hat
}


// calc curr vol frac
func vol_p(L float64, D float64) float64 {
    return (4.0/3.0)*Pi*math.Pow(D/2, 3) + Pi*math.Pow(D/2, 2)*(L-D)
}


// calc start and end points of rod core line segment
func rod_to_points(rod Particle, L_curr float64) ([3]float64, [3]float64, float64) {
    core_length := L_curr - (L_curr/asp)
    rod_vec := [3]float64{rod.x_pos, rod.y_pos, rod.z_pos}
    u_vec   := [3]float64{rod.ux, rod.uy, rod.uz}
    start := add(rod_vec, mult(mult(u_vec, -1),core_length/2.0) )
    end   := add(rod_vec, mult(u_vec,core_length/2.0) )
    return start, end, core_length
}


// write particle positions and unit vectors at some L
func write_position_file(filename string, parlist []Particle, L_curr float64) {
    position_file := get_write_file(filename)
    position_file.Write([]byte(strconv.FormatFloat(L_curr, 'f', -1, 64)) )
    position_file.Write([]byte("\n"))
    for _, rod := range parlist {
        // write id tab px tab py tab pz tab ux tab uy tab uz
        r_id := strconv.Itoa(rod.id)
        r_x := strconv.FormatFloat(rod.x_pos, 'f', -1, 64)
        r_y := strconv.FormatFloat(rod.y_pos, 'f', -1, 64)
        r_z := strconv.FormatFloat(rod.z_pos, 'f', -1, 64)
        r_ux := strconv.FormatFloat(rod.ux, 'f', -1, 64)
        r_uy := strconv.FormatFloat(rod.uy, 'f', -1, 64)
        r_uz := strconv.FormatFloat(rod.uz, 'f', -1, 64)
        _, err := position_file.Write([]byte(r_id + "\t" + r_x + "\t"+ r_y + "\t"+ r_z + "\t"+ r_ux + "\t"+ r_uy + "\t"+ r_uz + "\n"))
        if err != nil {
            fmt.Println("Couldn't write postion string: ", err)
        }
    }
    position_file.Close()

}


// write L, and each particle avg over unoverlap steps contact #
func write_contact_avg(filename string, contact_m map[int]int, L_curr float64, num_unoverlap_steps int) {
    avg_file := get_write_file(filename)
    // process and write map, currently have map(at some L)[id] -> total contacts over all num_unoverlap_steps
    // write L first
    avg_file.Write([]byte("L = "+strconv.FormatFloat(L_curr, 'f', -1, 64)) )
    avg_file.Write([]byte("\n"))
    for id, contact_count := range contact_m {
        avg := float64(contact_count) / float64(num_unoverlap_steps)
        //fmt.Println("writting ID, #contacts, avg: ", id, "    ", contact_count, "   ", avg)
        avg_str := strconv.Itoa(id) + "\t" + strconv.FormatFloat(avg, 'f', -1, 64)+"\n"
        // write id "tab" avg
        _, err := avg_file.Write([]byte(avg_str))
        if err != nil {
            fmt.Println("Couldn't write avg string: ", err)
        }
    }
    avg_file.Close()
    // close file, above make sure file opened as append style....
}


// read the postions and angles of all aprtcles form some data file
func read_particles(filename string) []Particle {
    var parlist []Particle // store pulled data in this
    curr_directory, err := os.Getwd()
    if err != nil {
        fmt.Println("Couldn't check directory: ", err)
        os.Exit(1)
    }
    data_file, err := os.Open(curr_directory+"/"+filename)
    if err != nil {
        fmt.Println("Couldn't open data file: ", err)
        os.Exit(1)
    }
    defer data_file.Close()
    scanner := bufio.NewScanner(data_file)
    for scanner.Scan() {
        string_data := strings.Fields(scanner.Text())
        id, _ := strconv.Atoi(string_data[0])
        x, _  := strconv.ParseFloat(string_data[1],64)
        y, _  := strconv.ParseFloat(string_data[2],64)
        z, _  := strconv.ParseFloat(string_data[3],64)
        ux, _ := strconv.ParseFloat(string_data[4],64)
        uy, _ := strconv.ParseFloat(string_data[5],64)
        uz, _ := strconv.ParseFloat(string_data[6],64)
        parlist = append(parlist, Particle{id, x,y,z,ux,uy,uz})
    }
    return parlist
}





// -- GENERAL UTILITY FUNCTIONS --





// open and return file via name
func get_write_file(filename string) *os.File {
    curr_directory, err := os.Getwd()
    if err != nil {
        fmt.Println("Couldn't check directory: ", err)
        os.Exit(1)
    }
    file, err := os.OpenFile(curr_directory+"/"+filename, os.O_APPEND|os.O_CREATE|os.O_WRONLY, 0644)
    if err != nil {
        fmt.Println("Couldn't open write file: ", err)
        os.Exit(1)
    }
    return file
}


// calc 3-D determinant
func det(a,b,c [3]float64) float64 {
    // a is 1st column, b 2nd, etc
    return a[0]*(b[1]*c[2] - c[1]*b[2]) - b[0]*(a[1]*c[2] - c[1]*a[2]) + c[0]*(a[1]*b[2] - b[1]*a[2])
}


// calc cross poduct of 3-D vects
func cross(a,b [3]float64) [3]float64 {
    var c [3]float64
    c[0] = a[1]*b[2] - a[2]*b[1]
    c[1] = -1*(a[0]*b[2] - a[2]*b[0])
    c[2] = a[0]*b[1] - a[1]*b[0]
    return c
}


// calc dot product of two 3-D vectors
func dot(a,b [3]float64) float64 {
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]
}


// length of 3-D vect
func norm(a [3]float64) float64 {
    return math.Sqrt(dot(a,a))
}


// add 3-D vectors
func add(a,b [3]float64) [3]float64 {
    c := [3]float64{a[0] + b[0], a[1] + b[1] , a[2] + b[2]}
    return c
}


// mult 3-D vect by scalar
func mult(a [3]float64, m float64) [3]float64 {
    b := [3]float64{m*a[0], m*a[1], m*a[2]}
    return b
}


func memberOf(p_test Particle, parlist []Particle) bool {
    is_member := false
    for _, rod := range parlist {
        if p_test.id == rod.id {
            is_member = true
        }
    }
    return is_member
}