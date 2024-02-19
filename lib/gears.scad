/* Library for Involute Gears, Screws and Racks

This library contains the following modules
- rack(modul, length, height, width, pressure_angle=20, helix_angle=0)
- spur_gear(modul, tooth_number, width, bore, pressure_angle=20, helix_angle=0, optimized=true)
- herringbone_gear(modul, tooth_number, width, bore, pressure_angle=20, helix_angle=0, optimized=true)
- rack_and_pinion (modul, rack_length, gear_teeth, rack_height, gear_bore, width, pressure_angle=20, helix_angle=0, together_built=true, optimized=true)
- ring_gear(modul, tooth_number, width, rim_width, pressure_angle=20, helix_angle=0)
- herringbone_ring_gear(modul, tooth_number, width, rim_width, pressure_angle=20, helix_angle=0)
- planetary_gear(modul, sun_teeth, planet_teeth, number_planets, width, rim_width, bore, pressure_angle=20, helix_angle=0, together_built=true, optimized=true)
- bevel_gear(modul, tooth_number,  partial_cone_angle, tooth_width, bore, pressure_angle=20, helix_angle=0)
- bevel_herringbone_gear(modul, tooth_number, partial_cone_angle, tooth_width, bore, pressure_angle=20, helix_angle=0)
- bevel_gear_pair(modul, gear_teeth, pinion_teeth, axis_angle=90, tooth_width, bore, pressure_angle = 20, helix_angle=0, together_built=true)
- bevel_herringbone_gear_pair(modul, gear_teeth, pinion_teeth, axis_angle=90, tooth_width, bore, pressure_angle = 20, helix_angle=0, together_built=true)
- worm(modul, thread_starts, length, bore, pressure_angle=20, lead_angle=10, together_built=true)
- worm_gear(modul, tooth_number, thread_starts, width, length, worm_bore, gear_bore, pressure_angle=20, lead_angle=0, optimized=true, together_built=true)

Examples of each module are commented out at the end of this file

Author:      Dr Jörg Janssen
Last Verified On:      1. June 2018
Version:    2.2
License:     Creative Commons - Attribution, Non Commercial, Share Alike

Permitted modules according to DIN 780:
0.05 0.06 0.08 0.10 0.12 0.16
0.20 0.25 0.3  0.4  0.5  0.6
0.7  0.8  0.9  1    1.25 1.5
2    2.5  3    4    5    6
8    10   12   16   20   25
32   40   50   60

*/

$fn = 96;
pi = 3.14159;
rad = 57.29578;
clearance = 0.2;
function grad(pressure_angle) = pressure_angle*rad;
function radian(pressure_angle) = pressure_angle/rad;
function polar_to_cartesian(polvect) = [polvect[0]*cos(polvect[1]), polvect[0]*sin(polvect[1])];
function ev(r,rho) = [r/cos(rho), grad(tan(rho)-radian(rho))];
function sphere_ev(theta0,theta) = 1/sin(theta0)*acos(cos(theta)/cos(theta0))-acos(tan(theta0)/tan(theta));
function sphere_to_cartesian(vect) = [vect[0]*sin(vect[1])*cos(vect[2]), vect[0]*sin(vect[1])*sin(vect[2]), vect[0]*cos(vect[1])];
function is_even(number) = (number == floor(number/2)*2) ? 1 : 0;
function ggt(a,b) = a%b == 0 ? b : ggt(b,a%b);
function spiral(a, r0, phi) = a*phi + r0;
module copier(vect, number, distance, winkel){
    for(i = [0:number-1]){
        translate(v=vect*distance*i) rotate(a=i*winkel, v = [0,0,1]) children(0);
    }
}
module rack(modul, length, height, width, pressure_angle = 20, helix_angle = 0) {
    modul=modul*(1-clearance);
    c = modul / 6;
    mx = modul/cos(helix_angle);
    a = 2*mx*tan(pressure_angle)+c*tan(pressure_angle);
    b = pi*mx/2-2*mx*tan(pressure_angle);
    x = width*tan(helix_angle);
    nz = ceil((length+abs(2*x))/(pi*mx));
    translate([-pi*mx*(nz-1)/2-a-b/2,-modul,0]){
        intersection(){
            copier([1,0,0], nz, pi*mx, 0){
                polyhedron(
                    points=[[0,-c,0], [a,2*modul,0], [a+b,2*modul,0], [2*a+b,-c,0], [pi*mx,-c,0], [pi*mx,modul-height,0], [0,modul-height,0],
                        [0+x,-c,width], [a+x,2*modul,width], [a+b+x,2*modul,width], [2*a+b+x,-c,width], [pi*mx+x,-c,width], [pi*mx+x,modul-height,width], [0+x,modul-height,width]],
                    faces=[[6,5,4,3,2,1,0],
                        [1,8,7,0],
                        [9,8,1,2],
                        [10,9,2,3],
                        [11,10,3,4],
                        [12,11,4,5],
                        [13,12,5,6],
                        [7,13,6,0],
                        [7,8,9,10,11,12,13]
                    ]
                );
            };
            translate([abs(x),-height+modul-0.5,-0.5]){
                cube([length,height+modul+1,width+1]);
            }
        };
    };
}
module spur_gear(modul, tooth_number, width, bore, pressure_angle = 20, helix_angle = 0, optimized = true) {
    d = modul * tooth_number;
    r = d / 2;
    alpha_spur = atan(tan(pressure_angle)/cos(helix_angle));
    db = d * cos(alpha_spur);
    rb = db / 2;
    da = (modul <1)? d + modul * 2.2 : d + modul * 2;
    ra = da / 2;
    c =  (tooth_number <3)? 0 : modul/6;
    df = d - 2 * (modul + c);
    rf = df / 2;
    rho_ra = acos(rb/ra);
    rho_r = acos(rb/r);
    phi_r = grad(tan(rho_r)-radian(rho_r));
    gamma = rad*width/(r*tan(90-helix_angle));
    step = rho_ra/16;
    tau = 360/tooth_number;
    r_hole = (2*rf - bore)/8;
    rm = bore/2+2*r_hole;
    z_hole = floor(2*pi*rm/(3*r_hole));
    optimized = (optimized && r >= width*1.5 && d > 2*bore);
    union(){
        rotate([0,0,-phi_r-90*(1-clearance)/tooth_number]){
            linear_extrude(height = width, twist = gamma){
                difference(){
                    union(){
                        tooth_width = (180*(1-clearance))/tooth_number+2*phi_r;
                        circle(rf);
                        for (rot = [0:tau:360]){
                            rotate (rot){
                                polygon(concat([[0,0]], [for (rho = [0:step:rho_ra]) polar_to_cartesian(ev(rb,rho))], [polar_to_cartesian(ev(rb,rho_ra))], [for (rho = [rho_ra:-step:0]) polar_to_cartesian([ev(rb,rho)[0], tooth_width-ev(rb,rho)[1]])]));
                            }
                        }
                    }
                    circle(r = rm+r_hole*1.49);
                }
            }
        }
        if (optimized) {
            linear_extrude(height = width){
                difference(){
                        circle(r = (bore+r_hole)/2);
                        circle(r = bore/2);
                    }
                }
            linear_extrude(height = (width-r_hole/2 < width*2/3) ? width*2/3 : width-r_hole/2){
                difference(){
                    circle(r=rm+r_hole*1.51);
                    union(){
                        circle(r=(bore+r_hole)/2);
                        for (i = [0:1:z_hole]){
                            translate(sphere_to_cartesian([rm,90,i*360/z_hole]))
                                circle(r = r_hole);
                        }
                    }
                }
            }
        }
        else {
            linear_extrude(height = width){
                difference(){
                    circle(r = rm+r_hole*1.51);
                    circle(r = bore/2);
                }
            }
        }
    }
}
module herringbone_gear(modul, tooth_number, width, bore, pressure_angle = 20, helix_angle=0, optimized=true){
    width = width/2;
    d = modul * tooth_number;
    r = d / 2;
    c =  (tooth_number <3)? 0 : modul/6;
    df = d - 2 * (modul + c);
    rf = df / 2;
    r_hole = (2*rf - bore)/8;
    rm = bore/2+2*r_hole;
    z_hole = floor(2*pi*rm/(3*r_hole));
    optimized = (optimized && r >= width*3 && d > 2*bore);
    translate([0,0,width]){
        union(){
            spur_gear(modul, tooth_number, width, 2*(rm+r_hole*1.49), pressure_angle, helix_angle, false);
            mirror([0,0,1]){
                spur_gear(modul, tooth_number, width, 2*(rm+r_hole*1.49), pressure_angle, helix_angle, false);
            }
        }
    }
    if (optimized) {
        linear_extrude(height = width*2){
            difference(){
                    circle(r = (bore+r_hole)/2);
                    circle(r = bore/2);
                }
            }
        linear_extrude(height = (2*width-r_hole/2 < 1.33*width) ? 1.33*width : 2*width-r_hole/2){
            difference(){
                circle(r=rm+r_hole*1.51);
                union(){
                    circle(r=(bore+r_hole)/2);
                    for (i = [0:1:z_hole]){
                        translate(sphere_to_cartesian([rm,90,i*360/z_hole])) circle(r = r_hole);
                    }
                }
            }
        }
    }
    else {
        linear_extrude(height = width*2){
            difference(){
                circle(r = rm+r_hole*1.51);
                circle(r = bore/2);
            }
        }
    }
}
module rack_and_pinion (modul, rack_length, gear_teeth, rack_height, gear_bore, width, pressure_angle=20, helix_angle=0, together_built=true, optimized=true) {
    distance = together_built? modul*gear_teeth/2 : modul*gear_teeth;
    rack(modul, rack_length, rack_height, width, pressure_angle, -helix_angle);
    translate([0,distance,0])
        rotate(a=360/gear_teeth)
            spur_gear (modul, gear_teeth, width, gear_bore, pressure_angle, helix_angle, optimized);
}
module ring_gear(modul, tooth_number, width, rim_width, pressure_angle = 20, helix_angle = 0) {
    ha = (tooth_number >= 20) ? 0.02 * atan((tooth_number/15)/pi) : 0.6;
    d = modul * tooth_number;
    r = d / 2;
    alpha_spur = atan(tan(pressure_angle)/cos(helix_angle));
    db = d * cos(alpha_spur);
    rb = db / 2;
    c = modul / 6;
    da = (modul <1)? d + (modul+c) * 2.2 : d + (modul+c) * 2;
    ra = da / 2;
    df = d - 2 * modul * ha;
    rf = df / 2;
    rho_ra = acos(rb/ra);
    rho_r = acos(rb/r);
    phi_r = grad(tan(rho_r)-radian(rho_r));
    gamma = rad*width/(r*tan(90-helix_angle));
    step = rho_ra/16;
    tau = 360/tooth_number;
    rotate([0,0,-phi_r-90*(1+clearance)/tooth_number]) linear_extrude(height = width, twist = gamma){
        difference(){
            circle(r = ra + rim_width);
            union(){
                tooth_width = (180*(1+clearance))/tooth_number+2*phi_r;
                circle(rf);
                for (rot = [0:tau:360]){
                    rotate (rot) {
                        polygon(concat([[0,0]], [for (rho = [0:step:rho_ra]) polar_to_cartesian(ev(rb,rho))], [polar_to_cartesian(ev(rb,rho_ra))], [for (rho = [rho_ra:-step:0])  polar_to_cartesian([ev(rb,rho)[0], tooth_width-ev(rb,rho)[1]])]));
                    }
                }
            }
        }
    }
}
module herringbone_ring_gear(modul, tooth_number, width, rim_width, pressure_angle = 20, helix_angle = 0) {
    width = width / 2;
    translate([0,0,width]) union(){
        ring_gear(modul, tooth_number, width, rim_width, pressure_angle, helix_angle);
        mirror([0,0,1])
            ring_gear(modul, tooth_number, width, rim_width, pressure_angle, helix_angle);
    }
}
module planetary_gear(modul, sun_teeth, planet_teeth, number_planets, width, rim_width, bore, pressure_angle=20, helix_angle=0, together_built=true, optimized=true){
    d_sun = modul*sun_teeth;
    d_planet = modul*planet_teeth;
    center_distance = modul*(sun_teeth +  planet_teeth) / 2;
    ring_teeth = sun_teeth + 2*planet_teeth;
    d_ring = modul*ring_teeth;
    rotate = is_even(planet_teeth);
    n_max = floor(180/asin(modul*(planet_teeth)/(modul*(sun_teeth +  planet_teeth))));
    rotate([0,0,180/sun_teeth*rotate]){
        herringbone_gear (modul, sun_teeth, width, bore, pressure_angle, -helix_angle, optimized);
    }
    if (together_built){
        if(number_planets==0){
            list = [ for (n=[2 : 1 : n_max]) if ((((ring_teeth+sun_teeth)/n)==floor((ring_teeth+sun_teeth)/n))) n];
            number_planets = list[0];
             center_distance = modul*(sun_teeth + planet_teeth)/2;
            for(n=[0:1:number_planets-1]){
                translate(sphere_to_cartesian([center_distance,90,360/number_planets*n])) rotate([0,0,n*360*d_sun/d_planet]) herringbone_gear (modul, planet_teeth, width, 0, pressure_angle, helix_angle);
            }
       }
       else{
            center_distance = modul*(sun_teeth + planet_teeth)/2;
            for(n=[0:1:number_planets-1]){
                translate(sphere_to_cartesian([center_distance,90,360/number_planets*n])) rotate([0,0,n*360*d_sun/(d_planet)]) herringbone_gear (modul, planet_teeth, width, 0, pressure_angle, helix_angle);
            }
        }
    }
    else{
        planet_distance = ring_teeth*modul/2+rim_width+d_planet;
        for(i=[-(number_planets-1):2:(number_planets-1)]){
            translate([planet_distance, d_planet*i,0]) herringbone_gear (modul, planet_teeth, width, 0, pressure_angle, helix_angle);
        }
    }
    herringbone_ring_gear (modul, ring_teeth, width, rim_width, pressure_angle, helix_angle);
}
module bevel_gear(modul, tooth_number, partial_cone_angle, tooth_width, bore, pressure_angle = 20, helix_angle=0) {
    d_outside = modul * tooth_number;
    r_outside = d_outside / 2;
    rg_outside = r_outside/sin(partial_cone_angle);
    rg_inside = rg_outside - tooth_width;
    r_inside = r_outside*rg_inside/rg_outside;
    alpha_spur = atan(tan(pressure_angle)/cos(helix_angle));
    delta_b = asin(cos(alpha_spur)*sin(partial_cone_angle));
    da_outside = (modul <1)? d_outside + (modul * 2.2) * cos(partial_cone_angle): d_outside + modul * 2 * cos(partial_cone_angle);
    ra_outside = da_outside / 2;
    delta_a = asin(ra_outside/rg_outside);
    c = modul / 6;
    df_outside = d_outside - (modul +c) * 2 * cos(partial_cone_angle);
    rf_outside = df_outside / 2;
    delta_f = asin(rf_outside/rg_outside);
    rkf = rg_outside*sin(delta_f);
    height_f = rg_outside*cos(delta_f);
    height_k = (rg_outside-tooth_width)/cos(partial_cone_angle);
    rk = (rg_outside-tooth_width)/sin(partial_cone_angle);
    rfk = rk*height_k*tan(delta_f)/(rk+height_k*tan(delta_f));
    height_fk = rk*height_k/(height_k*tan(delta_f)+rk);
    phi_r = sphere_ev(delta_b, partial_cone_angle);
    gamma_g = 2*atan(tooth_width*tan(helix_angle)/(2*rg_outside-tooth_width));
    gamma = 2*asin(rg_outside/r_outside*sin(gamma_g/2));
    step = (delta_a - delta_b)/16;
    tau = 360/tooth_number;
    start = (delta_b > delta_f) ? delta_b : delta_f;
    mirrpoint = (180*(1-clearance))/tooth_number+2*phi_r;
    rotate([0,0,phi_r+90*(1-clearance)/tooth_number]){
        translate([0,0,height_f]) rotate(a=[0,180,0]){
            union(){
                translate([0,0,height_f]) rotate(a=[0,180,0]){
                    difference(){
                        linear_extrude(height=height_f-height_fk, scale=rfk/rkf) circle(rkf*1.001);
                        translate([0,0,-1]){
                            cylinder(h = height_f-height_fk+2, r = bore/2);
                        }
                    }
                }
                for (rot = [0:tau:360]){
                    rotate (rot) {
                        union(){
                            if (delta_b > delta_f){
                                flankpoint_under = 1*mirrpoint;
                                flankpoint_over = sphere_ev(delta_f, start);
                                polyhedron(
                                    points = [
                                        sphere_to_cartesian([rg_outside, start*1.001, flankpoint_under]),
                                        sphere_to_cartesian([rg_inside, start*1.001, flankpoint_under+gamma]),
                                        sphere_to_cartesian([rg_inside, start*1.001, mirrpoint-flankpoint_under+gamma]),
                                        sphere_to_cartesian([rg_outside, start*1.001, mirrpoint-flankpoint_under]),
                                        sphere_to_cartesian([rg_outside, delta_f, flankpoint_under]),
                                        sphere_to_cartesian([rg_inside, delta_f, flankpoint_under+gamma]),
                                        sphere_to_cartesian([rg_inside, delta_f, mirrpoint-flankpoint_under+gamma]),
                                        sphere_to_cartesian([rg_outside, delta_f, mirrpoint-flankpoint_under])
                                    ],
                                    faces = [[0,1,2],[0,2,3],[0,4,1],[1,4,5],[1,5,2],[2,5,6],[2,6,3],[3,6,7],[0,3,7],[0,7,4],[4,6,5],[4,7,6]],
                                    convexity =1
                                );
                            }
                            for (delta = [start:step:delta_a-step]){
                                flankpoint_under = sphere_ev(delta_b, delta);
                                flankpoint_over = sphere_ev(delta_b, delta+step);
                                polyhedron(
                                    points = [
                                        sphere_to_cartesian([rg_outside, delta, flankpoint_under]),
                                        sphere_to_cartesian([rg_inside, delta, flankpoint_under+gamma]),
                                        sphere_to_cartesian([rg_inside, delta, mirrpoint-flankpoint_under+gamma]),
                                        sphere_to_cartesian([rg_outside, delta, mirrpoint-flankpoint_under]),
                                        sphere_to_cartesian([rg_outside, delta+step, flankpoint_over]),
                                        sphere_to_cartesian([rg_inside, delta+step, flankpoint_over+gamma]),
                                        sphere_to_cartesian([rg_inside, delta+step, mirrpoint-flankpoint_over+gamma]),
                                        sphere_to_cartesian([rg_outside, delta+step, mirrpoint-flankpoint_over])
                                    ],
                                    faces = [[0,1,2],[0,2,3],[0,4,1],[1,4,5],[1,5,2],[2,5,6],[2,6,3],[3,6,7],[0,3,7],[0,7,4],[4,6,5],[4,7,6]],
                                    convexity =1
                                );
                            }
                        }
                    }
                }
            }
        }
    }
}
module bevel_herringbone_gear(modul, tooth_number, partial_cone_angle, tooth_width, bore, pressure_angle = 20, helix_angle=0){
    tooth_width = tooth_width / 2;
    d_outside = modul * tooth_number;
    r_outside = d_outside / 2;
    rg_outside = r_outside/sin(partial_cone_angle);
    c = modul / 6;
    df_outside = d_outside - (modul +c) * 2 * cos(partial_cone_angle);
    rf_outside = df_outside / 2;
    delta_f = asin(rf_outside/rg_outside);
    height_f = rg_outside*cos(delta_f);
    gamma_g = 2*atan(tooth_width*tan(helix_angle)/(2*rg_outside-tooth_width));
    gamma = 2*asin(rg_outside/r_outside*sin(gamma_g/2));
    height_k = (rg_outside-tooth_width)/cos(partial_cone_angle);
    rk = (rg_outside-tooth_width)/sin(partial_cone_angle);
    rfk = rk*height_k*tan(delta_f)/(rk+height_k*tan(delta_f));
    height_fk = rk*height_k/(height_k*tan(delta_f)+rk);
    modul_inside = modul*(1-tooth_width/rg_outside);
    union(){
        bevel_gear(modul, tooth_number, partial_cone_angle, tooth_width, bore, pressure_angle, helix_angle);
        translate([0,0,height_f-height_fk]) rotate(a=-gamma,v=[0,0,1]) bevel_gear(modul_inside, tooth_number, partial_cone_angle, tooth_width, bore, pressure_angle, -helix_angle);
    }
}
module spiral_bevel_gear(modul, tooth_number, partial_cone_angle, tooth_width, bore, pressure_angle = 20, helix_angle=30){
    steps = 16;
    b = tooth_width / steps;
    d_outside = modul * tooth_number;
    r_outside = d_outside / 2;
    rg_outside = r_outside/sin(partial_cone_angle);
    rg_center = rg_outside-tooth_width/2;
    a=tan(helix_angle)/rg_center;
    union(){
    for(i=[0:1:steps-1]){
        r = rg_outside-i*b;
        helix_angle = a*r;
        modul_r = modul-b*i/rg_outside;
        translate([0,0,b*cos(partial_cone_angle)*i])
            rotate(a=-helix_angle*i,v=[0,0,1])
                bevel_gear(modul_r, tooth_number, partial_cone_angle, b, bore, pressure_angle, helix_angle);
        }
    }
}
module bevel_gear_pair(modul, gear_teeth, pinion_teeth, axis_angle=90, tooth_width, gear_bore, pinion_bore, pressure_angle=20, helix_angle=0, together_built=true){
    r_gear = modul*gear_teeth/2;
    delta_gear = atan(sin(axis_angle)/(pinion_teeth/gear_teeth+cos(axis_angle)));
    delta_pinion = atan(sin(axis_angle)/(gear_teeth/pinion_teeth+cos(axis_angle)));
    rg = r_gear/sin(delta_gear);
    c = modul / 6;
    df_pinion = pi*rg*delta_pinion/90 - 2 * (modul + c);
    rf_pinion = df_pinion / 2;
    delta_f_pinion = rf_pinion/(pi*rg) * 180;
    rkf_pinion = rg*sin(delta_f_pinion);
    height_f_pinion = rg*cos(delta_f_pinion);
    df_gear = pi*rg*delta_gear/90 - 2 * (modul + c);
    rf_gear = df_gear / 2;
    delta_f_gear = rf_gear/(pi*rg) * 180;
    rkf_gear = rg*sin(delta_f_gear);
    height_f_gear = rg*cos(delta_f_gear);
    rotate = is_even(pinion_teeth);
    rotate([0,0,180*(1-clearance)/gear_teeth*rotate])
        bevel_gear(modul, gear_teeth, delta_gear, tooth_width, gear_bore, pressure_angle, helix_angle);
    if (together_built)
        translate([-height_f_pinion*cos(90-axis_angle),0,height_f_gear-height_f_pinion*sin(90-axis_angle)])
            rotate([0,axis_angle,0])
                bevel_gear(modul, pinion_teeth, delta_pinion, tooth_width, pinion_bore, pressure_angle, -helix_angle);
    else
        translate([rkf_pinion*2+modul+rkf_gear,0,0])
            bevel_gear(modul, pinion_teeth, delta_pinion, tooth_width, pinion_bore, pressure_angle, -helix_angle);
 }
module bevel_herringbone_gear_pair(modul, gear_teeth, pinion_teeth, axis_angle=90, tooth_width, gear_bore, pinion_bore, pressure_angle = 20, helix_angle=10, together_built=true){
    r_gear = modul*gear_teeth/2;
    delta_gear = atan(sin(axis_angle)/(pinion_teeth/gear_teeth+cos(axis_angle)));
    delta_pinion = atan(sin(axis_angle)/(gear_teeth/pinion_teeth+cos(axis_angle)));
    rg = r_gear/sin(delta_gear);
    c = modul / 6;
    df_pinion = pi*rg*delta_pinion/90 - 2 * (modul + c);
    rf_pinion = df_pinion / 2;
    delta_f_pinion = rf_pinion/(pi*rg) * 180;
    rkf_pinion = rg*sin(delta_f_pinion);
    height_f_pinion = rg*cos(delta_f_pinion);
    df_gear = pi*rg*delta_gear/90 - 2 * (modul + c);
    rf_gear = df_gear / 2;
    delta_f_gear = rf_gear/(pi*rg) * 180;
    rkf_gear = rg*sin(delta_f_gear);
    height_f_gear = rg*cos(delta_f_gear);
    rotate = is_even(pinion_teeth);
    rotate([0,0,180*(1-clearance)/gear_teeth*rotate])
        bevel_herringbone_gear(modul, gear_teeth, delta_gear, tooth_width, gear_bore, pressure_angle, helix_angle);
    if (together_built)
        translate([-height_f_pinion*cos(90-axis_angle),0,height_f_gear-height_f_pinion*sin(90-axis_angle)])
            rotate([0,axis_angle,0])
                bevel_herringbone_gear(modul, pinion_teeth, delta_pinion, tooth_width, pinion_bore, pressure_angle, -helix_angle);
    else
        translate([rkf_pinion*2+modul+rkf_gear,0,0])
            bevel_herringbone_gear(modul, pinion_teeth, delta_pinion, tooth_width, pinion_bore, pressure_angle, -helix_angle);
}
module worm(modul, thread_starts, length, bore, pressure_angle=20, lead_angle, together_built=true){
    c = modul / 6;
    r = modul*thread_starts/(2*sin(lead_angle));
    rf = r - modul - c;
    a = modul*thread_starts/(90*tan(pressure_angle));
    tau_max = 180/thread_starts*tan(pressure_angle);
    gamma = -rad*length/((rf+modul+c)*tan(lead_angle));
    step = tau_max/16;
    if (together_built) {
        rotate([0,0,tau_max]){
            linear_extrude(height = length, center = false, convexity = 10, twist = gamma){
                difference(){
                    union(){
                        for(i=[0:1:thread_starts-1]){
                            polygon(
                                concat(
                                    [[0,0]],
                                    [for (tau = [0:step:tau_max])
                                        polar_to_cartesian([spiral(a, rf, tau), tau+i*(360/thread_starts)])],
                                    [for (tau = [tau_max:step:180/thread_starts])
                                        polar_to_cartesian([spiral(a, rf, tau_max), tau+i*(360/thread_starts)])],
                                    [for (tau = [180/thread_starts:step:(180/thread_starts+tau_max)])
                                        polar_to_cartesian([spiral(a, rf, 180/thread_starts+tau_max-tau), tau+i*(360/thread_starts)])]
                                )
                            );
                        }
                        circle(rf);
                    }
                    circle(bore/2);
                }
            }
        }
    }
    else {
        difference(){
            union(){
                translate([1,r*1.5,0]){
                    rotate([90,0,90])
                        worm(modul, thread_starts, length, bore, pressure_angle, lead_angle, together_built=true);
                }
                translate([length+1,-r*1.5,0]){
                    rotate([90,0,-90])
                        worm(modul, thread_starts, length, bore, pressure_angle, lead_angle, together_built=true);
                    }
                }
            translate([length/2+1,0,-(r+modul+1)/2]){
                    cube([length+2,3*r+2*(r+modul+1),r+modul+1], center = true);
                }
        }
    }
}
module worm_gear(modul, tooth_number, thread_starts, width, length, worm_bore, gear_bore, pressure_angle=20, lead_angle, optimized=true, together_built=true, show_spur=1, show_worm=1){
    c = modul / 6;
    r_worm = modul*thread_starts/(2*sin(lead_angle));
    r_gear = modul*tooth_number/2;
    rf_worm = r_worm - modul - c;
    gamma = -90*width*sin(lead_angle)/(pi*r_gear);
    tooth_distance = modul*pi/cos(lead_angle);
    x = is_even(thread_starts)? 0.5 : 1;
    if (together_built) {
        if(show_worm)
        translate([r_worm,(ceil(length/(2*tooth_distance))-x)*tooth_distance,0])
            rotate([90,180/thread_starts,0])
                worm(modul, thread_starts, length, worm_bore, pressure_angle, lead_angle, together_built);
        if(show_spur)
        translate([-r_gear,0,-width/2])
            rotate([0,0,gamma])
                spur_gear (modul, tooth_number, width, gear_bore, pressure_angle, -lead_angle, optimized);
    }
    else {
        if(show_worm)
        worm(modul, thread_starts, length, worm_bore, pressure_angle, lead_angle, together_built);
        if(show_spur)
        translate([-2*r_gear,0,0])
            spur_gear (modul, tooth_number, width, gear_bore, pressure_angle, -lead_angle, optimized);
    }
}

// rack(modul=1, length=30, height=5, width=5, pressure_angle=20, helix_angle=20);

// spur_gear (modul=1, tooth_number=30, width=5, bore=4, pressure_angle=20, helix_angle=20, optimized=true);

// herringbone_gear (modul=1, tooth_number=30, width=5, bore=4, pressure_angle=20, helix_angle=30, optimized=true);

// rack_and_pinion (modul=1, rack_length=50, gear_teeth=30, rack_height=4, gear_bore=4, width=5, pressure_angle=20, helix_angle=0, together_built=true, optimized=true);

// ring_gear (modul=1, tooth_number=30, width=5, rim_width=3, pressure_angle=20, helix_angle=20);

// herringbone_ring_gear (modul=1, tooth_number=30, width=5, rim_width=3, pressure_angle=20, helix_angle=30);

// planetary_gear(modul=1.2, sun_teeth=9, planet_teeth=6, number_planets=5, width=19.5, rim_width=1, bore=7.1, pressure_angle=20, helix_angle=0, together_built=true, optimized=false);

// bevel_gear(modul=1, tooth_number=30,  partial_cone_angle=45, tooth_width=5, bore=4, pressure_angle=20, helix_angle=20);

// bevel_herringbone_gear(modul=1, tooth_number=30, partial_cone_angle=45, tooth_width=5, bore=4, pressure_angle=20, helix_angle=30);

// bevel_gear_pair(modul=1, gear_teeth=30, pinion_teeth=11, axis_angle=100, tooth_width=5, bore=4, pressure_angle = 20, helix_angle=20, together_built=true);

// bevel_herringbone_gear_pair(modul=1, gear_teeth=30, pinion_teeth=11, axis_angle=100, tooth_width=5, bore=4, pressure_angle = 20, helix_angle=30, together_built=true);

// worm(modul=1, thread_starts=2, length=15, bore=4, pressure_angle=20, lead_angle=10, together_built=true);

// worm_gear(modul=1, tooth_number=30, thread_starts=2, width=8, length=20, worm_bore=4, gear_bore=4, pressure_angle=20, lead_angle=10, optimized=1, together_built=1, show_spur=1, show_worm=1);
